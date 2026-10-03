"""Deterministic crisprprimer computation, distributed by BioV.

Authors: Xiangchao Gan <gan@njau.edu.cn> & Ziya Tang <tang@stu.njau.edu.cn>
License: MIT
"""

import json
import logging
import re
from collections.abc import Sequence
from functools import partial
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO, Literal, cast

import fsspec
import pandas as pd
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction, nt_search

from biov import BioDataFrame, read_fasta, read_gff3, settings
from biov.file_sources import copy_stream
from biov.io._preprocess import preprocessing
from biov.software import run_software

from .id_converter import guess_id_system
from .nuclease import Nuclease, Nucleases, SpCas9
from .score import cfd_score, pair_score

logger = logging.getLogger("crisprprimer")

BASE_PATH = Path(__file__).parent
RESTRICTION_ENZYMES = json.loads((BASE_PATH / "restriction_enzymes.json").read_text())
GENE_FAMILIES = {
    "EL5": {
        "Os02g0559800",
        "Os02g0560200",
        "Os02g0560600",
        "Os02g0561000",
        "Os02g0561400",
        "Os02g0561800",
        "LOC_Os02g35329",
        "LOC_Os02g35365",
        "LOC_Os02g35347",
    },
}
_BLAT_COLUMNS = [  # Same as output from https://genome.ucsc.edu/cgi-bin/hgBlat
    "matches",
    "misMatches",
    "repMatches",
    "nCount",
    "qNumInsert",
    "qBaseInsert",
    "tNumInsert",
    "tBaseInsert",
    "strand",
    "qName",
    "qSize",
    "qStart",
    "qEnd",
    "tName",
    "tSize",
    "tStart",
    "tEnd",
    "blockCount",
    "blockSizes",
    "qStarts",
    "tStarts",
]


def _parse_region(reg: str) -> tuple[str, int, int] | None:
    m = re.match(r"(\w+):(\d+)(?:-|\.\.)(\d+)", reg)
    if m:
        chrom, start, end = m.groups()
        return chrom, int(start) - 1, int(end)
    return None


def get_pairs(
    spacers: BioDataFrame,
    window_size: tuple[int, int],
    nuclease: Nuclease,
) -> BioDataFrame:
    """Enumerate the existing opposite-strand pairs within the window.

    Returns:
        Pairs with coordinates and existing calculated annotations.
    """
    df_plus = spacers[spacers["strand"] == "+"]
    df_minus = spacers[spacers["strand"] == "-"]
    records = []
    for (_, plus), (_, minus) in product(df_plus.iterrows(), df_minus.iterrows()):
        if nuclease.pam_side == "3prime":
            start = minus["start"]
            end = plus["end"]
            margin = abs(min(nuclease.cut_sites))
            spacer1 = minus
            spacer2 = plus
        else:
            start = plus["start"]
            end = minus["end"]
            margin = max(nuclease.cut_sites)
            spacer1 = plus
            spacer2 = minus
        if window_size[0] < (ws := end - start - 2 * margin) < window_size[1]:
            records.append(
                {
                    "start": start,
                    "end": end,
                    "spacer1": spacer1["ID"],
                    "spacer2": spacer2["ID"],
                    "spacer1_seq": spacer1["spacer"],
                    "spacer2_seq": spacer2["spacer"],
                    "spacer1_with_linker": "gcggtctcaggcg"
                    + spacer1["spacer"]
                    + "gttttagagctagaaatagcaag",
                    "spacer2_with_linker": "ttggtctctaaac"
                    + reverse_complement(spacer2["spacer"])
                    + "cacacaagcgacagc",
                    "window_size": ws,
                }
            )
    pairs = BioDataFrame(records)
    pairs["seqid"] = spacers["seqid"].iloc[0]
    return pairs


def get_extra_bases(nuclease: Nuclease) -> int:
    """Compute the original search margin outside the target CDS.

    Returns:
        Number of additional bases.
    """
    if nuclease.pam_side == "3prime":
        return 2 * (abs(min(nuclease.cut_sites)) + len(nuclease.pam)) + 3
    else:
        return 2 * max(nuclease.cut_sites) + 3


def _exclude_gc(spacers: BioDataFrame) -> BioDataFrame:
    _gc_fraction = spacers["spacer"].apply(gc_fraction)
    return cast(BioDataFrame, spacers[(_gc_fraction >= 0.2) & (_gc_fraction <= 0.8)])


def _exclude_polyT(
    spacers: BioDataFrame,
) -> BioDataFrame:
    return cast(BioDataFrame, spacers[~spacers["spacer"].str.contains("TTTT")])


def _exclude_restriction_enzymes(
    spacers: BioDataFrame,
    restriction_enzymes: Sequence[str],
    flanking: tuple[tuple[str, str], tuple[str, str]],
    strand: Literal["+", "-"] = "+",
) -> BioDataFrame:
    for enzyme in restriction_enzymes:
        flanking5 = spacers["strand"].apply(lambda x: flanking[int(x == strand)][0])
        flanking3 = spacers["strand"].apply(lambda x: flanking[int(x == strand)][1])
        spacer_with_flanking = flanking5 + spacers["spacer"] + flanking3
        enzyme_seq = RESTRICTION_ENZYMES[enzyme]
        rc_enzyme_seq = reverse_complement(enzyme_seq)
        spacers = cast(
            BioDataFrame,
            spacers[
                (
                    (~spacer_with_flanking.str.contains(enzyme_seq).astype(bool))
                    & (~spacer_with_flanking.str.contains(rc_enzyme_seq).astype(bool))
                )
            ],
        )
    return spacers


def _blat(
    ref: str,
    seq: str,
    cache_dir: Path | None = None,
) -> pd.DataFrame:
    """Run the existing BLAT comparison through BioV's native software executor.

    Returns:
        Native PSL columns, optionally reused from the caller's reference cache.

    Raises:
        ValueError: If native execution points at a different host.
    """
    if settings.execution_host is not None:
        raise ValueError(
            "crisprprimer BLAT uses local files; run the complete Python script "
            "on the execution host instead of configuring remote native execution"
        )
    cached = pd.DataFrame(columns=_BLAT_COLUMNS)
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        subdir = cache_dir / f"5-mer={seq[:5]}"
        if subdir.exists():
            cached = pd.read_parquet(subdir)
            if seq in cached["qName"].to_numpy():
                return cast(pd.DataFrame, cached[cached["qName"] == seq].copy())
    with TemporaryDirectory(prefix="crisprprimer-blat-") as directory:
        work = Path(directory)
        query = work / "query.fa"
        query.write_text(f">{seq}\n{seq}\n")
        output = work / "hits.psl"
        reference = ref
        if "://" in ref or ref.endswith(".gz"):
            reference = str(work / "reference.fa")
            source, options = preprocessing(ref)
            with (
                fsspec.open(source, "rb", compression="infer", **options) as reader,
                Path(reference).open("wb") as writer,
            ):
                copy_stream(cast(BinaryIO, reader), writer)
        run_software(
            "blat",
            (
                "-noHead",
                "-minMatch=0",
                "-minScore=18",
                "-stepSize=5",
                "-fine",
                reference,
                str(query),
                str(output),
            ),
            cwd=Path.cwd(),
        ).check_returncode()
        result = pd.read_table(output, names=_BLAT_COLUMNS)
    if cache_dir is not None and not result.empty:
        combined = pd.concat([cached, result], ignore_index=True)
        combined["5-mer"] = seq[:5]
        combined.to_parquet(
            cache_dir,
            partition_cols=["5-mer"],
            existing_data_behavior="delete_matching",
        )
    return result


def get_pam(
    row: pd.Series, seqs: dict[str, SeqRecord], nuclease: Nuclease = SpCas9
) -> str:
    """Extract the PAM beside a native PSL alignment.

    Returns:
        Uppercase oriented PAM sequence.

    Raises:
        ValueError: If the referenced record has no sequence.
    """
    tName = cast(str, row["tName"])
    tStart = row["tStart"]
    tEnd = row["tEnd"]
    strand = row["strand"]
    sequence = seqs[tName].seq
    if sequence is None:
        raise ValueError(f"Reference record {tName} has no sequence")
    if nuclease.pam_side == "3prime":
        if strand == "+":
            return str(
                sequence[
                    tEnd + nuclease.spacer_range[1] : tEnd
                    + nuclease.spacer_range[1]
                    + len(nuclease.pam)
                ]
            ).upper()
        else:
            return str(
                sequence[
                    tStart - nuclease.spacer_range[1] - len(nuclease.pam) : tStart
                    - nuclease.spacer_range[1]
                ].reverse_complement()
            ).upper()
    else:
        if strand == "+":
            return str(
                sequence[
                    tStart - nuclease.spacer_range[0] : tStart
                    - nuclease.spacer_range[0]
                    + len(nuclease.pam)
                ]
            ).upper()
        else:
            return str(
                sequence[
                    tEnd + nuclease.spacer_range[0] - len(nuclease.pam) : tEnd
                    + nuclease.spacer_range[0]
                ].reverse_complement()
            ).upper()


def get_protospacer(row: pd.Series, seqs: dict[str, SeqRecord]) -> str:
    """Extract the existing aligned protospacer context.

    Returns:
        Uppercase oriented protospacer sequence.

    Raises:
        ValueError: If the referenced record has no sequence.
    """
    tName = cast(str, row["tName"])
    sequence = seqs[tName].seq
    if sequence is None:
        raise ValueError(f"Reference record {tName} has no sequence")
    strand = row["strand"]
    if row["blockCount"] == 2:
        blockSizes = [
            int(size) for size in cast(str, row["blockSizes"]).split(",")[:-1]
        ]
        qStarts = [int(start) for start in cast(str, row["qStarts"]).split(",")[:-1]]
        tStarts = [int(start) for start in cast(str, row["tStarts"]).split(",")[:-1]]
        # get biggest block's qStart & tStart
        if blockSizes[0] > blockSizes[1]:
            qStart = qStarts[0]
            qEnd = qStart + blockSizes[0]
            tStart = tStarts[0]
            tEnd = tStart + blockSizes[0]
        else:
            qStart = qStarts[1]
            qEnd = qStarts[1] + blockSizes[1]
            tStart = tStarts[1]
            tEnd = tStart + blockSizes[1]
    else:
        tStart = row["tStart"]
        tEnd = row["tEnd"]

        qStart = row["qStart"]
        qEnd = row["qEnd"]
    if strand == "+":
        return str(sequence[tStart - qStart : tEnd - qEnd + 20]).upper()
    else:
        return str(
            sequence[tStart - qStart : tEnd - qEnd + 20].reverse_complement()
        ).upper()


def _off_target(
    ref: str,
    seqs: dict[str, SeqRecord],
    gff: BioDataFrame,
    locus_id: str,
    query: str,
    blat_cache_dir: Path | None = None,
    nuclease: Nuclease = SpCas9,
):
    for genes in GENE_FAMILIES.values():
        if locus_id in genes:
            genes_df = gff[gff["ID"].isin(list(genes))]
            break
    else:
        genes_df = gff.query("ID == @locus_id").copy()
    # query string not in all members of one gene family
    qs = " and ".join(
        [
            f"(tName != '{gene['seqid']}' or tStart >= {gene['end']} or tEnd <= {gene['start']})"
            for _, gene in genes_df.iterrows()
        ]
        + ["blockCount <= 2"]
    )
    blat_result = _blat(ref, query, cache_dir=blat_cache_dir)
    candidates = blat_result.query(qs).copy()
    if len(candidates) == 0:
        return False
    candidates["pam"] = candidates.apply(
        partial(get_pam, seqs=seqs, nuclease=nuclease), axis=1
    )
    candidates["protospacer"] = candidates.apply(
        partial(get_protospacer, seqs=seqs), axis=1
    )
    candidates = candidates[
        ~candidates["protospacer"].str.contains("N")
        & ~candidates["pam"].str.contains("N")
    ]
    if str(nuclease.pam) == "NGG":
        candidates["score"] = candidates.apply(
            lambda row: cfd_score(row["qName"], row["protospacer"], row["pam"][1:]),
            axis=1,
        )
        return (candidates["score"] > 0.2).any()
    else:
        return (
            len(
                candidates.query(
                    f"matches == {nuclease.spacer_length} and blockCount == 1"
                )
            )
            > 1
        )


def _exclude_off_targets(
    spacers: BioDataFrame,
    ref: str,
    seqs: dict[str, SeqRecord],
    gff: BioDataFrame,
    blat_cache_dir: Path | None = None,
    nuclease: Nuclease = SpCas9,
) -> BioDataFrame:
    return cast(
        BioDataFrame,
        spacers[
            ~spacers.apply(
                lambda row: _off_target(
                    ref,
                    seqs,
                    gff,
                    row["locus_id"],
                    row["spacer"],
                    blat_cache_dir,
                    nuclease,
                ),
                axis=1,
            )
        ],
    )


def _on_target(
    ind: str,
    ind_seqs: dict[str, SeqRecord],
    query: str,
    ind_blat_cache_dir: Path | None = None,
    nuclease: Nuclease = SpCas9,
):
    result = _blat(ind, query, cache_dir=ind_blat_cache_dir).query(
        f"matches == {nuclease.spacer_length} and blockCount == 1"
    )
    if len(result) == 0:
        return False
    pam = get_pam(result.iloc[0], ind_seqs, nuclease)
    _, *pos = nt_search(pam, str(nuclease.pam))
    return len(result) == 1 and len(pos) == 1 and pos[0] == 0


def _on_target_in_ind(
    ind: str,
    ind_seqs: dict[str, SeqRecord],
    spacers: BioDataFrame,
    ind_blat_cache_dir: Path | None = None,
    nuclease: Nuclease = SpCas9,
) -> BioDataFrame:
    return cast(
        BioDataFrame,
        spacers[
            spacers.apply(
                lambda row: _on_target(
                    ind, ind_seqs, row["spacer"], ind_blat_cache_dir, nuclease
                ),
                axis=1,
            )
        ],
    )


def _crispr_for_one_region(
    region: str,
    annotation: BioDataFrame,
    system: Nucleases,
    window_size: tuple[int, int],
    seqs: dict[str, SeqRecord],
    ref: str,
    restriction_enzymes: Sequence[str],
    flanking: tuple[tuple[str, str], tuple[str, str]],
    strand: Literal["+", "-"] = "+",
    blocks_dir: Path | None = None,
    ind: str | None = None,
    ind_seqs: dict[str, SeqRecord] | None = None,
    ref_blat_cache_dir: Path | None = None,
    ind_blat_cache_dir: Path | None = None,
    processing: Literal["pre", "post"] = "pre",
):
    parsed = _parse_region(region)
    if parsed is not None:
        seqid, start, end = parsed
        cds = annotation.query(
            "type == 'CDS' and seqid == @seqid and start >= @start and end <= @end"
        )
    else:
        id_system = guess_id_system(region)
        if id_system is None:
            raise ValueError(f"Your input id {region} is not RAP id or MSU id.")
        if id_system != "MSU":
            raise ValueError(
                "Use an MSU identifier or genomic region; RAP annotation lookup is not implemented"
            )
        cds = annotation.query("type == 'CDS' and ID.str.startswith(@region)")
    if len(cds) == 0:
        logger.warning(f"Your input id {region} could not find any CDS.")
        return
    seqid = cast(str, cds["seqid"].iloc[0])
    extra_bases = get_extra_bases(system.model)
    size = window_size[1] + extra_bases
    start = max(0, cds["start"].min() - size)
    end = min(cds["end"].max() + size, len(seqs[seqid]))
    sequence = seqs[seqid].seq
    if sequence is None:
        raise ValueError(f"Reference record {seqid} has no sequence")
    spacers = system.model.find_spacers_on(sequence[start:end].upper())
    if spacers.empty:
        return
    spacers["nuclease"] = system.name
    spacers["seqid"] = seqid
    spacers["locus_id"] = region
    spacers[["start", "end"]] += start
    spacers["ID"] = (
        seqid
        + ":"
        + (spacers["start"] + 1).astype(str)
        + ".."
        + spacers["end"].astype(str)
        + " ("
        + spacers["strand"]
        + "strand)"
    )
    spacers = _exclude_gc(spacers)
    if spacers.empty:
        return
    spacers = _exclude_polyT(spacers)
    if spacers.empty:
        return
    spacers = _exclude_restriction_enzymes(
        spacers, restriction_enzymes, flanking, strand
    )
    if spacers.empty:
        return
    if processing == "pre":
        spacers = _exclude_off_targets(
            spacers, ref, seqs, annotation, ref_blat_cache_dir, system.model
        )
        if spacers.empty:
            return
        if ind_seqs is not None and ind is not None:
            spacers = _on_target_in_ind(
                ind, ind_seqs, spacers, ind_blat_cache_dir, system.model
            )
            if spacers.empty:
                return
    if blocks_dir is not None:
        d = blocks_dir / "spacers" / f"seqid={seqid}"
        d.mkdir(parents=True, exist_ok=True)
        spacers.to_parquet(d / f"{region}.parquet")
    if system.name != "SpCas9":
        return
    pairs = get_pairs(spacers, window_size, system.model)
    if pairs.empty:
        return
    pairs["locus_id"] = region
    pairs["nuclease"] = system.name
    pairs["overlap_with_all_CDS"] = pairs.apply(
        lambda row: ((cds["start"] <= row["end"]) & (cds["end"] >= row["start"])).all(),
        axis=1,
    )
    pairs["score"] = pairs.apply(
        partial(pair_score, cds=cds, dist=abs(system.model.cut_sites[1])),
        axis=1,
    )
    if processing == "post":
        for i, row in pairs.iterrows():
            if not _off_target(
                ref,
                seqs,
                annotation,
                cast(str, row["locus_id"]),
                cast(str, row["spacer1_seq"]),
                ref_blat_cache_dir,
            ) and not _off_target(
                ref,
                seqs,
                annotation,
                cast(str, row["locus_id"]),
                cast(str, row["spacer2_seq"]),
                ref_blat_cache_dir,
            ):
                if ind_seqs is not None and ind is not None:
                    if _on_target(
                        ind, ind_seqs, cast(str, row["spacer1_seq"]), ind_blat_cache_dir
                    ) and _on_target(
                        ind, ind_seqs, cast(str, row["spacer2_seq"]), ind_blat_cache_dir
                    ):
                        pairs = pairs.iloc[[i]]
                        break
                else:
                    pairs = pairs.iloc[[i]]
                    break
        else:
            return
    if blocks_dir is not None:
        d = blocks_dir / "pairs" / f"seqid={seqid}"
        d.mkdir(parents=True, exist_ok=True)
        pairs.to_parquet(d / f"{region}.parquet")


def crisprprimer(
    ref: Literal[
        "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz",
        "https://rice.uga.edu/osa1r7_download/osa1_r7.asm.fa.gz",
        "https://rice.uga.edu/osa1r7_download/osa1_r7.asm.chrs.fa.gz",
    ] = "https://rice.uga.edu/osa1r7_download/osa1_r7.asm.chrs.fa.gz",
    gff: Literal[
        "https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.gff3.gz",
        (
            "tar://IRGSP-1.0_representative/transcripts.gff::"
            "https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/"
            "IRGSP-1.0_representative_2025-03-19.tar.gz"
        ),
    ] = "https://rice.uga.edu/osa1r7_download/osa1_r7.all_models.gff3.gz",
    regions: list[str] | None = None,
    window_size: tuple[int, int] = (37, 70),
    system: Nucleases = Nucleases.SpCas9,
    flanking: tuple[tuple[str, str], tuple[str, str]] = (
        ("GGCG", "GTTT"),
        ("TGTG", "GTTT"),
    ),
    restriction_enzymes: Sequence[str] = (
        "BbsI",
        "BsaI",
    ),
    blocks_dir: Path | None = None,
    strand: Literal["+", "-"] = "+",
    ind: str | None = None,
    ref_blat_cache_dir: Path | None = None,
    ind_blat_cache_dir: Path | None = None,
    processing: Literal["pre", "post"] = "pre",
):
    """Design CRISPR primers.

    Parameters
    ----------
    ref : str
        The reference genome file.
    gff : str
        The GFF file.
    regions : list[str] | None, optional
        The region to search for CRISPR primers, Locus_id (in GFF's attributes) is also accepted. If None provided, all gene in GFF will process.
        Format: "chr02:21250083..21251343" or "LOC_Os02g35329". RAP annotation lookup is not implemented.
    window_size : tuple[int, int], optional
        The window size to search for CRISPR primers.
        Default: (37, 70).
    system : System, optional
        The system of the sequence. Default: System.CAS9.
    flanking : Sequence[tuple[str, str]], optional
        The flanking sequences, each tuple is a pair of flanking sequences. Default: (("GGCG", "GTTT"), ("TGTG", "GTTT")). It makes pairs of spacer different among strands.
    restriction_enzymes : Sequence[str], optional
        The restriction enzymes to exclude. Default: ("BbsI", "BsaI").
    blocks_dir : Path | None, optional
        The directory to store the intermediate results. If None, the results will not be stored.
    strand : Literal["+", "-"], optional
        The strand to search for CRISPR primers. Default: "+".
    ind : str | None, optional
        The individual genome sequence to search for CRISPR primers. Default: None.
    ref_blat_cache_dir : Path | None, optional
        The directory to store the BLAT results for reference genome. If None, the results will not be stored.
    ind_blat_cache_dir : Path | None, optional
        The directory to store the BLAT results for cultivar. If None, the results will not be stored.
    processing : Literal["pre", "post"], optional
        Whether to process the spacers before or after pairing, post-processing would be much faster but would not filter out all off-targets. Default: "pre".

    Notes:
    -----
    Double-stranded
    ```
    5'                                                          ┌────> spacer2 ─────┐ 3'
    CCNNNN NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN NNNNGG
           └────────────────────────────── window ──────────────────────────────┘^└┬┘└┬┘
                                                                   cutting point─┤ │ PAM
    3' ┌───── spacer1 <────┐                                                     v └dist
    GGNNNN NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN NNNNCC
    ```

    Cas12a Double-stranded
    ```
    5'  ┌───── 18 bp ────┐                                   ┌─────── 23 bp ───────┐  3'
    TTTVNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNNNNNBAAA
    └┬─┘└── spacer1 (21 bp) ─┘  └───────── window ─────────┘
    PAM
    3'  ┌─────── 23 bp ───────┐                               ┌────── spacer2 ─────┐
    AAABNNNNNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNNNNVTTT
    ```
    """
    seqs: dict[str, SeqRecord] = read_fasta(ref)

    # build a map from transcript id to locus id
    annotation = read_gff3(gff)
    if regions is None:
        regions = cast(
            list[str],
            annotation.query("type == 'gene' and seqid in @seqs")["ID"].to_list(),
        )
    if ind is not None:
        ind_seqs = read_fasta(ind)
    else:
        ind_seqs = None
    for region in regions:
        _crispr_for_one_region(
            region,
            annotation,
            system,
            window_size,
            seqs,
            ref,
            restriction_enzymes,
            flanking,
            strand,
            blocks_dir,
            ind,
            ind_seqs,
            ref_blat_cache_dir,
            ind_blat_cache_dir,
            processing=processing,
        )
