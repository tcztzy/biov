"""Fast sequence search command line tool."""

import platform
import re
import shutil
import subprocess  # noqa: S404
from contextlib import ExitStack
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Annotated, Any

import fsspec
import pandas as pd
from Bio.Data.IUPACData import (
    protein_letters,
    unambiguous_dna_letters,
    unambiguous_rna_letters,
)
from click import Choice
from typeguard import typechecked
from typer import Argument, Option, Typer

from ..config import settings

_FSSPEC_URL_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9+\-+.]*(::[A-Za-z0-9+\-+.]+)*://")

_SUB_RANGE_PATTERN = re.compile(
    r"(?P<path>.+\.2bit|\.nib)(?::(?P<seqid>[a-zA-Z]\w+))?:(?P<start>\d+)-(?P<end>\d+)$"
)

_PSL_HEAD = """psLayout version 3

match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
---------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
_RSYNC_SUPPORTED = platform.system() == "Darwin" or (
    platform.system() == "Linux" and platform.machine() == "x86_64"
)
_SUBDIR = (
    f"{'macOSX' if platform.system() == 'Darwin' else 'linux'}.{platform.machine()}"
)
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

blat_app = Typer()


def _handle_fsspec_path(path: str, stack: ExitStack) -> str:
    """Handle fsspec path.

    Args:
        path: the path need to handle
        stack: used for create temporary file

    Returns:
        original path if not fsspec else temporary file name.
    """
    if _FSSPEC_URL_PATTERN.match(path):
        *protocols, path = path.split("::")
        if len(protocols) == 0 or protocols[-1] != "filecache":
            protocols = [*protocols, "filecache"]
        db = stack.enter_context(NamedTemporaryFile())
        with fsspec.open("::".join([*protocols, path]), compression="infer") as x:
            shutil.copyfileobj(x, db)  # pyright: ignore
        path = db.name
    return path


def _handle_seq(path: str, stack: ExitStack) -> str:
    """Handle sequence as query.

    Args:
        path: the path need to handle
        stack: used for create temporary file

    Returns:
        original path if not seq else temporary file name.
    """
    if (
        _SUB_RANGE_PATTERN.match(path) is None
        and (
            len(path) >= 256  # Maximum file name length
            or not Path(path).exists()
        )
        and re.match(
            rf"^([{unambiguous_dna_letters}]+|[{unambiguous_rna_letters}]+|[{protein_letters}]+)$",
            path,
            flags=re.IGNORECASE,
        )
    ):
        f = stack.enter_context(NamedTemporaryFile("w+t", suffix=".fa"))
        name = path if len(path) <= 23 else f"{path[:11]}_{len(path)}"
        f.write(f">{name}\n{path}")
        f.flush()
        path = f.name
    return path


@blat_app.command(
    epilog="To filter PSL files to the best hits (e.g.  minimum ID > 90% or 'only best match'), you can use the commands pslReps, pslCDnaFilter or pslUniq."
)
def blat_cli(
    database: Annotated[
        str,
        Argument(
            help="a .fa, .nib or .2bit file, or a list of these files with one file name per line.",
            show_default=False,
        ),
    ],
    query: Annotated[
        list[str],
        Argument(
            help="Same as database but additionally support directly providing sequence. Last one is output for backward compatibility.",
            show_default=False,
        ),
    ],
    *,
    t: Annotated[
        str,
        Option("-t", click_type=Choice(["dna", "prot", "dnax"]), help="Database type."),
    ] = "dna",
    q: Annotated[
        str,
        Option(
            "-q",
            click_type=Choice(["dna", "rna", "prot", "dnax", "rnax"]),
            help="Query type.",
        ),
    ] = "dna",
    prot: Annotated[
        bool, Option("-prot", help="Synonymous with -t=prot -q=prot.")
    ] = False,
    ooc: Annotated[
        str | None,
        Option(
            "-ooc",
            help="Use overused tile file N.ooc. N should correspond to the tileSize",
        ),
    ] = None,
    tileSize: Annotated[
        int | None,
        Option(
            "-tileSize",
            help="Sets the size of match that triggers an alignment. Usually between 8 and 12.",
            show_default="11 for DNA and 5 for protein",
        ),
    ] = None,
    stepSize: Annotated[
        int | None,
        Option("-stepSize", help="Spacing between tiles.", show_default="tileSize"),
    ] = None,
    oneOff: Annotated[
        str,
        Option(
            "-oneOff",
            click_type=Choice(["0", "1"]),
            help="If set to 1, this allows one mismatch in tile and still triggers an alignment.",
        ),
    ] = "0",
    minMatch: Annotated[
        int | None,
        Option(
            "-minMatch",
            help="Sets the number of tile matches.  Usually set from 2 to 4.",
            show_default="2 for nucleotide, 1 for protein",
        ),
    ] = None,
    minScore: Annotated[
        int,
        Option(
            "-minScore",
            help="Sets minimum score.  This is the matches minus the mismatches minus some sort of gap penalty.",
        ),
    ] = 30,
    minIdentity: Annotated[
        int | None,
        Option(
            "-minIdentity",
            help="Sets minimum sequence identity (in percent).",
            show_default="90 for nucleotide searches, 25 for protein or translated protein searches",
        ),
    ] = None,
    maxGap: Annotated[
        int,
        Option(
            "-maxGap",
            help="Sets the size of maximum gap between tiles in a clump.  Usually set from 0 to 3.  Only relevant for -minMatch > 1.",
        ),
    ] = 2,
    noHead: Annotated[
        bool,
        Option(
            "-noHead",
            help="Suppresses .psl header (so it's just a tab-separated file).",
        ),
    ] = False,
    makeOoc: Annotated[
        str | None,
        Option(
            "-makeOoc",
            help="Make overused tile file. Target needs to be complete genome.",
        ),
    ] = None,
    repMatch: Annotated[
        int | None,
        Option(
            "-repMatch",
            help="Sets the number of repetitions of a tile allowed before it is marked as overused. Typically this is 256 for -tileSize 12, 1024 for tile size 11, 4096 for tile size 10. Typically comes into play only with makeOoc. Also affected by stepSize: when stepSize is halved, repMatch is doubled to compensate.",
        ),
    ] = None,
    noSimpRepMask: Annotated[
        bool, Option("-noSimpRepMask", help="Suppresses simple repeat masking.")
    ] = False,
    mask: Annotated[
        str | None,
        Option(
            "-mask",
            help="Mask out repeats. Alignments won't be started in masked region but may extend through it in nucleotide searches. Masked areas are ignored entirely in protein or translated searches.",
            click_type=Choice(["lower", "upper", "out", "file.out"]),
        ),
    ] = None,
    qMask: Annotated[
        str | None,
        Option(
            "-qMask",
            help="Mask out repeats in query sequence. Similar to -mask above, but for query rather than target sequence.",
            click_type=Choice(["lower", "upper", "out", "file.out"]),
        ),
    ] = None,
    repeats: Annotated[
        str | None,
        Option(
            "-repeats",
            help="Type is same as mask types above.  Repeat bases will not be masked in any way, but matches in repeat areas will be reported separately from matches in other areas in the psl output.",
        ),
    ] = None,
    minRepDivergence: Annotated[
        int,
        Option(
            "-minRepDivergence",
            help="Minimum percent divergence of repeats to allow them to be unmasked. Only relevant for masking using RepeatMasker .out files.",
        ),
    ] = 15,
    dots: Annotated[
        int | None,
        Option(
            "-dots", help="Output dot every N sequences to show program's progress."
        ),
    ] = None,
    trimT: Annotated[bool, Option("-trimT", help="Trim leading poly-T")] = False,
    noTrimA: Annotated[
        bool, Option("-noTrimA", help="Don't trim trailing poly-A.")
    ] = False,
    trimHardA: Annotated[
        bool,
        Option(
            "-trimHardA",
            help="Remove poly-A tail from qSize as well as alignments in psl output.",
        ),
    ] = False,
    fastMap: Annotated[
        bool,
        Option(
            "-fastMap",
            help="Run for fast DNA/DNA remapping - not allowing introns, requiring high %ID, Query size must not exceed 5000.",
        ),
    ] = False,
    out: Annotated[
        str,
        Option(
            "-out",
            help="Controls output file format.",
            click_type=Choice(
                [
                    "psl",
                    "pslx",
                    "axt",
                    "maf",
                    "sim4",
                    "wublast",
                    "blast",
                    "blast8",
                    "blast9",
                ]
            ),
        ),
    ] = "psl",
    fine: Annotated[
        bool,
        Option(
            "-fine",
            help="For high-quality mRNAs, look harder for small initial and terminal exons. Not recommended for ESTs.",
        ),
    ] = False,
    maxIntron: Annotated[
        int, Option("-maxIntron", help="Sets maximum intron size.")
    ] = 750000,
    extendThroughN: Annotated[
        bool,
        Option(
            "-extendThroughN",
            help="Allows extension of alignment through large blocks of Ns.",
        ),
    ] = False,
):
    """Standalone BLAT v. 39x1 fast sequence search command line tool."""
    *query, output = query
    kwargs: dict[str, Any] = {}
    if prot:
        kwargs["prot"] = prot
    elif t != "dna":
        kwargs["t"] = t
    elif q != "dna":
        kwargs["q"] = q
    if ooc is not None:
        kwargs["ooc"] = ooc
    if tileSize is not None:
        kwargs["tileSize"] = tileSize
    if stepSize is not None:
        kwargs["stepSize"] = stepSize
    if oneOff != "0":
        kwargs["oneOff"] = oneOff
    if minMatch is not None:
        kwargs["minMatch"] = minMatch
    if minScore != 30:
        kwargs["minScore"] = minScore
    if minIdentity is not None:
        kwargs["minIdentity"] = minIdentity
    if maxGap != 2:
        kwargs["maxGap"] = maxGap
    if makeOoc is not None:
        kwargs["makeOoc"] = makeOoc
    if repMatch is not None:
        kwargs["repMatch"] = repMatch
    if noSimpRepMask:
        kwargs["noSimpRepMask"] = noSimpRepMask
    if mask is not None:
        kwargs["mask"] = mask
    if qMask is not None:
        kwargs["qMask"] = qMask
    if repeats is not None:
        kwargs["repeats"] = repeats
    if minRepDivergence != 15:
        kwargs["minRepDivergence"] = minRepDivergence
    if dots is not None:
        kwargs["dots"] = dots
    if trimT:
        kwargs["trimT"] = trimT
    if noTrimA:
        kwargs["noTrimA"] = noTrimA
    if trimHardA:
        kwargs["trimHardA"] = trimHardA
    if fastMap:
        kwargs["fastMap"] = fastMap
    if out != "psl":
        kwargs["out"] = out
    if fine:
        kwargs["fine"] = fine
    if maxIntron != 750000:
        kwargs["maxIntron"] = maxIntron
    if extendThroughN:
        kwargs["extendThroughN"] = extendThroughN
    o = blat(database, query, **kwargs)
    body = o.to_csv(sep="\t", header=False)
    content = body if noHead else f"{_PSL_HEAD}{body}"
    with fsspec.open(output, "wt") as f:
        f.write(content)  # pyright: ignore


@typechecked
def blat(
    database: str,
    query: str | list[str],
    **kwargs: Any,
) -> pd.DataFrame:
    """Standalone BLAT v. 39x1 fast sequence search command line tool.

    Args:
        database: a .fa, .nib or .2bit file, or a list of these files with one file name per line.
        query: Same as database but additionally support directly providing sequence.
        **kwargs: extra options will pass directly to original blat executable.

    Returns:
        a DataFrame object

    Raises:
        FileNotFoundError: if could not found blat in $BIOV_HOME/.bin:$PATH and could not rsync from ucsc.
        RuntimeError: if error occurred while rsync.
    """
    if isinstance(query, str):
        query = [query]
    bin_dir = settings.home / ".bin"
    bin_dir.mkdir(exist_ok=True, parents=True)
    blat_exec = shutil.which("blat", path=bin_dir)
    if blat_exec is None:
        blat_exec = shutil.which("blat")
    if blat_exec is None and _RSYNC_SUPPORTED:
        rsync = shutil.which("rsync")
        if rsync is None:
            raise FileNotFoundError("Sorry, biov could not find blat.")
        bin_dir = settings.home / ".bin"
        bin_dir.mkdir(exist_ok=True, parents=True)

        try:
            # ruff: noqa: S603
            subprocess.run(
                [
                    rsync,
                    "-aP",
                    f"rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/{_SUBDIR}/blat/blat",
                    bin_dir,
                ],
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError("Sorry, biov could not rsync blat from ucsc.") from e
        return blat(database, query, **kwargs)

    if blat_exec is None:
        raise FileNotFoundError("Sorry, biov could not find blat.")

    with ExitStack() as stack:
        database = _handle_fsspec_path(database, stack)
        queries = []
        for _q in query:
            _q = _handle_fsspec_path(_q, stack)
            _q = _handle_seq(_q, stack)
            queries.append(_q)
        if len(queries) > 1:
            names_file = stack.enter_context(NamedTemporaryFile("w+t"))
            names_file.writelines([f"{_q}\n" for _q in queries])
            query = [names_file.name]
        else:
            query = queries
        f = stack.enter_context(NamedTemporaryFile(mode="w+t"))
        subprocess.run(
            [
                blat_exec,
                "-noHead",
                *[
                    f"-{k}" if isinstance(v, bool) and v else f"-{k}={v}"
                    for k, v in kwargs.items()
                    if v is not None
                    if k != "noHead"
                ],
                database,
                *query,
                f.name,
            ]
        )
        return pd.read_table(
            f,
            names=_BLAT_COLUMNS,
        )
