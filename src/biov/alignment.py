"""Align paired FASTQ files with BioV-managed BWA and sort the resulting BAM."""

import tempfile
from pathlib import Path

import pysam

from .config import settings
from .software import run_software


def align_paired_reads(
    reference_fasta: Path,
    fastq_r1: Path,
    fastq_r2: Path,
    bam_path: Path,
    *,
    mem_args: tuple[str, ...] = (),
    sort_args: tuple[str, ...] = (),
    build_index: bool = True,
) -> None:
    """Index a reference, align paired reads, and write a sorted local BAM.

    Args:
        reference_fasta: Local reference FASTA; BWA index files are written beside it.
        fastq_r1: Local first-mate FASTQ, optionally compressed.
        fastq_r2: Local second-mate FASTQ, optionally compressed.
        bam_path: Destination BAM; its parent directory must exist.
        mem_args: Additional native BWA-MEM options.
        sort_args: Native samtools sort options, passed through pysam.
        build_index: Write a BAM index; disable for query-name sorting.

    Raises:
        ValueError: If configured SSH dispatch would separate local temporary files.

    BWA uses the configured locked Pixi environment. The bundled environment
    currently supports Linux. This API owns local files: run the complete Python
    script on another host instead of dispatching its individual BWA commands.
    """
    if settings.execution_host is not None:
        raise ValueError(
            "Paired alignment uses local files; run the complete Python script "
            "on the execution host instead of configuring SSH command dispatch"
        )
    reference_fasta = reference_fasta.resolve()
    fastq_r1 = fastq_r1.resolve()
    fastq_r2 = fastq_r2.resolve()
    bam_path = bam_path.resolve()
    run_software(
        "bwa", ("index", str(reference_fasta)), cwd=Path.cwd()
    ).check_returncode()
    with tempfile.TemporaryDirectory() as temp_dir:
        sam_path = Path(temp_dir) / f"{bam_path.stem}.sam"
        run_software(
            "bwa",
            (
                "mem",
                *mem_args,
                "-o",
                str(sam_path),
                str(reference_fasta),
                str(fastq_r1),
                str(fastq_r2),
            ),
            cwd=Path.cwd(),
        ).check_returncode()
        pysam.sort(*sort_args, "-o", str(bam_path), str(sam_path), catch_stdout=False)
    if build_index:
        pysam.index(str(bam_path), catch_stdout=False)
