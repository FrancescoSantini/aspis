#!/usr/bin/env python3
"""Contract tests for real-project preflight validation."""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "workflow/scripts/validate_project_inputs.py"


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def touch(path: Path, text: str = "placeholder\n") -> Path:
    return write(path, text)


def make_fastq(path: Path) -> Path:
    return write(path, "@r1\nACGTACGT\n+\nFFFFFFFF\n")


def run_preflight(config: Path, assay: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), "--config", str(config), "--assay", assay],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def rnaseq_intake(tmp: Path, rows: list[tuple[str, str, str, str]]) -> Path:
    lines = ["library_id\tproject\tinput_1\tinput_2\tassay_hint\tcondition\ttime_h"]
    for library_id, condition, r1, r2 in rows:
        lines.append(f"{library_id}\tRNASEQ_TEST\t{r1}\t{r2}\trnaseq\t{condition}\t24")
    return write(tmp / "rnaseq_intake.tsv", "\n".join(lines) + "\n")


def smallrna_intake(tmp: Path, rows: list[tuple[str, str, str, str]]) -> Path:
    lines = ["library_id\tproject\tinput_1\tinput_2\tassay_hint\tcondition\ttime_h"]
    for library_id, condition, r1, r2 in rows:
        lines.append(f"{library_id}\tSMALLRNA_TEST\t{r1}\t{r2}\tsmallrna\t{condition}\t24")
    return write(tmp / "smallrna_intake.tsv", "\n".join(lines) + "\n")


def rnaseq_config(tmp: Path, intake: Path, genome: Path, annotation: Path) -> Path:
    return write(
        tmp / "rnaseq.yaml",
        f"""
intake: {intake}
design:
  condition_col: condition
  control_label: control
  min_condition_groups: 2
rnaseq_alignment:
  run: true
  aligner: star
  reference_fasta: {genome}
  star_genome_dir: {tmp / "star"}
  annotation_gtf: {annotation}
rnaseq_quantification:
  run: true
  reference_fasta: {genome}
  annotation_gtf: {annotation}
rnaseq_differential:
  run: true
  levels: [gene, transcript]
  condition_col: condition
  control_label: control
  contrast_by: [time_h]
  min_replicates_per_group: 2
""".lstrip(),
    )


def smallrna_config(tmp: Path, intake: Path, mirbase: Path, contaminants: Path, genome: Path, annotation: Path) -> Path:
    return write(
        tmp / "smallrna.yaml",
        f"""
intake: {intake}
design:
  condition_col: condition
  control_label: control
  min_condition_groups: 2
smallrna:
  run: true
  preprocess_run: true
  depletion_run: true
  alignment_run: true
  quantification_run: true
  differential_run: true
  reference_run: true
  mirbase_fasta: {mirbase}
  build_contaminant_index: true
  contaminant_fasta: {contaminants}
  build_bowtie_index: true
  residual_run: true
  build_residual_genome_index: true
  residual_genome_fasta: {genome}
  residual_annotation_gtf: {annotation}
  condition_col: condition
  control_label: control
  contrast_by: [time_h]
  min_replicates_per_group: 2
""".lstrip(),
    )


def assert_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        raise AssertionError(f"{label} should have passed\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")


def assert_failure_contains(result: subprocess.CompletedProcess[str], label: str, expected: str) -> None:
    if result.returncode == 0:
        raise AssertionError(f"{label} should have failed\nSTDOUT:\n{result.stdout}")
    output = result.stdout + result.stderr
    if expected not in output:
        raise AssertionError(f"{label} did not mention {expected!r}\n{output}")


def main() -> int:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        genome = touch(tmp / "genome.fa")
        annotation = touch(tmp / "annotation.gtf")
        mirbase = touch(tmp / "mature.fa")
        contaminants = touch(tmp / "contaminants.fa")
        r1 = [make_fastq(tmp / f"sample_{idx}_R1.fastq") for idx in range(4)]
        r2 = [make_fastq(tmp / f"sample_{idx}_R2.fastq") for idx in range(4)]

        valid_rnaseq = rnaseq_intake(
            tmp,
            [
                ("control_1", "control", str(r1[0]), str(r2[0])),
                ("control_2", "control", str(r1[1]), str(r2[1])),
                ("treated_1", "treated", str(r1[2]), str(r2[2])),
                ("treated_2", "treated", str(r1[3]), str(r2[3])),
            ],
        )
        assert_success(run_preflight(rnaseq_config(tmp, valid_rnaseq, genome, annotation), "rnaseq"), "valid RNA-seq")

        low_rep_rnaseq = rnaseq_intake(
            tmp,
            [
                ("control_only", "control", str(r1[0]), str(r2[0])),
                ("treated_only", "treated", str(r1[1]), str(r2[1])),
            ],
        )
        assert_failure_contains(
            run_preflight(rnaseq_config(tmp, low_rep_rnaseq, genome, annotation), "rnaseq"),
            "low-replicate RNA-seq",
            "has 1 sample(s)",
        )

        missing_fastq = rnaseq_intake(
            tmp,
            [
                ("control_1", "control", str(tmp / "missing_R1.fastq"), str(r2[0])),
                ("control_2", "control", str(r1[1]), str(r2[1])),
                ("treated_1", "treated", str(r1[2]), str(r2[2])),
                ("treated_2", "treated", str(r1[3]), str(r2[3])),
            ],
        )
        assert_failure_contains(
            run_preflight(rnaseq_config(tmp, missing_fastq, genome, annotation), "rnaseq"),
            "missing FASTQ",
            "input_1 does not exist",
        )

        valid_smallrna = smallrna_intake(
            tmp,
            [
                ("control_1", "control", str(r1[0]), ""),
                ("control_2", "control", str(r1[1]), ""),
                ("treated_1", "treated", str(r1[2]), ""),
                ("treated_2", "treated", str(r1[3]), ""),
            ],
        )
        assert_success(
            run_preflight(smallrna_config(tmp, valid_smallrna, mirbase, contaminants, genome, annotation), "smallrna"),
            "valid smallRNA",
        )

        paired_smallrna = smallrna_intake(
            tmp,
            [
                ("control_1", "control", str(r1[0]), str(r2[0])),
                ("control_2", "control", str(r1[1]), ""),
                ("treated_1", "treated", str(r1[2]), ""),
                ("treated_2", "treated", str(r1[3]), ""),
            ],
        )
        assert_failure_contains(
            run_preflight(smallrna_config(tmp, paired_smallrna, mirbase, contaminants, genome, annotation), "smallrna"),
            "paired smallRNA",
            "single-end",
        )

    print("project preflight contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
