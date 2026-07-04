"""Validate the native rMATS DTU contract with a fake rmats.py executable."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("tmp_validation/rmats_contract")
INPUT = BASE / "input"


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def read_tsv(path: Path, required: set[str] | None = None) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is empty")
        if required:
            missing = required - set(reader.fieldnames)
            if missing:
                raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(completed.stdout + completed.stderr)


def write_common_inputs() -> dict[str, Path]:
    if BASE.exists():
        shutil.rmtree(BASE)
    INPUT.mkdir(parents=True)
    gtf = INPUT / "annotation.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\ttoy\tgene\t1\t200\t.\t+\t.\tgene_id "GENE1"; gene_name "Gene One";',
                'chr1\ttoy\ttranscript\t1\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1";',
                'chr1\ttoy\texon\t1\t80\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "1";',
                'chr1\ttoy\texon\t121\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX1"; exon_number "2";',
                'chr1\ttoy\ttranscript\t1\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2";',
                'chr1\ttoy\texon\t1\t50\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; exon_number "1";',
                'chr1\ttoy\texon\t151\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "TX2"; exon_number "2";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    samples = INPUT / "samples.tsv"
    write_tsv(
        samples,
        ["library_id", "assay", "project", "layout", "condition", "time_h"],
        [
            {"library_id": "ctrl_1", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "paired", "condition": "control", "time_h": "24"},
            {"library_id": "ctrl_2", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "paired", "condition": "control", "time_h": "24"},
            {"library_id": "treat_1", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "paired", "condition": "treated", "time_h": "24"},
            {"library_id": "treat_2", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "paired", "condition": "treated", "time_h": "24"},
        ],
    )
    transcript_counts = INPUT / "transcript_counts.tsv"
    write_tsv(
        transcript_counts,
        ["transcript_id", "ctrl_1", "ctrl_2", "treat_1", "treat_2"],
        [
            {"transcript_id": "TX1", "ctrl_1": "80", "ctrl_2": "78", "treat_1": "25", "treat_2": "28"},
            {"transcript_id": "TX2", "ctrl_1": "20", "ctrl_2": "22", "treat_1": "75", "treat_2": "72"},
        ],
    )
    transcript_metadata = INPUT / "transcript_metadata.tsv"
    write_tsv(
        transcript_metadata,
        ["transcript_id", "gene_id"],
        [
            {"transcript_id": "TX1", "gene_id": "GENE1"},
            {"transcript_id": "TX2", "gene_id": "GENE1"},
        ],
    )
    aligned_samples = INPUT / "aligned_samples.tsv"
    aligned_rows = []
    for sample_id in ["ctrl_1", "ctrl_2", "treat_1", "treat_2"]:
        bam = INPUT / f"{sample_id}.bam"
        bam.write_text("", encoding="utf-8")
        aligned_rows.append({"library_id": sample_id, "bam": str(bam), "layout": "paired"})
    write_tsv(aligned_samples, ["library_id", "bam", "layout"], aligned_rows)
    return {
        "gtf": gtf,
        "samples": samples,
        "transcript_counts": transcript_counts,
        "transcript_metadata": transcript_metadata,
        "aligned_samples": aligned_samples,
    }


def write_fake_rmats() -> Path:
    helper = INPUT / "fake_rmats.py"
    helper.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "def value(flag):",
                "    return Path(args[args.index(flag) + 1])",
                "outdir = value('--od')",
                "tmpdir = value('--tmp')",
                "if tmpdir.exists() and any(tmpdir.iterdir()):",
                "    sys.stderr.write('stale rmats temp directory was not cleaned\\n')",
                "    sys.exit(8)",
                "outdir.mkdir(parents=True, exist_ok=True)",
                "(outdir / 'SE.MATS.JC.txt').write_text(",
                "    'ID\\tGeneID\\tgeneSymbol\\tPValue\\tFDR\\tIncLevelDifference\\n'",
                "    'EVENT1\\tGENE1\\tGene One\\t0.001\\t0.01\\t0.35\\n'",
                "    'EVENT2\\tGENE1\\tGene One\\t0.5\\t0.8\\t-0.10\\n',",
                "    encoding='utf-8',",
                ")",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return helper


def write_failing_rmats() -> Path:
    helper = INPUT / "fake_rmats_fail.py"
    helper.write_text(
        "\n".join(
            [
                "import sys",
                "sys.stderr.write('rmats simulated alignment error\\n')",
                "sys.exit(9)",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return helper


def main() -> int:
    snakefile = Path("Snakefile").read_text(encoding="utf-8")
    if "--rmats-extra-args {params.rmats_extra_args:q}" in snakefile:
        raise ValueError("Snakefile emits --rmats-extra-args even when the configured value is empty")
    if 'rmats_extra_args=optional_shell_arg("--rmats-extra-args"' not in snakefile:
        raise ValueError("Snakefile does not guard empty rMATS extra args with optional_shell_arg")

    paths = write_common_inputs()
    fake_rmats = write_fake_rmats()
    failing_rmats = write_failing_rmats()
    dtu_dir = BASE / "dtu"
    run_command(
        [
            sys.executable,
            "workflow/scripts/plan_rnaseq_dtu.py",
            "--samples",
            str(paths["samples"]),
            "--transcript-counts",
            str(paths["transcript_counts"]),
            "--transcript-metadata",
            str(paths["transcript_metadata"]),
            "--annotation-gtf",
            str(paths["gtf"]),
            "--output",
            str(dtu_dir / "dtu_plan.tsv"),
            "--done",
            str(dtu_dir / "dtu.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "rMATS",
            "--candidate-methods",
            "rMATS",
            "--contrast-by",
            "time_h",
            "--min-replicates",
            "1",
        ]
    )
    plan_rows = read_tsv(dtu_dir / "dtu_plan.tsv", {"method", "status", "contrast_id"})
    contrast_id = "treated_vs_control__time_h_24"
    if plan_rows[0]["method"] != "rMATS" or plan_rows[0]["status"] != "ready" or plan_rows[0]["contrast_id"] != contrast_id:
        raise ValueError(f"unexpected rMATS plan row: {plan_rows}")

    blocked_dir = BASE / "dtu_blocked"
    run_command(
        [
            sys.executable,
            "workflow/scripts/run_rnaseq_dtu_methods.py",
            "--plan",
            str(dtu_dir / "dtu_plan.tsv"),
            "--samples",
            str(paths["samples"]),
            "--aligned-samples",
            str(paths["aligned_samples"]),
            "--transcript-counts",
            str(paths["transcript_counts"]),
            "--transcript-metadata",
            str(paths["transcript_metadata"]),
            "--annotation-gtf",
            str(paths["gtf"]),
            "--outdir",
            str(blocked_dir / "methods"),
            "--manifest",
            str(blocked_dir / "dtu_method_manifest.tsv"),
            "--done",
            str(blocked_dir / "dtu_methods.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "rMATS",
            "--methods",
            "rMATS",
            "--contrast-id",
            contrast_id,
            "--rmats-executable",
            f'"{sys.executable}" "{fake_rmats}"',
        ]
    )
    blocked_rows = read_tsv(blocked_dir / "dtu_method_manifest.tsv", {"method", "status", "reason"})
    if blocked_rows[0]["status"] != "blocked" or "read_length" not in blocked_rows[0]["reason"]:
        raise ValueError(f"rMATS missing read-length was not blocked clearly: {blocked_rows}")

    failed_dir = BASE / "dtu_failed"
    run_command(
        [
            sys.executable,
            "workflow/scripts/run_rnaseq_dtu_methods.py",
            "--plan",
            str(dtu_dir / "dtu_plan.tsv"),
            "--samples",
            str(paths["samples"]),
            "--aligned-samples",
            str(paths["aligned_samples"]),
            "--transcript-counts",
            str(paths["transcript_counts"]),
            "--transcript-metadata",
            str(paths["transcript_metadata"]),
            "--annotation-gtf",
            str(paths["gtf"]),
            "--outdir",
            str(failed_dir / "methods"),
            "--manifest",
            str(failed_dir / "dtu_method_manifest.tsv"),
            "--done",
            str(failed_dir / "dtu_methods.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "rMATS",
            "--methods",
            "rMATS",
            "--contrast-id",
            contrast_id,
            "--rmats-executable",
            f'"{sys.executable}" "{failing_rmats}"',
            "--rmats-read-length",
            "100",
            "--rmats-lib-type",
            "fr-unstranded",
        ]
    )
    failed_rows = read_tsv(failed_dir / "dtu_method_manifest.tsv", {"method", "status", "reason", "stderr"})
    if failed_rows[0]["status"] != "failed" or "rmats simulated alignment error" not in failed_rows[0]["reason"]:
        raise ValueError(f"rMATS failure was not preserved with stderr context: {failed_rows}")
    failed_done = read_tsv(failed_dir / "dtu_methods.done", {"status", "failed", "reason"})
    if failed_done[0]["status"] != "failed" or failed_done[0]["failed"] != "1":
        raise ValueError(f"rMATS failed done status was not written: {failed_done}")
    run_command(
        [
            sys.executable,
            "workflow/scripts/merge_status_manifests.py",
            "--kind",
            "dtu",
            "--manifest",
            str(failed_dir / "merged_dtu_method_manifest.tsv"),
            "--done",
            str(failed_dir / "merged_dtu_methods.done"),
            str(failed_dir / "dtu_method_manifest.tsv"),
        ]
    )
    merged_done = read_tsv(failed_dir / "merged_dtu_methods.done", {"status", "failed", "reason"})
    if merged_done[0]["status"] != "failed" or merged_done[0]["failed"] != "1":
        raise ValueError(f"failed rMATS status was not preserved by DTU merge: {merged_done}")

    method_dir = dtu_dir / "methods" / "rmats" / contrast_id
    stale_tmp = method_dir / "tmp"
    stale_tmp.mkdir(parents=True, exist_ok=True)
    (stale_tmp / "old.rmats").write_text("stale\n", encoding="utf-8")
    stale_output = method_dir / "rmats_output"
    stale_output.mkdir(parents=True, exist_ok=True)
    (stale_output / "old.txt").write_text("stale\n", encoding="utf-8")
    run_command(
        [
            sys.executable,
            "workflow/scripts/run_rnaseq_dtu_methods.py",
            "--plan",
            str(dtu_dir / "dtu_plan.tsv"),
            "--samples",
            str(paths["samples"]),
            "--aligned-samples",
            str(paths["aligned_samples"]),
            "--transcript-counts",
            str(paths["transcript_counts"]),
            "--transcript-metadata",
            str(paths["transcript_metadata"]),
            "--annotation-gtf",
            str(paths["gtf"]),
            "--outdir",
            str(dtu_dir / "methods"),
            "--manifest",
            str(dtu_dir / "dtu_method_manifest.tsv"),
            "--done",
            str(dtu_dir / "dtu_methods.done"),
            "--project",
            "ASPIS_CONTRACT",
            "--method",
            "rMATS",
            "--methods",
            "rMATS",
            "--contrast-id",
            contrast_id,
            "--rmats-executable",
            f'"{sys.executable}" "{fake_rmats}"',
            "--rmats-read-length",
            "100",
            "--rmats-lib-type",
            "fr-unstranded",
        ]
    )
    method_rows = read_tsv(
        dtu_dir / "dtu_method_manifest.tsv",
        {"method", "status", "reason", "transcript_results", "summary", "standardized_results", "standardized_status", "standardized_result_count"},
    )
    row = method_rows[0]
    if row["method"] != "rMATS" or row["status"] != "completed" or row["standardized_status"] != "ok":
        raise ValueError(f"rMATS method row was not completed and standardized: {method_rows}")
    if "aligned BAMs" not in row["reason"]:
        raise ValueError(f"rMATS reason did not describe native mode: {row}")
    if (stale_output / "old.txt").exists():
        raise ValueError("rMATS stale output directory was not cleaned before rerun")
    standardized = read_tsv(Path(row["standardized_results"]), {"method", "feature_id", "gene_id", "event_type", "delta_psi", "pvalue", "padj", "direction"})
    if len(standardized) != 2 or standardized[0]["method"] != "rMATS":
        raise ValueError(f"standardized rMATS rows were not preserved: {standardized}")
    if standardized[0]["feature_id"] != "EVENT1" or standardized[0]["event_type"] != "SE" or standardized[0]["direction"] != "increased_usage":
        raise ValueError(f"standardized rMATS row lost event fields: {standardized}")
    summary = read_tsv(Path(row["summary"]), {"n_events", "n_significant", "event_mode"})
    if summary[0]["n_events"] != "2" or summary[0]["n_significant"] != "1" or summary[0]["event_mode"] != "junction":
        raise ValueError(f"rMATS summary was not populated: {summary}")

    run_command(
        [
            sys.executable,
            "workflow/scripts/render_rnaseq_dtu_plots.py",
            "--method-manifest",
            str(dtu_dir / "dtu_method_manifest.tsv"),
            "--outdir",
            str(dtu_dir / "plots"),
            "--manifest",
            str(dtu_dir / "plots" / "dtu_plot_manifest.tsv"),
            "--done",
            str(dtu_dir / "plots" / "dtu_plots.done"),
            "--top-n",
            "5",
            "--max-points",
            "100",
        ]
    )
    plot_rows = read_tsv(dtu_dir / "plots" / "dtu_plot_manifest.tsv", {"method", "status", "usage_plot", "feature_plot"})
    if plot_rows[0]["method"] != "rMATS" or plot_rows[0]["status"] != "ok":
        raise ValueError(f"rMATS plots were not ok: {plot_rows}")
    usage_svg = Path(plot_rows[0]["usage_plot"]).read_text(encoding="utf-8")
    if "Top rMATS genes: event detail" not in usage_svg or "Gene One" not in usage_svg or "SE #EVENT1" not in usage_svg:
        raise ValueError(f"rMATS delta PSI plot was not rendered correctly: {plot_rows}")
    feature_svg = Path(plot_rows[0]["feature_plot"]).read_text(encoding="utf-8")
    if "Ranked rMATS event candidates" not in feature_svg or "SE skipped exon" not in feature_svg:
        raise ValueError(f"rMATS candidate plot did not explain event types: {plot_rows}")
    print("rMATS contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
