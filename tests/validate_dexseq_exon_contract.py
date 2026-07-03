#!/usr/bin/env python3
"""Validate the native exon-bin DEXSeq DTU contract with fake helper tools."""

from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path


BASE = Path("tmp_validation/dexseq_exon_contract")
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
                'chr1\ttoy\tgene\t1\t200\t.\t+\t.\tgene_id "GENE1"; gene_biotype "protein_coding";',
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
            {"library_id": "ctrl_1", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "single", "condition": "control", "time_h": "24"},
            {"library_id": "ctrl_2", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "single", "condition": "control", "time_h": "24"},
            {"library_id": "treat_1", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "single", "condition": "treated", "time_h": "24"},
            {"library_id": "treat_2", "assay": "rnaseq", "project": "ASPIS_CONTRACT", "layout": "single", "condition": "treated", "time_h": "24"},
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
        aligned_rows.append({"library_id": sample_id, "bam": str(bam), "layout": "single"})
    write_tsv(aligned_samples, ["library_id", "bam", "layout"], aligned_rows)
    return {
        "gtf": gtf,
        "samples": samples,
        "transcript_counts": transcript_counts,
        "transcript_metadata": transcript_metadata,
        "aligned_samples": aligned_samples,
    }


def write_fake_helpers() -> dict[str, Path]:
    prepare = INPUT / "fake_dexseq_prepare_annotation.py"
    prepare.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "out = Path(sys.argv[-1])",
                "out.parent.mkdir(parents=True, exist_ok=True)",
                "out.write_text(",
                "    'chr1\\tfake\\texonic_part\\t1\\t80\\t.\\t+\\t.\\tgene_id \"GENE1\"; exonic_part_number \"001\";\\n'",
                "    'chr1\\tfake\\texonic_part\\t121\\t200\\t.\\t+\\t.\\tgene_id \"GENE1\"; exonic_part_number \"002\";\\n',",
                "    encoding='utf-8',",
                ")",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    count = INPUT / "fake_dexseq_count.py"
    count.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "bam = Path(sys.argv[-2]).name",
                "out = Path(sys.argv[-1])",
                "out.parent.mkdir(parents=True, exist_ok=True)",
                "if bam.startswith('ctrl'):",
                "    rows = [('GENE1:E001', 80), ('GENE1:E002', 20)]",
                "else:",
                "    rows = [('GENE1:E001', 25), ('GENE1:E002', 75)]",
                "out.write_text(''.join(f'{feature}\\t{count}\\n' for feature, count in rows) + '__no_feature\\t0\\n', encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    runner = INPUT / "fake_dexseq_exon_runner.py"
    runner.write_text(
        "\n".join(
            [
                "import sys",
                "from pathlib import Path",
                "args = sys.argv[1:]",
                "def value(flag):",
                "    return Path(args[args.index(flag) + 1])",
                "gene = value('--gene-results')",
                "feature = value('--feature-results')",
                "summary = value('--summary')",
                "gene.parent.mkdir(parents=True, exist_ok=True)",
                "gene.write_text('gene_id\\tn_features\\tmin_pvalue\\tmin_padj\\tstatus\\nGENE1\\t2\\t0.001\\t0.01\\tok\\n', encoding='utf-8')",
                "feature.write_text(",
                "    'gene_id\\tfeature_id\\tstatistic\\tlog2_fold_change\\tpvalue\\tpadj\\tevent_type\\tmean_usage_control\\tmean_usage_test\\tdelta_usage\\tstatus\\n'",
                "    'GENE1\\t\"GENE1\"\"E001\"\\t12.5\\t-1.2\\t0.001\\t0.01\\texon_bin_usage\\tNA\\tNA\\tNA\\tok\\n'",
                "    'GENE1\\t\"GENE1\"\"E002\"\\t11.7\\t1.4\\t0.002\\t0.02\\texon_bin_usage\\tNA\\tNA\\tNA\\tok\\n',",
                "    encoding='utf-8',",
                ")",
                "summary.write_text('status\\treason\\tn_input_exon_bins\\tn_tested_genes\\tn_usage_exon_bins\\tcontrol_label\\ttest_label\\nok\\ttrue exon-bin fake DEXSeq\\t2\\t1\\t2\\tcontrol\\ttreated\\n', encoding='utf-8')",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {"prepare": prepare, "count": count, "runner": runner}


def main() -> int:
    snakefile = Path("Snakefile").read_text(encoding="utf-8")
    if 'method="DRIMSeq|DEXSeq|DEXSeqExon|SUPPA2|rMATS"' not in snakefile:
        raise ValueError("Snakefile method wildcard constraint does not allow DEXSeqExon")

    paths = write_common_inputs()
    helpers = write_fake_helpers()
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
            "DEXSeqExon",
            "--candidate-methods",
            "DEXSeqExon",
            "--contrast-by",
            "time_h",
            "--min-replicates",
            "1",
        ]
    )
    plan_rows = read_tsv(dtu_dir / "dtu_plan.tsv", {"method", "status", "contrast_id"})
    contrast_id = "treated_vs_control__time_h_24"
    if plan_rows[0]["method"] != "DEXSeqExon" or plan_rows[0]["status"] != "ready" or plan_rows[0]["contrast_id"] != contrast_id:
        raise ValueError(f"unexpected DEXSeqExon plan row: {plan_rows}")

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
            "DEXSeqExon",
            "--methods",
            "DEXSeqExon",
            "--contrast-id",
            contrast_id,
            "--rscript",
            sys.executable,
            "--dexseq-exon-script",
            str(helpers["runner"]),
            "--dexseq-prepare-annotation-command",
            f'"{sys.executable}" "{helpers["prepare"]}"',
            "--dexseq-count-command",
            f'"{sys.executable}" "{helpers["count"]}"',
            "--dtu-min-count",
            "1",
            "--dtu-min-samples",
            "1",
            "--dtu-min-gene-count",
            "1",
            "--dtu-min-transcripts-per-gene",
            "2",
        ]
    )
    method_rows = read_tsv(
        dtu_dir / "dtu_method_manifest.tsv",
        {"method", "status", "reason", "transcript_results", "standardized_results", "standardized_status", "standardized_result_count"},
    )
    row = method_rows[0]
    if row["method"] != "DEXSeqExon" or row["status"] != "completed" or row["standardized_status"] != "ok":
        raise ValueError(f"DEXSeqExon method row was not completed and standardized: {method_rows}")
    if "exon-bin" not in row["reason"]:
        raise ValueError(f"DEXSeqExon reason did not describe exon-bin mode: {row}")

    method_dir = dtu_dir / "methods" / "dexseq_exon" / contrast_id
    counts = read_tsv(method_dir / "dexseq_exon_counts.tsv", {"exon_bin_id", "ctrl_1", "treat_1"})
    if {item["exon_bin_id"] for item in counts} != {"GENE1:E001", "GENE1:E002"}:
        raise ValueError(f"DEXSeqExon count matrix lost exon-bin identifiers: {counts}")
    metadata = read_tsv(method_dir / "dexseq_exon_metadata.tsv", {"exon_bin_id", "gene_id"})
    if {item["gene_id"] for item in metadata} != {"GENE1"}:
        raise ValueError(f"DEXSeqExon metadata did not map exon bins to genes: {metadata}")
    standardized = read_tsv(Path(row["standardized_results"]), {"method", "feature_id", "gene_id", "event_type", "delta_psi", "pvalue", "padj"})
    exon_rows = [item for item in standardized if item["event_type"] == "exon_bin_usage"]
    if len(exon_rows) != 2 or {item["feature_id"] for item in exon_rows} != {'GENE1"E001', 'GENE1"E002'}:
        raise ValueError(f"standardized DEXSeqExon exon-bin rows were not preserved: {standardized}")

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
    plot_rows = read_tsv(dtu_dir / "plots" / "dtu_plot_manifest.tsv", {"method", "status", "overview_plot", "usage_plot", "feature_plot"})
    if plot_rows[0]["method"] != "DEXSeqExon" or plot_rows[0]["status"] != "ok":
        raise ValueError(f"DEXSeqExon plots were not ok: {plot_rows}")
    usage_svg = Path(plot_rows[0]["usage_plot"]).read_text(encoding="utf-8")
    if "Top DEXSeqExon genes: exon-bin detail" not in usage_svg or "exon bin E002" not in usage_svg or "log2FC" not in usage_svg:
        raise ValueError(f"DEXSeqExon usage plot did not include exon-bin features: {plot_rows}")
    if '""' in usage_svg:
        raise ValueError(f"DEXSeqExon usage plot exposed raw quoted DEXSeq bin identifiers: {plot_rows}")
    feature_svg = Path(plot_rows[0]["feature_plot"]).read_text(encoding="utf-8")
    if "Ranked DEXSeqExon exon-bin candidates" not in feature_svg or "exon bin E002" not in feature_svg:
        raise ValueError(f"DEXSeqExon candidate plot did not include cross-gene exon-bin features: {plot_rows}")
    print("DEXSeqExon contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
