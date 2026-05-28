#!/usr/bin/env python3
"""Contract test for the isoform-switch event report renderer."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    script = repo / "workflow" / "scripts" / "render_isoform_switch_report.py"
    with tempfile.TemporaryDirectory(prefix="aspis_isoform_switch_") as tmp_text:
        tmp = Path(tmp_text)
        manifest = tmp / "isoform_switch_manifest.tsv"
        detailed = tmp / "contrast" / "isoform_switch_detailed.tsv"
        consequences = tmp / "contrast" / "switch_consequences.tsv"
        nt_fasta = tmp / "contrast" / "switch_nt.fa"
        aa_fasta = tmp / "contrast" / "switch_aa.fa"
        metadata = tmp / "transcript_metadata.tsv"
        gtf = tmp / "annotated.gtf"
        annotations = tmp / "domains.tsv"
        outdir = tmp / "report"

        write(
            manifest,
            "\t".join(["contrast_id", "status", "reason", "detailed", "consequences", "nt_fasta", "aa_fasta"])
            + "\n"
            + "\t".join(["treated_vs_control", "ok", "", str(detailed), str(consequences), str(nt_fasta), str(aa_fasta)])
            + "\n",
        )
        write(
            detailed,
            "\t".join(["isoform_id", "gene_id", "gene_name", "dIF", "isoform_switch_q_value"])
            + "\n"
            + "\t".join(["tx_known", "geneA", "GeneA", "-0.42", "0.01"])
            + "\n"
            + "\t".join(["tx_novel", "geneA", "GeneA", "0.42", "0.01"])
            + "\n"
            + "\t".join(["tx_context", "geneA", "GeneA", "0.03", "0.8"])
            + "\n",
        )
        write(
            consequences,
            "\t".join(["gene_id", "isoform_id", "consequence", "category"])
            + "\n"
            + "\t".join(["geneA", "tx_novel", "ORF_seq_similarity", "coding"])
            + "\n",
        )
        write(
            nt_fasta,
            ">tx_known\nATGGCCATGGCC\n>tx_novel\nATGGCCATGGCCCCCGGG\n",
        )
        write(
            aa_fasta,
            ">tx_known\nMAMA\n>tx_novel\nMAMAPG\n",
        )
        write(
            metadata,
            "\t".join(
                [
                    "transcript_id",
                    "gene_id",
                    "transcript_discovery_class",
                    "transcript_novelty",
                    "transcript_plot_group",
                    "gffcompare_class_code",
                ]
            )
            + "\n"
            + "\t".join(["tx_known", "geneA", "reference_compatible", "known", "known_compatible", "="])
            + "\n"
            + "\t".join(["tx_novel", "geneA", "novel_isoform", "novel", "novel_isoform", "j"])
            + "\n"
            + "\t".join(["tx_context", "geneA", "reference_compatible", "known", "known_compatible", "c"])
            + "\n",
        )
        write(
            gtf,
            'chr1\tASPIS\texon\t1\t30\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_known"; gene_name "GeneA";\n'
            'chr1\tASPIS\texon\t80\t120\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_known"; gene_name "GeneA";\n'
            'chr1\tASPIS\texon\t1\t30\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_novel"; gene_name "GeneA";\n'
            'chr1\tASPIS\texon\t150\t210\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_novel"; gene_name "GeneA";\n'
            'chr1\tASPIS\texon\t1\t30\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_context"; gene_name "GeneA";\n',
        )
        write(
            annotations,
            "\t".join(["isoform_id", "source", "feature_type", "feature_id", "feature_name", "start_aa", "end_aa", "description"])
            + "\n"
            + "\t".join(["tx_novel", "toy", "domain", "PF00001", "Catalytic fold", "2", "5", "toy switched domain"])
            + "\n",
        )

        command = [
            sys.executable,
            str(script),
            "--manifest",
            str(manifest),
            "--transcript-metadata",
            str(metadata),
            "--annotated-gtf",
            str(gtf),
            "--outdir",
            str(outdir),
            "--candidate-table",
            str(outdir / "switch_candidates.tsv"),
            "--event-summary",
            str(outdir / "switch_event_summary.tsv"),
            "--sequence-table",
            str(outdir / "switch_sequence_summary.tsv"),
            "--functional-annotation-table",
            str(outdir / "functional_annotation_summary.tsv"),
            "--plot-manifest",
            str(outdir / "switch_plot_manifest.tsv"),
            "--html",
            str(outdir / "index.html"),
            "--done",
            str(outdir / "report.done"),
            "--padj",
            "0.05",
            "--dif",
            "0.1",
            "--top-n",
            "5",
            "--functional-annotation-tables",
            str(annotations),
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode:
            sys.stderr.write(completed.stdout)
            sys.stderr.write(completed.stderr)
            return completed.returncode

        candidates = read_rows(outdir / "switch_candidates.tsv")
        events = read_rows(outdir / "switch_event_summary.tsv")
        sequences = read_rows(outdir / "switch_sequence_summary.tsv")
        annotation_rows = read_rows(outdir / "functional_annotation_summary.tsv")
        plots = read_rows(outdir / "switch_plot_manifest.tsv")
        assert len(events) == 1, events
        assert {row["switch_role"] for row in candidates} >= {"switch_in", "switch_out"}, candidates
        assert any(row["isoform_id"] == "tx_novel" and row["aa_sequence"] == "MAMAPG" for row in sequences)
        assert any(row["feature_name"] == "Catalytic fold" for row in annotation_rows)
        assert Path(plots[0]["plot_svg"]).exists(), plots
        assert Path(plots[0]["event_html"]).exists(), plots
        assert (outdir / "index.html").exists()
        assert (outdir / "report.done").exists()
    print("isoform_switch_report\tok\tcandidate, sequence, annotation, SVG, and HTML outputs present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
