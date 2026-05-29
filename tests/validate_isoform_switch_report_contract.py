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
        ncrna_annotations = tmp / "ncrna_annotations.tsv"
        interproscan = tmp / "interproscan.tsv"
        pfam_domtblout = tmp / "pfam.domtblout"
        coding_potential = tmp / "coding_potential.tsv"
        signalp = tmp / "signalp.tsv"
        deeptmhmm = tmp / "deeptmhmm.gff3"
        deeploc = tmp / "deeploc2.tsv"
        iupred = tmp / "iupred2a.tsv"
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
            "\t".join(
                [
                    "isoform_id",
                    "gene_id",
                    "gene_name",
                    "dIF",
                    "isoform_switch_q_value",
                    "isoform_fraction_control",
                    "isoform_fraction_test",
                    "ORF_length",
                    "NMD_status",
                    "coding_potential",
                ]
            )
            + "\n"
            + "\t".join(["tx_known", "geneA", "GeneA", "-0.42", "0.01", "0.72", "0.30", "4", "not_NMD", "coding"])
            + "\n"
            + "\t".join(["tx_novel", "geneA", "GeneA", "0.42", "0.01", "0.28", "0.70", "6", "NMD_sensitive", "coding"])
            + "\n"
            + "\t".join(["tx_context", "geneA", "GeneA", "0.03", "0.8", "0.05", "0.06", "3", "not_NMD", "coding"])
            + "\n"
            + "\t".join(["lnc_short", "geneB", "LncB", "-0.31", "0.02", "0.65", "0.34", "", "", "noncoding"])
            + "\n"
            + "\t".join(["lnc_long", "geneB", "LncB", "0.31", "0.02", "0.35", "0.66", "", "", "noncoding"])
            + "\n"
            + "\t".join(["pg_short", "geneC", "PseudoC", "-0.24", "0.03", "0.62", "0.38", "2", "", "low_potential"])
            + "\n"
            + "\t".join(["pg_long", "geneC", "PseudoC", "0.24", "0.03", "0.38", "0.62", "2", "", "low_potential"])
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
            ">tx_known\nATGGCCATGGCC\n>tx_novel\nATGGCCATGGCCCCCGGG\n>lnc_short\nATGATGATG\n>lnc_long\nATGATGATGCCCCCCCC\n>pg_short\nATGCGT\n>pg_long\nATGCGTAAAAAA\n",
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
                    "gene_biotype",
                    "transcript_biotype",
                ]
            )
            + "\n"
            + "\t".join(["tx_known", "geneA", "reference_compatible", "known", "known_compatible", "=", "protein_coding", "protein_coding"])
            + "\n"
            + "\t".join(["tx_novel", "geneA", "novel_isoform", "novel", "novel_isoform", "j", "protein_coding", "protein_coding"])
            + "\n"
            + "\t".join(["tx_context", "geneA", "reference_compatible", "known", "known_compatible", "c", "protein_coding", "protein_coding"])
            + "\n"
            + "\t".join(["lnc_short", "geneB", "reference_compatible", "known", "known_compatible", "=", "lncRNA", "lncRNA"])
            + "\n"
            + "\t".join(["lnc_long", "geneB", "novel_isoform", "novel", "novel_isoform", "j", "lncRNA", "lncRNA"])
            + "\n"
            + "\t".join(["pg_short", "geneC", "reference_compatible", "known", "known_compatible", "=", "processed_pseudogene", "processed_pseudogene"])
            + "\n"
            + "\t".join(["pg_long", "geneC", "novel_isoform", "novel", "novel_isoform", "j", "processed_pseudogene", "processed_pseudogene"])
            + "\n",
        )
        write(
            gtf,
            'chr1\tASPIS\texon\t1\t6\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_known"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\texon\t80\t85\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_known"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\tCDS\t1\t6\t.\t+\t0\tgene_id "geneA"; transcript_id "tx_known"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\tCDS\t80\t85\t.\t+\t0\tgene_id "geneA"; transcript_id "tx_known"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\texon\t1\t6\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_novel"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\texon\t150\t161\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_novel"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\tCDS\t1\t6\t.\t+\t0\tgene_id "geneA"; transcript_id "tx_novel"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\tCDS\t150\t161\t.\t+\t0\tgene_id "geneA"; transcript_id "tx_novel"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr1\tASPIS\texon\t1\t6\t.\t+\t.\tgene_id "geneA"; transcript_id "tx_context"; gene_name "GeneA"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";\n'
            'chr2\tASPIS\texon\t10\t18\t.\t+\t.\tgene_id "geneB"; transcript_id "lnc_short"; gene_name "LncB"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";\n'
            'chr2\tASPIS\texon\t10\t18\t.\t+\t.\tgene_id "geneB"; transcript_id "lnc_long"; gene_name "LncB"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";\n'
            'chr2\tASPIS\texon\t80\t86\t.\t+\t.\tgene_id "geneB"; transcript_id "lnc_long"; gene_name "LncB"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";\n'
            'chr3\tASPIS\texon\t20\t25\t.\t+\t.\tgene_id "geneC"; transcript_id "pg_short"; gene_name "PseudoC"; gene_biotype "processed_pseudogene"; transcript_biotype "processed_pseudogene";\n'
            'chr3\tASPIS\texon\t20\t25\t.\t+\t.\tgene_id "geneC"; transcript_id "pg_long"; gene_name "PseudoC"; gene_biotype "processed_pseudogene"; transcript_biotype "processed_pseudogene";\n'
            'chr3\tASPIS\texon\t110\t115\t.\t+\t.\tgene_id "geneC"; transcript_id "pg_long"; gene_name "PseudoC"; gene_biotype "processed_pseudogene"; transcript_biotype "processed_pseudogene";\n',
        )
        write(
            annotations,
            "\t".join(["isoform_id", "source", "feature_type", "feature_id", "feature_name", "start_aa", "end_aa", "description"])
            + "\n"
            + "\t".join(["tx_novel", "toy", "domain", "PF00001", "Catalytic fold", "2", "5", "toy switched domain"])
            + "\n",
        )
        write(
            interproscan,
            "\t".join(
                [
                    "tx_novel",
                    "md5",
                    "6",
                    "Pfam",
                    "PF00002",
                    "Transferase signature",
                    "3",
                    "6",
                    "1.0E-8",
                    "T",
                    "2026-01-01",
                    "IPR000002",
                    "Transferase domain",
                ]
            )
            + "\n",
        )
        write(
            pfam_domtblout,
            "# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias hmm from hmm to ali from ali to env from env to acc description\n"
            "PF00003 PF00003.1 120 tx_novel - 6 1e-20 80.0 0.0 1 1 1e-20 1e-20 80.0 0.0 1 8 2 5 2 5 0.99 Zinc finger domain\n",
        )
        write(
            coding_potential,
            "\t".join(["ID", "mRNA_size", "ORF_size", "coding_prob", "coding_label"])
            + "\n"
            + "\t".join(["tx_novel", "18", "18", "0.93", "coding"])
            + "\n",
        )
        write(
            signalp,
            "\t".join(["# ID", "Prediction", "OTHER", "SP(Sec/SPI)", "CS Position"])
            + "\n"
            + "\t".join(["tx_novel", "SP", "0.01", "0.99", "CS pos: 4-5. Pr: 0.90"])
            + "\n",
        )
        write(
            deeptmhmm,
            "##gff-version 3\n"
            "tx_novel\tDeepTMHMM\tTMhelix\t2\t5\t0.95\t.\t.\tID=tx_novel.tm1\n",
        )
        write(
            deeploc,
            "\t".join(["Protein_ID", "Localizations", "Signals", "Score"])
            + "\n"
            + "\t".join(["tx_novel", "Nucleus", "No signal", "0.77"])
            + "\n",
        )
        write(
            iupred,
            "\t".join(["protein_id", "position", "aa", "iupred2"])
            + "\n"
            + "\t".join(["tx_novel", "1", "M", "0.20"])
            + "\n"
            + "\t".join(["tx_novel", "2", "A", "0.82"])
            + "\n"
            + "\t".join(["tx_novel", "3", "M", "0.73"])
            + "\n"
            + "\t".join(["tx_novel", "4", "A", "0.30"])
            + "\n",
        )
        write(
            ncrna_annotations,
            "\t".join(
                [
                    "transcript_id",
                    "gene_id",
                    "chrom",
                    "start",
                    "end",
                    "source",
                    "feature_type",
                    "feature_id",
                    "feature_name",
                    "description",
                ]
            )
            + "\n"
            + "\t".join(
                [
                    "lnc_long",
                    "geneB",
                    "chr2",
                    "80",
                    "86",
                    "phyloP",
                    "conserved_exon",
                    "cons_exon_1",
                    "Conserved exon",
                    "overlaps gained lncRNA exon",
                ]
            )
            + "\n"
            + "\t".join(
                [
                    "lnc_long",
                    "geneB",
                    "chr2",
                    "82",
                    "84",
                    "RBPDB",
                    "rbp_motif",
                    "RBP1",
                    "RBP motif",
                    "motif in gained lncRNA exon",
                ]
            )
            + "\n"
            + "\t".join(
                [
                    "lnc_long",
                    "geneB",
                    "chr2",
                    "83",
                    "85",
                    "GENCODE",
                    "host_small_rna",
                    "sno1",
                    "embedded snoRNA",
                    "embedded small-RNA locus in gained lncRNA exon",
                ]
            )
            + "\n"
            + "\t".join(
                [
                    "lnc_long",
                    "geneB",
                    "chr2",
                    "80",
                    "86",
                    "GENCODE",
                    "antisense_overlap",
                    "GENEA-AS",
                    "antisense overlap",
                    "resource-backed antisense overlap",
                ]
            )
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
            "--ncrna-switch-table",
            str(outdir / "ncrna_switch_interpretation.tsv"),
            "--coding-switch-summary",
            str(outdir / "coding_switch_summary.tsv"),
            "--sequence-table",
            str(outdir / "switch_sequence_summary.tsv"),
            "--functional-annotation-table",
            str(outdir / "functional_annotation_summary.tsv"),
            "--plot-manifest",
            str(outdir / "switch_plot_manifest.tsv"),
            "--external-tool-manifest",
            str(outdir / "external_tool_manifest.tsv"),
            "--plots-pdf",
            str(outdir / "switch_plots.pdf"),
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
            ",".join(
                [
                    str(annotations),
                    str(interproscan),
                    str(pfam_domtblout),
                    str(coding_potential),
                    str(signalp),
                    str(deeptmhmm),
                    str(deeploc),
                    str(iupred),
                ]
            ),
            "--ncrna-annotation-tables",
            str(ncrna_annotations),
        ]
        completed = subprocess.run(command, check=False, capture_output=True, text=True)
        if completed.returncode:
            sys.stderr.write(completed.stdout)
            sys.stderr.write(completed.stderr)
            return completed.returncode

        candidates = read_rows(outdir / "switch_candidates.tsv")
        events = read_rows(outdir / "switch_event_summary.tsv")
        ncrna_rows = read_rows(outdir / "ncrna_switch_interpretation.tsv")
        coding_rows = read_rows(outdir / "coding_switch_summary.tsv")
        sequences = read_rows(outdir / "switch_sequence_summary.tsv")
        annotation_rows = read_rows(outdir / "functional_annotation_summary.tsv")
        plots = read_rows(outdir / "switch_plot_manifest.tsv")
        assert len(events) == 3, events
        assert {row["switch_biotype_class"] for row in events} >= {"coding", "noncoding"}, events
        assert coding_rows, "coding_switch_summary.tsv is empty"
        assert coding_rows[0]["gene_id"] == "geneA", coding_rows
        assert coding_rows[0]["coding_priority_tier"] == "high", coding_rows
        assert "protein_domain_gain_loss" in coding_rows[0]["coding_priority_reasons"], coding_rows
        assert coding_rows[0]["gained_domain"], coding_rows
        assert coding_rows[0]["gained_signal_peptide"], coding_rows
        assert coding_rows[0]["gained_transmembrane_region"], coding_rows
        assert "NMD_status_change" in coding_rows[0]["coding_priority_reasons"], coding_rows
        assert any(row["gene_id"] == "geneA" and row["coding_priority_score"] == coding_rows[0]["coding_priority_score"] for row in events), events
        assert {row["switch_role"] for row in candidates} >= {"switch_in", "switch_out"}, candidates
        assert any(row["gene_biotype"] == "lncRNA" and row["switch_biotype_class"] == "noncoding" for row in candidates)
        assert candidates[0]["switch_rank"] == "1"
        assert candidates[0]["reason_selected"]
        assert any(row["gene_id"] == "geneB" and row["interpretation_label"] == "noncoding_structure_change" for row in ncrna_rows), ncrna_rows
        assert any(row["gene_id"] == "geneB" and row["transcript_length_change"] for row in ncrna_rows), ncrna_rows
        pseudogene_events = [row for row in events if row["gene_id"] == "geneC"]
        assert pseudogene_events, events
        assert pseudogene_events[0]["switch_interpretation_label"] == "pseudogene_transcript_architecture_change", pseudogene_events
        pseudogene_rows = [row for row in ncrna_rows if row["gene_id"] == "geneC"]
        assert pseudogene_rows, ncrna_rows
        assert all(row["pseudogene_caution"] == "interpret_as_transcript_architecture_not_protein_consequence" for row in pseudogene_rows), pseudogene_rows
        assert all(row["interpretation_label"] == "pseudogene_transcript_architecture_change" for row in pseudogene_rows), pseudogene_rows
        lnc_long_rows = [row for row in ncrna_rows if row["gene_id"] == "geneB" and row["isoform_id"] == "lnc_long"]
        assert lnc_long_rows, ncrna_rows
        assert "cons_exon_1" in lnc_long_rows[0]["conserved_exon_change"], lnc_long_rows
        assert "RBP1" in lnc_long_rows[0]["motif_change"], lnc_long_rows
        assert "sno1" in lnc_long_rows[0]["host_smallrna_change"], lnc_long_rows
        assert "GENEA-AS" in lnc_long_rows[0]["resource_antisense_overlap"], lnc_long_rows
        assert "conserved_exon" in lnc_long_rows[0]["ncrna_resource_annotations"], lnc_long_rows
        assert any(row["isoform_id"] == "tx_novel" and row["aa_sequence"] == "MAMAPG" for row in sequences)
        assert any(row["isoform_id"] == "tx_novel" and row["gained_exon_coordinates"] for row in sequences)
        assert any(row["isoform_id"] == "tx_novel" and row["affected_aa_sequence"] for row in sequences)
        assert any(row["feature_name"] == "Catalytic fold" for row in annotation_rows)
        assert any(row["source"] == "interproscan:Pfam" and row["feature_id"] == "IPR000002" for row in annotation_rows)
        assert any(row["source"] == "hmmer_domtblout" and row["feature_id"] == "PF00003.1" for row in annotation_rows)
        assert any(row["feature_type"] == "coding_potential" and row["feature_name"] == "coding" for row in annotation_rows)
        assert any(row["feature_type"] == "signal_peptide" and row["end_aa"] == "4" for row in annotation_rows)
        assert any(row["feature_type"] == "transmembrane" and row["start_aa"] == "2" for row in annotation_rows)
        assert any(row["feature_type"] == "localization" and row["feature_name"] == "Nucleus" for row in annotation_rows)
        assert any(row["feature_type"] == "disorder" and row["start_aa"] == "2" and row["end_aa"] == "3" for row in annotation_rows)
        assert Path(plots[0]["plot_svg"]).exists(), plots
        assert Path(plots[0]["event_html"]).exists(), plots
        assert Path(plots[0]["nt_fasta"]).exists(), plots
        assert Path(plots[0]["aa_fasta"]).exists(), plots
        assert (outdir / "switch_plots.pdf").exists()
        assert (outdir / "external_tool_manifest.tsv").exists()
        assert (outdir / "index.html").exists()
        assert (outdir / "report.done").exists()
    print("isoform_switch_report\tok\tcandidate, sequence, annotation, SVG, and HTML outputs present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
