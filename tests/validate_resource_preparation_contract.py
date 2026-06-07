#!/usr/bin/env python3
"""Contract test for offline ASPIS resource preparation helpers."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FEATURE_SCRIPT = REPO / "workflow/scripts/prepare_feature_set_resources.py"


def write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        return [{key: value or "" for key, value in row.items()} for row in reader]


def assert_contains_row(rows: list[dict[str, str]], label: str, **expected: str) -> None:
    for row in rows:
        if all(row.get(key) == value for key, value in expected.items()):
            return
    raise AssertionError(f"{label} missing row {expected}; observed rows: {rows}")


def main() -> int:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        gtf = write(
            tmp / "toy.gtf",
            'chr1\ttoy\tgene\t1\t100\t.\t+\t.\tgene_id "GENE1"; gene_name "BRCA1"; gene_biotype "protein_coding"; gene_synonym "RNF53|BRCC1";\n'
            'chr1\ttoy\ttranscript\t1\t100\t.\t+\t.\tgene_id "GENE1"; transcript_id "ENST000001"; transcript_name "BRCA1-201";\n',
        )
        custom = write(
            tmp / "sets.tsv",
            "set_id\tfeature_id\tsource\tcollection\tresource_version\tdescription\n"
            "SYMBOL_SET\tBRCA1\tcustom\tsymbols\ttoy\tSymbol-mapped feature\n"
            "TRANSCRIPT_SET\tENST000001\tcustom\ttranscripts\ttoy\tTranscript-mapped feature\n"
            "MISSING_SET\tMISSING\tcustom\tmissing\ttoy\tUnmapped feature\n",
        )
        outdir = tmp / "resource_bundle"
        fragment = outdir / "aspis_feature_sets.yaml"
        command = [
            sys.executable,
            str(FEATURE_SCRIPT),
            "--gtf",
            str(gtf),
            "--outdir",
            str(outdir),
            "--resource-version",
            "toy-release",
            "--custom-table",
            str(custom),
            "--config-fragment",
            str(fragment),
        ]
        result = subprocess.run(command, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if result.returncode != 0:
            raise AssertionError(f"resource preparation failed\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

        normalized = read_tsv(outdir / "custom_sets.tsv")
        assert_contains_row(normalized, "normalized custom table", set_id="SYMBOL_SET", feature_id="GENE1")
        assert_contains_row(normalized, "normalized custom table", set_id="TRANSCRIPT_SET", feature_id="GENE1")
        if any(row.get("set_id") == "MISSING_SET" for row in normalized):
            raise AssertionError("unmapped feature should have been dropped from normalized output")

        unmapped = read_tsv(outdir / "unmapped_features.tsv")
        assert_contains_row(unmapped, "unmapped table", set_id="MISSING_SET", original_feature_id="MISSING")

        gene_map = read_tsv(outdir / "gene_id_map.tsv")
        assert_contains_row(gene_map, "gene map", gene_id="GENE1", gene_name="BRCA1")

        identifier_map = read_tsv(outdir / "gene_identifier_map.tsv")
        assert_contains_row(identifier_map, "identifier map", alias_id="BRCA1", alias_type="gene_symbol", gene_id="GENE1")
        assert_contains_row(identifier_map, "identifier map", alias_id="ENST000001", alias_type="transcript_id", gene_id="GENE1")

        transcript_map = read_tsv(outdir / "transcript_to_gene_map.tsv")
        assert_contains_row(transcript_map, "transcript map", transcript_id="ENST000001", gene_id="GENE1")

        provenance = read_tsv(outdir / "resource_provenance.tsv")
        assert_contains_row(provenance, "provenance", resource_id="gene_identifier_map")
        assert_contains_row(provenance, "provenance", resource_id="transcript_to_gene_map")
        assert_contains_row(provenance, "provenance", resource_id="custom_sets")

        summary = read_tsv(outdir / "resource_summary.tsv")
        assert_contains_row(summary, "summary", resource_id="custom_sets", n_unmapped_or_ambiguous="1")

        fragment_text = fragment.read_text(encoding="utf-8")
        if "report_feature_set_tables" not in fragment_text or "gene_identifier_map" not in fragment_text:
            raise AssertionError(f"config fragment misses expected entries:\n{fragment_text}")

    print("resource preparation contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())