#!/usr/bin/env python3
"""Validate report-facing feature display label fallbacks."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path("workflow/scripts").resolve()))

from display_labels import feature_display_label, gene_display_label  # noqa: E402
from render_rnaseq_differential_summary import top_feature_table  # noqa: E402


def assert_contains(text: str, expected: str) -> None:
    if expected not in text:
        raise AssertionError(f"Expected {expected!r} in text:\n{text}")


def assert_not_contains(text: str, unexpected: str) -> None:
    if unexpected in text:
        raise AssertionError(f"Unexpected {unexpected!r} in text:\n{text}")


def main() -> int:
    assert gene_display_label("ENSG00000141510", "TP53") == "TP53 (ENSG00000141510)"
    assert gene_display_label("MSTRG.11769", "NA") == "MSTRG.11769"
    assert feature_display_label(
        {
            "Geneid": "ENSG00000141510",
            "gene_name": "TP53",
            "log2FoldChange": "2.5",
            "padj": "1e-8",
        },
        "Geneid",
    ) == "TP53 (ENSG00000141510)"
    assert feature_display_label(
        {
            "transcript_id": "ENST00000269305",
            "gene_id": "ENSG00000141510",
            "gene_name": "TP53",
        },
        "transcript_id",
    ) == "TP53 (ENSG00000141510) | ENST00000269305"
    assert feature_display_label(
        {
            "feature_display": "TP53 (ENSG00000141510)",
            "Geneid": "ENSG00000141510",
            "gene_name": "WRONG",
        },
        "Geneid",
    ) == "TP53 (ENSG00000141510)"

    html = top_feature_table(
        [
            {
                "Geneid": "ENSG00000141510",
                "gene_name": "TP53",
                "log2FoldChange": "2.5",
                "padj": "1e-8",
            },
            {
                "Geneid": "MSTRG.11769",
                "gene_name": "NA",
                "log2FoldChange": "-1.3",
                "padj": "0.02",
            },
        ],
        top_n=10,
    )
    assert_contains(html, "TP53 (ENSG00000141510)")
    assert_contains(html, "MSTRG.11769")
    assert_not_contains(html, "NA (MSTRG.11769)")
    print("report display label contract ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
