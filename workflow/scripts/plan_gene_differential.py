#!/usr/bin/env python3
"""Compatibility entry point for gene-level DESeq2 contrast planning."""

from __future__ import annotations

import argparse

from plan_feature_differential import add_common_arguments, add_feature_arguments, run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(
        parser,
        counts_flag="--gene-counts",
        counts_help="featureCounts gene_counts.tsv",
    )
    add_feature_arguments(
        parser,
        default_level="gene",
        default_feature_id_column="Geneid",
        default_count_metadata_columns=["Geneid", "Chr", "Start", "End", "Strand", "Length"],
        default_matrix_label="Count matrix",
    )
    return parser.parse_args()


def main() -> int:
    return run(parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
