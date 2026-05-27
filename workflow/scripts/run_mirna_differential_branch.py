#!/usr/bin/env python3
"""Run planned miRNA-level DESeq2 contrasts for a smallRNA branch."""

from __future__ import annotations

import argparse

from run_feature_differential_branch import FeatureRunSpec, add_common_arguments, run


SPEC = FeatureRunSpec(
    label="mirna",
    counts_attr="mirna_counts",
    metadata_attr="mirna_metadata",
    feature_id_column="Geneid",
    count_matrix_label="miRNA count matrix",
    no_rows_message="miRNA differential plan has no rows",
    failed_message="miRNA DESeq2 failed for contrast(s)",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="miRNA differential contrast plan TSV")
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--mirna-counts", required=True, help="smallRNA miRNA count matrix TSV")
    parser.add_argument("--mirna-metadata", required=True, help="smallRNA miRNA metadata TSV")
    add_common_arguments(parser, "Minimum total count per miRNA")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(run(parse_args(), SPEC))
