#!/usr/bin/env python3
"""Run planned gene-level DESeq2 contrasts for an RNA-seq branch."""

from __future__ import annotations

import argparse

from run_feature_differential_branch import FeatureRunSpec, add_common_arguments, run


SPEC = FeatureRunSpec(
    label="gene",
    counts_attr="gene_counts",
    metadata_attr="gene_metadata",
    feature_id_column="Geneid",
    count_matrix_label="Count matrix",
    no_rows_message="Gene differential plan has no rows",
    failed_message="DESeq2 failed for contrast(s)",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Gene differential contrast plan TSV")
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--gene-counts", required=True, help="featureCounts gene_counts.tsv")
    parser.add_argument("--gene-metadata", required=True, help="featureCounts gene_metadata.tsv")
    add_common_arguments(parser, "Minimum total count per gene")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(run(parse_args(), SPEC))
