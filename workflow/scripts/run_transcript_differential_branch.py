#!/usr/bin/env python3
"""Run planned transcript-level DESeq2 contrasts for an RNA-seq branch."""

from __future__ import annotations

import argparse

from run_feature_differential_branch import FeatureRunSpec, add_common_arguments, run


SPEC = FeatureRunSpec(
    label="transcript",
    counts_attr="transcript_counts",
    metadata_attr="transcript_metadata",
    feature_id_column="transcript_id",
    count_matrix_label="Transcript count matrix",
    no_rows_message="Transcript differential plan has no rows",
    failed_message="Transcript DESeq2 failed for contrast(s)",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, help="Transcript differential contrast plan TSV")
    parser.add_argument("--samples", required=True, help="Branch samples.tsv")
    parser.add_argument("--transcript-counts", required=True, help="StringTie transcript_counts.tsv")
    parser.add_argument("--transcript-metadata", required=True, help="StringTie transcript_metadata.tsv")
    add_common_arguments(parser, "Minimum total count per transcript")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(run(parse_args(), SPEC))
