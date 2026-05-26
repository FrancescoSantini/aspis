#!/usr/bin/env python3
"""Compatibility entry point for transcript-level DESeq2 contrast planning."""

from __future__ import annotations

import argparse

from plan_feature_differential import add_common_arguments, add_feature_arguments, run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(
        parser,
        counts_flag="--transcript-counts",
        counts_help="StringTie transcript count matrix",
    )
    add_feature_arguments(
        parser,
        default_level="transcript",
        default_feature_id_column="transcript_id",
        default_count_metadata_columns=["transcript_id"],
        default_matrix_label="Transcript count matrix",
    )
    return parser.parse_args()


def main() -> int:
    return run(parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
