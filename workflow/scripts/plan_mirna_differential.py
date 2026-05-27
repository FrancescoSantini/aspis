#!/usr/bin/env python3
"""Plan miRNA-level DESeq2 contrasts from smallRNA featureCounts matrices."""

from __future__ import annotations

import argparse

from plan_feature_differential import add_common_arguments, add_feature_arguments, run


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(
        parser,
        counts_flag="--mirna-counts",
        counts_help="smallRNA miRNA count matrix TSV",
    )
    add_feature_arguments(
        parser,
        default_level="mirna",
        default_feature_id_column="Geneid",
        default_count_metadata_columns=["Geneid", "Chr", "Start", "End", "Strand", "Length"],
        default_matrix_label="miRNA count matrix",
    )
    parser.set_defaults(assay="smallrna")
    return parser.parse_args()


def main() -> int:
    return run(parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
