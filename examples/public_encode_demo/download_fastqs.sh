#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
manifest="$root/examples/public_encode_demo/encode_fastq_manifest.tsv"
outdir="$root/examples/public_encode_demo/data"
mkdir -p "$outdir"

tail -n +2 "$manifest" |
while IFS= read -r row; do
    condition=$(printf '%s\n' "$row" | cut -f1)
    assay=$(printf '%s\n' "$row" | cut -f3)
    accession=$(printf '%s\n' "$row" | cut -f6)
    size=$(printf '%s\n' "$row" | cut -f9)
    md5=$(printf '%s\n' "$row" | cut -f10)
    url=$(printf '%s\n' "$row" | cut -f11)
    target="$outdir/${accession}.fastq.gz"
    echo "Downloading $accession ($condition, $assay, ${size} bytes)"
    curl --fail --location --retry 5 --retry-delay 10 --continue-at - --output "$target" "$url"
    printf '%s  %s\n' "$md5" "$target" | md5sum --check --status || {
        echo "Checksum failed for $target" >&2
        exit 1
    }
done

echo "All ENCODE FASTQs downloaded and verified in $outdir"
