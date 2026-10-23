# ASPIS deterministic synthetic project

This fixture is a redistribution-safe technical test of a matched RNA-seq and
smallRNA-seq project. It contains two conditions, two matched biospecimens per
condition, and both assays for every biospecimen. The inserted counts encode
opposite shifts for `syn_gene_a`/`syn_gene_b` and
`syn-miR-1-5p`/`syn-miR-2-3p`. SmallRNA libraries also contain fixed contaminant
and residual-genome reads.

The fixture is not a biological benchmark. Its purpose is to test installation,
routing, paired-assay integration, and the propagation of known signals through
the workflow.

Regenerate and validate the committed data:

```bash
python3 examples/synthetic_project/generate_fixture.py
python3 tests/validate_synthetic_project.py
```

Run the project after installing the full ASPIS environment:

```bash
snakemake results/synthetic_project/index.done \
  --configfile examples/synthetic_project/aspis.synthetic.yaml \
  --cores 4 \
  --rerun-incomplete
```

`truth_rnaseq_counts.tsv` and `truth_smallrna_counts.tsv` describe the intended
input signals. `checksums.sha256` records every generated reference and FASTQ.

For a clean local benchmark, first remove only this fixture's generated run
namespaces, regenerate the fixture, and use the benchmark helper:

```bash
rm -rf results/synthetic_project work/synthetic_project meta/synthetic_project
python3 examples/synthetic_project/generate_fixture.py
python3 workflow/scripts/benchmark_snakemake.py \
  --target results/synthetic_project/index.done \
  --configfile examples/synthetic_project/aspis.synthetic.yaml \
  --cores 4 \
  --label synthetic_project_local \
  --mode fresh \
  --rerun-incomplete
```

See `docs/benchmarking.md` for interpretation and the separate HPC procedure.
