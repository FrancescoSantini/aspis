# Local Smoke Tests

Run local checks from the repository root inside the `aspis-smk9` environment.
These tests use tiny fixtures under `tests/` and do not require G100.

```bash
conda activate aspis-smk9
cd /mnt/d/wdir/aspis
```

Use a dry-run first when changing rule targets or config contracts:

```bash
MODE=dry-run bash tests/run_local_smokes.sh
```

Run the technical smoke suite before opening a PR:

```bash
bash tests/run_local_smokes.sh
```

The suite covers the default materialization/QC path, HISAT2 alignment, STAR
alignment, RNA-seq quantification, differential-layer planning, and a
DESeq2-only differential test. It is a local confidence gate only; it does not replace a deliberate G100 SLURM smoke
run after local contracts are stable.
