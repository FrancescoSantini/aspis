# Independent public biological demonstration dataset

The selected dataset is the ENCODE RNA Evaluation Study comparison of human
GM12878 and K562 cells. It is independent of the datasets used to develop
ASPIS. The subset contains two independent cell-growth biospecimens per cell
line, with matched long total RNA-seq and small RNA-seq from each biospecimen.

| Cell line | Biospecimens | total RNA-seq | small RNA-seq |
| --- | --- | --- | --- |
| GM12878 | ENCBS089RNA, ENCBS090RNA | ENCSR000AEC | ENCSR000AEJ |
| K562 | ENCBS087RNA, ENCBS088RNA | ENCSR000AEL | ENCSR000AES |

The four experiments use rRNA-depleted total RNA, the same Gingeras laboratory
and Illumina HiSeq 2000 platform. Long-RNA libraries are paired-end 101 nt and
small-RNA libraries are single-end 101 nt. The public files total approximately
77.4 GiB compressed.

This is a cell-line contrast, not a treatment experiment. The intake therefore
uses `condition=GM12878|K562`, `treatment=none`, and matched ENCODE biospecimen
accessions. The small-RNA libraries are broad `<200 nt` libraries pretreated
with tobacco acid pyrophosphatase; they must not be described as miRNA-only
libraries.

Download and verify the raw FASTQs:

```bash
bash examples/public_encode_demo/download_fastqs.sh
```

The script supports interrupted transfers through `curl --continue-at -` and
verifies the ENCODE MD5 checksums. `encode_fastq_manifest.tsv` preserves the
source accessions, URLs, sizes, layouts and checksums. `intake.tsv` is ready for
ASPIS after download.

Before analysis, provide a consistent GRCh38 FASTA/GTF, human mature-miRNA
FASTA, contaminant reference and corresponding indexes in a project config.
The adapter setting should be taken from the ENCODE library metadata or
confirmed from the FASTQ sequence content; it is intentionally not guessed in
this example.

Primary records:

- https://www.encodeproject.org/experiments/ENCSR000AEC/
- https://www.encodeproject.org/experiments/ENCSR000AEJ/
- https://www.encodeproject.org/experiments/ENCSR000AEL/
- https://www.encodeproject.org/experiments/ENCSR000AES/
