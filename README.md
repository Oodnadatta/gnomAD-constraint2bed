# gnomAD-constraint2bed

Author: **Anne-Sophie Denommé-Pichon**

Build BED files of constrained genes (missense Z score and LOEUF) from the **gnomAD v4.1.2** gene constraint table (GRCh38), using all transcripts:

- genes intolerant to **missense** variants (`mis.z_score` above a threshold);
- genes intolerant to **loss-of-function** variants (`LOEUF` below a threshold).

Default thresholds: **2** for `mis.z_score` and **0.45** for `LOEUF`. Both can be changed
with command-line options (see [Usage](#usage)).

## Requirements

Python 3.9 or later. No external package is needed (standard library only).
The code follows PEP 8 and PEP 257, and is type-annotated (PEP 484); it is checked with
flake8, pydocstyle and `mypy --strict`.

## Input

The gnomAD v4.1.2 constraint metrics table:

```
https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.2/constraint/gnomad.v4.1.2.constraint_metrics.tsv.bgz
```

The script reads the compressed `.tsv.bgz` file or the decompressed `.tsv` file.
The table is about 1.5 GB once decompressed; reading it takes a few minutes.

## Usage

```bash
python3 gnomad_constraint_to_bed_asdp.py gnomad.v4.1.2.constraint_metrics.tsv.bgz
python3 gnomad_constraint_to_bed_asdp.py gnomad.v4.1.2.constraint_metrics.tsv.bgz \
    --missense-z 4 --loeuf 0.45 --outdir results
```

| Option | Default | Meaning |
|---|---|---|
| `--missense-z` | `2` | keep genes with `mis.z_score` **above** this value |
| `--loeuf` | `0.45` | keep genes with `LOEUF` (`lof.oe_ci.upper`) **below** this value |
| `--outdir` | current directory | output directory |

The default LOEUF threshold of 0.45 is the one recommended by gnomAD for v4
(it replaces the 0.35 threshold used with gnomAD v2).

## Output

Two BED files, named after the thresholds used:

```
gnomad_v4.1.2_GRCh38_misZ_gt<threshold>.bed
gnomad_v4.1.2_GRCh38_LOEUF_lt<threshold>.bed
```

Five tab-separated columns, no header:

| Column | Content |
|---|---|
| 1 | chromosome (`chr1`, …, `chrX`, `chrY`) |
| 2 | start (0-based, BED convention) |
| 3 | end |
| 4 | gene symbol (Ensembl gene ID when the gene has no symbol) |
| 5 | best score among the transcripts kept: highest `mis.z_score`, or lowest `LOEUF` |

The files are not sorted. To sort them, for tools that require it:

```bash
sort -k1,1 -k2,2n file.bed > file.sorted.bed
```

A summary is printed at the end: number of Ensembl transcripts read and excluded (per file),
number of readthrough genes split, number of genes in each file.

## Method

- **Ensembl transcripts only.** The table holds Ensembl (`ENST…`) and RefSeq transcripts,
  but only the Ensembl transcripts have coordinates.
- **All transcripts, not only MANE Select.** A gene is kept as soon as **at least one**
  of the gene's transcripts passes the threshold, so that exons present only in
  alternative transcripts (including MANE Plus Clinical) are not missed.
- **Whole transcripts.** For each gene, only the transcripts that pass the threshold are used.
  The region starts at the first base of these transcripts and ends at the last base of
  these transcripts: UTRs, exons **and introns** are included.
- **Genes without a symbol** are named by their Ensembl gene ID (`ENSG…`).
- **Genes on chrX and chrY** (pseudoautosomal regions) get one line per chromosome.
- **Excluded transcripts.** Transcripts flagged in `gene_flags` are excluded from both files:
  - `low_exome_mapping_quality`: in these regions (mostly segmental duplications) reads map
    poorly, gnomAD observes fewer variants than really exist, and genes look constrained
    artificially;
  - `low_exome_coverage`: too few exomes are well covered, so the observed number of variants
    is underestimated in the same way.

  Transcripts flagged in `constraint_flags` are excluded from one file only:
  - `outlier_mis`: excluded from the missense file (the observed number of missense variants
    is an outlier compared with the expected number);
  - `outlier_lof`: excluded from the LOEUF file (same, for loss-of-function variants).

  The excluded flags are listed at the top of the script (`EXCLUDED_FLAGS`,
  `EXCLUDED_FLAGS_MISSENSE`, `EXCLUDED_FLAGS_LOEUF`).
- **Readthrough transcripts.** A readthrough transcript spanning two neighbouring genes
  (e.g. `ANKHD1-EIF4EBP3`) is assigned to each of the two genes (`ANKHD1` and `EIF4EBP3`),
  so that every name in the BED file is a gene symbol. A symbol is split only if each part
  is itself a gene symbol of the table, so that symbols such as `NKX2-1` or `H3-3A` are kept.

## Limitations

- For a gene split from a readthrough transcript, the region covers the whole readthrough
  transcript and therefore overlaps the neighbouring gene.
- The score in column 5 comes from the most constrained transcript of the gene, which is not
  necessarily the MANE Select transcript shown first on the gnomAD website. This is intended:
  all transcripts are considered so that another transcript of interest, constrained even
  when the MANE Select transcript is not, is not discarded.
- Transcripts without a LOEUF value (no loss-of-function variant expected, e.g. very short
  transcripts) cannot appear in the LOEUF file.
- Genes are grouped by symbol: two distinct Ensembl genes sharing the same symbol on the
  same chromosome would be merged into one region.
- Column and file names are specific to gnomAD **v4.1.2**.

## GRCh37

The output is in **GRCh38**. If you need GRCh37 coordinates, the BED files can be converted
with a liftover tool, such as UCSC liftOver (https://genome.ucsc.edu/cgi-bin/hgLiftOver)
or CrossMap. Some regions may fail to convert.

## Citation

If you use this script, please cite it (see `CITATION.cff`) and cite gnomAD as
recommended on https://gnomad.broadinstitute.org.

## License

MIT License: see `LICENSE`. You may use, modify and share this script, provided that the
copyright notice with the author's name is kept.
