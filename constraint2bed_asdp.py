#!/usr/bin/env python3
# Copyright (c) 2026 Anne-Sophie Denommé-Pichon
# License: MIT (see LICENSE)
# If you use this script, please cite it (see CITATION.cff in the repository).
r"""Build BED files of constrained genes from the gnomAD v4.1.2 table.

The gnomAD v4.1.2 constraint metrics table (GRCh38) gives two BED files:

  gnomad_v4.1.2_GRCh38_misZ_gt<threshold>.bed
      genes with at least one Ensembl transcript with mis.z_score > threshold
  gnomad_v4.1.2_GRCh38_LOEUF_lt<threshold>.bed
      genes with at least one Ensembl transcript with LOEUF < threshold

Default thresholds: 2 for mis.z_score and 0.45 for LOEUF. They can be
changed with the --missense-z and --loeuf options (see Usage below).

Rules:
  - Only Ensembl transcripts (identifiers starting with ENST) are used:
    they are the only transcripts with coordinates.
  - A gene is kept if at least one of the gene's transcripts passes the
    threshold.
  - For each gene, only the transcripts that pass the threshold are used.
    The region starts at the first base of these transcripts and ends at
    the last base of these transcripts (UTRs, exons and introns included).
  - Genes without a symbol ("NA") are named by their Ensembl gene ID.
  - Genes on both chrX and chrY get one line per chromosome.
  - Transcripts flagged "low_exome_mapping_quality" or "low_exome_coverage"
    (gene_flags) are excluded from both files. Transcripts flagged
    "outlier_mis" (constraint_flags) are excluded from the missense file
    only, and transcripts flagged "outlier_lof" from the LOEUF file only.
  - Readthrough transcripts (e.g. ANKHD1-EIF4EBP3) are assigned to each of
    the two genes (ANKHD1 and EIF4EBP3), so that every name in the BED
    files is a gene symbol.

Columns: chromosome, start (0-based), end, gene, best score among the
transcripts kept.

Usage:
  python3 constraint2bed_asdp.py TABLE
  python3 constraint2bed_asdp.py TABLE \
      --missense-z 4 --loeuf 0.45 --outdir results

TABLE is gnomad.v4.1.2.constraint_metrics.tsv, plain or compressed (.bgz).

Author: Anne-Sophie Denommé-Pichon. License: MIT.
If you use this script, please cite it (see CITATION.cff).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections.abc import Callable
from pathlib import Path
from typing import TextIO

AUTHOR = "Anne-Sophie Denommé-Pichon"

GENE_COLUMN = "gene"
TRANSCRIPT_COLUMN = "transcript"
CHROMOSOME_COLUMN = "chromosome"
START_COLUMN = "start_position"
END_COLUMN = "end_position"
MISSENSE_Z_COLUMN = "mis.z_score"
LOEUF_COLUMN = "lof.oe_ci.upper"
GENE_ID_COLUMN = "gene_id"
GENE_FLAGS_COLUMN = "gene_flags"
CONSTRAINT_FLAGS_COLUMN = "constraint_flags"

# Types used in the annotations
Flags = dict[str, list[str]]
# Transcript passing a threshold: (gene, chromosome, start, end, score)
Transcript = tuple[str, str, int, int, float]
# Gene regions, by (gene, chromosome), for example:
# {("GENE_X", "chr1"): {"start": 800, "end": 5000, "best_score": 0.3}}
Regions = dict[tuple[str, str], dict[str, float]]
ScoreComparison = Callable[[float, float], bool]

# Transcripts with one of these flags are excluded from both BED files
EXCLUDED_FLAGS: Flags = {
    GENE_FLAGS_COLUMN: ["low_exome_mapping_quality", "low_exome_coverage"],
    CONSTRAINT_FLAGS_COLUMN: [
        # "outlier_syn",
    ],
}

# Transcripts with one of these flags are excluded from the missense BED only
EXCLUDED_FLAGS_MISSENSE: Flags = {
    CONSTRAINT_FLAGS_COLUMN: ["outlier_mis"],
}

# Transcripts with one of these flags are excluded from the LOEUF BED only
EXCLUDED_FLAGS_LOEUF: Flags = {
    CONSTRAINT_FLAGS_COLUMN: ["outlier_lof"],
}


def parse_arguments() -> argparse.Namespace:
    """Read the command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "table", type=Path,
        help="gnomAD v4.1.2 constraint metrics table (.tsv or .tsv.bgz)",
    )
    parser.add_argument(
        "--missense-z", type=float, default=2.0,
        help="keep genes with mis.z_score above this value (default: 2)",
    )
    parser.add_argument(
        "--loeuf", type=float, default=0.45,
        help="keep genes with LOEUF below this value (default: 0.45)",
    )
    parser.add_argument(
        "--outdir", type=Path, default=Path("."),
        help="output directory (default: current directory)",
    )
    return parser.parse_args()


def open_table(path: Path) -> TextIO:
    """Open the table as text, whether it is compressed or not."""
    with open(path, "rb") as input_file:  # "rb" = read binary
        # gzip/bgzip signature
        is_compressed = input_file.read(2) == b"\x1f\x8b"
    if is_compressed:
        # "rt" = read text
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, encoding="utf-8", newline="")


def parse_float(value: str) -> float | None:
    """Return the value as a number, or None if it is not a number (NA)."""
    try:
        number = float(value)
    except ValueError:
        return None
    if number != number:  # NaN
        return None
    return number


def has_excluded_flag(row: dict[str, str], excluded_flags: Flags) -> bool:
    """Tell whether the transcript carries one of the excluded flags."""
    for column, flags in excluded_flags.items():
        for flag in flags:
            if flag in row[column]:
                return True
    return False


def higher_is_better(new_score: float, old_score: float) -> bool:
    """Compare missense Z scores: the highest value is kept."""
    return new_score > old_score


def lower_is_better(new_score: float, old_score: float) -> bool:
    """Compare LOEUF values: the lowest value is kept."""
    return new_score < old_score


def update_region(
    regions: Regions,
    gene: str,
    chromosome: str,
    start: int,
    end: int,
    score: float,
    is_better: ScoreComparison,
) -> None:
    """Extend the gene region with this transcript; keep the best score."""
    gene_chr = (gene, chromosome)
    if gene_chr not in regions:
        regions[gene_chr] = {"start": start, "end": end, "best_score": score}
        return
    region = regions[gene_chr]
    region["start"] = min(region["start"], start)
    region["end"] = max(region["end"], end)
    if is_better(score, region["best_score"]):
        region["best_score"] = score


def readthrough_parts(gene: str, known_symbols: set[str]) -> list[str]:
    """Split a readthrough symbol (e.g. "ANKHD1-EIF4EBP3") into genes.

    The symbol is split only if every part is itself a gene symbol of the
    table, so that symbols such as NKX2-1 or H3-3A are left untouched.
    Return the list of genes, or [gene] if the symbol is not split.
    """
    if "-" not in gene:
        return [gene]
    parts = gene.split("-")
    if all(part in known_symbols for part in parts):
        return parts
    return [gene]


def write_bed(regions: Regions, path: Path) -> int:
    """Write one BED line per gene and chromosome; return the line count."""
    with open(path, "w", encoding="utf-8") as output:
        for (gene, chromosome), region in regions.items():
            fields = [
                chromosome,
                str(region["start"]),
                str(region["end"]),
                gene,
                f"{region['best_score']:.4g}",
            ]
            output.write("\t".join(fields) + "\n")
    return len(regions)


def main() -> None:
    """Read the gnomAD table and write the two BED files."""
    arguments = parse_arguments()
    arguments.outdir.mkdir(parents=True, exist_ok=True)

    missense_transcripts: list[Transcript] = []
    loeuf_transcripts: list[Transcript] = []
    known_symbols: set[str] = set()
    ensembl_transcripts = 0
    excluded_transcripts = 0
    excluded_missense_only = 0
    excluded_loeuf_only = 0

    with open_table(arguments.table) as table:
        reader = csv.DictReader(table, delimiter="\t", quoting=csv.QUOTE_NONE)

        required = [
            GENE_COLUMN, TRANSCRIPT_COLUMN, CHROMOSOME_COLUMN, START_COLUMN,
            END_COLUMN, MISSENSE_Z_COLUMN, LOEUF_COLUMN, GENE_ID_COLUMN,
            GENE_FLAGS_COLUMN, CONSTRAINT_FLAGS_COLUMN,
        ]
        columns = reader.fieldnames or []  # None if the table is empty
        missing = [name for name in required if name not in columns]
        if missing:
            sys.exit(f"Missing columns: {', '.join(missing)}")

        for row in reader:
            known_symbols.add(row[GENE_COLUMN])

            # Ensembl transcripts only (RefSeq ones have no coordinates)
            if not row[TRANSCRIPT_COLUMN].startswith("ENST"):
                continue
            if row[CHROMOSOME_COLUMN] == "NA":
                continue
            ensembl_transcripts += 1

            if has_excluded_flag(row, EXCLUDED_FLAGS):
                excluded_transcripts += 1
                continue

            gene = row[GENE_COLUMN]
            if gene == "NA":  # no gene symbol: use the Ensembl gene ID
                gene = row[GENE_ID_COLUMN]
            chromosome = row[CHROMOSOME_COLUMN]
            # gnomAD is 1-based, BED is 0-based
            start = int(row[START_COLUMN]) - 1
            end = int(row[END_COLUMN])

            if has_excluded_flag(row, EXCLUDED_FLAGS_MISSENSE):
                excluded_missense_only += 1
            else:
                missense_z = parse_float(row[MISSENSE_Z_COLUMN])
                if (missense_z is not None
                        and missense_z > arguments.missense_z):
                    missense_transcripts.append(
                        (gene, chromosome, start, end, missense_z)
                    )

            if has_excluded_flag(row, EXCLUDED_FLAGS_LOEUF):
                excluded_loeuf_only += 1
            else:
                loeuf = parse_float(row[LOEUF_COLUMN])
                if loeuf is not None and loeuf < arguments.loeuf:
                    loeuf_transcripts.append(
                        (gene, chromosome, start, end, loeuf)
                    )

    # Build the gene regions once the whole table is read: readthrough
    # symbols can only be recognised when all gene symbols are known.
    missense_regions: Regions = {}
    loeuf_regions: Regions = {}
    split_readthroughs: set[str] = set()
    comparisons: list[tuple[list[Transcript], Regions, ScoreComparison]] = [
        (missense_transcripts, missense_regions, higher_is_better),
        (loeuf_transcripts, loeuf_regions, lower_is_better),
    ]
    for transcripts, regions, is_better in comparisons:
        for gene, chromosome, start, end, score in transcripts:
            parts = readthrough_parts(gene, known_symbols)
            if len(parts) > 1:
                split_readthroughs.add(gene)
            for part in parts:
                update_region(
                    regions, part, chromosome, start, end, score, is_better
                )

    # File names are built from the thresholds
    # (":g" writes 2.0 as 2 and 0.45 as 0.45)
    missense_file = (
        f"gnomad_v4.1.2_GRCh38_misZ_gt{arguments.missense_z:g}.bed"
    )
    loeuf_file = f"gnomad_v4.1.2_GRCh38_LOEUF_lt{arguments.loeuf:g}.bed"
    missense_count = write_bed(
        missense_regions, arguments.outdir / missense_file
    )
    loeuf_count = write_bed(loeuf_regions, arguments.outdir / loeuf_file)

    print(f"gnomAD-constraint2bed, by {AUTHOR}"
          " - please cite this script if you use it (see CITATION.cff)")
    print(f"Ensembl transcripts read       : {ensembl_transcripts}")
    print(f"  excluded from both BED files : {excluded_transcripts}")
    print(f"  excluded from missense only  : {excluded_missense_only}")
    print(f"  excluded from LOEUF only     : {excluded_loeuf_only}")
    print(f"Readthrough genes split        : {len(split_readthroughs)}")
    print(f"Genes with mis.z_score > {arguments.missense_z:g}"
          f" : {missense_count} -> {missense_file}")
    print(f"Genes with LOEUF < {arguments.loeuf:g}"
          f" : {loeuf_count} -> {loeuf_file}")


if __name__ == "__main__":
    main()
