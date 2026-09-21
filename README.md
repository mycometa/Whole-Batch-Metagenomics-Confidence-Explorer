# Whole-Batch-Metagenomics-Confidence-Explorer
# A sibling to the Fungal Panel Batch Confidence Explorer, but scoped to the
# ENTIRE wf-metagenomics counts table — every organism, every kingdom, not
# just a curated fungal reporting panel. Works from just a MinKNOW HTML
# report and the counts CSV; an optional position-mapping spreadsheet adds
# barcode-level read/QC context and DNA-quantity/PCR commentary on top.
#
# INPUT FILES
#   1. MinKNOW run report HTML (required)
#   2. wf-metagenomics counts CSV (required) — species x sample matrix, with
#      taxonomy rank columns (superkingdom..genus) and a "total" column, and
#      a sample column whose name contains "control" for the negative control.
#   3. Position-mapping spreadsheet .xlsx (optional) — columns Position,
#      Swab ID, and (if present) Pre-PCR (ng/uL) / Post-PCR (ng/uL). Rows
#      with a blank Position are treated as auxiliary/lot-marker rows and
#      skipped. Providing this file lets the app join the HTML's per-barcode
#      read/QC stats onto each sample, and enables PCR efficiency commentary.
#      Without it, the app still works fully, using classified-read counts
#      (computed directly from the CSV) as the depth metric throughout, and
#      samples are identified by their CSV column (swab ID) rather than a
#      barcode number.
#
# KEY ASSUMPTIONS
#   - If a position-mapping file is provided, barcode number = pooling
#     position number, in order (barcode01 = position 1, etc.) — same
#     convention used throughout this app family.
#   - "Notable" calls (shown in search, verdicts, and the export) are those
#     meeting BOTH a configurable minimum read count AND minimum % of that
#     sample's classified reads, to keep the full ~2,000-species table
#     tractable — adjustable in the sidebar.
#   - The negative control column is identified by column name containing
#     "control" (case-insensitive).
#
# REQUIRED PACKAGES
#   install.packages(c("shiny","bslib","DT","ggplot2","dplyr","tidyr",
#                       "stringr","readxl","scales","officer"))
#
# RUN
#   shiny::runApp("app.R")
