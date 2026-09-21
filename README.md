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
library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(scales)
library(officer)

options(shiny.maxRequestSize = 300 * 1024^2)

# -----------------------------------------------------------------------------
# Parsing helpers
# -----------------------------------------------------------------------------

parse_html_barcode_stats <- function(path) {
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  pattern <- paste0(
    '\\{"barcode":\\s*"barcode(\\d+)",\\s*"total_bases":\\s*(\\d+),\\s*',
    '"passed_bases_percent":\\s*([\\d.]+),\\s*"total_reads":\\s*(\\d+),\\s*',
    '"passed_reads_percent":\\s*([\\d.]+)\\}'
  )
  m <- str_match_all(txt, pattern)[[1]]
  if (nrow(m) == 0) stop("Could not find per-barcode statistics in this HTML report.")
  tibble(
    barcode = str_pad(m[, 2], 2, pad = "0"),
    total_bases = as.numeric(m[, 3]),
    passed_bases_percent = as.numeric(m[, 4]),
    total_reads = as.numeric(m[, 5]),
    passed_reads_percent = as.numeric(m[, 6])
  ) %>% distinct(barcode, .keep_all = TRUE) %>% arrange(barcode)
}

# Extract the full reportData JSON object embedded in a MinKNOW HTML report,
# for the exported report's Run Information / Run Configuration sections.
parse_html_run_info <- function(path) {
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  marker <- "const reportData="
  pos <- str_locate(txt, fixed(marker))
  if (is.na(pos[1, 1])) return(NULL)
  remainder <- substring(txt, pos[1, "end"] + 1)
  chars_vec <- strsplit(remainder, "")[[1]]
  depth_vec <- cumsum(ifelse(chars_vec == "{", 1L, ifelse(chars_vec == "}", -1L, 0L)))
  end_offset <- which(depth_vec == 0)[1]
  if (is.na(end_offset)) return(NULL)
  blob <- paste(chars_vec[1:end_offset], collapse = "")
  tryCatch(jsonlite::fromJSON(blob, simplifyVector = FALSE), error = function(e) NULL)
}

titlevalue_to_df <- function(lst) {
  if (is.null(lst) || length(lst) == 0) return(data.frame(Setting = character(0), Value = character(0)))
  data.frame(
    Setting = sapply(lst, function(x) as.character(x$title)),
    Value = sapply(lst, function(x) as.character(x$value)),
    stringsAsFactors = FALSE
  )
}

format_iso_time <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return("\u2014")
  str_replace(str_replace(x, "T", " "), "\\.\\d+Z$", " UTC")
}

# Optional position-mapping spreadsheet: Position, Swab ID, and (if present)
# Pre-PCR / Post-PCR columns. Rows with a blank Position are auxiliary rows
# (lot codes etc.) and are skipped.
parse_position_map <- function(path) {
  raw <- read_excel(path, col_names = TRUE)
  names(raw) <- str_squish(names(raw))
  pos_col  <- names(raw)[str_detect(names(raw), regex("^Position$", ignore_case = TRUE))][1]
  swab_col <- names(raw)[str_detect(names(raw), regex("Swab", ignore_case = TRUE))][1]
  pre_col  <- names(raw)[str_detect(names(raw), regex("Pre-?PCR", ignore_case = TRUE))][1]
  post_col <- names(raw)[str_detect(names(raw), regex("Post-?PCR", ignore_case = TRUE))][1]
  if (is.na(pos_col) || is.na(swab_col)) {
    stop("Could not identify Position / Swab ID columns in the position-mapping spreadsheet.")
  }

  # Row-by-row so an optional GTX case code can be picked up: rows with a
  # blank Position are auxiliary (e.g. a kit-lot or case-code marker row
  # immediately below the sample it belongs to) rather than genuine samples.
  # If the spreadsheet has no such rows, "case" simply stays NA for every
  # sample and case-grouping features are skipped gracefully.
  n <- nrow(raw)
  samples <- list()
  current <- NULL
  for (i in seq_len(n)) {
    pos_val <- raw[[pos_col]][i]
    if (!is.na(suppressWarnings(as.numeric(pos_val)))) {
      if (!is.null(current)) samples[[length(samples) + 1]] <- current
      current <- list(
        position = suppressWarnings(as.numeric(pos_val)),
        swab_id_raw = as.character(raw[[swab_col]][i]),
        pre_pcr = if (!is.na(pre_col)) suppressWarnings(as.numeric(raw[[pre_col]][i])) else NA_real_,
        post_pcr = if (!is.na(post_col)) suppressWarnings(as.numeric(raw[[post_col]][i])) else NA_real_,
        case = NA_character_
      )
    } else if (!is.null(current)) {
      swab_val <- as.character(raw[[swab_col]][i])
      if (!is.na(swab_val) && nzchar(swab_val) && is.na(current$case)) current$case <- swab_val
    }
  }
  if (!is.null(current)) samples[[length(samples) + 1]] <- current

  bind_rows(lapply(samples, as_tibble)) %>%
    arrange(position) %>%
    mutate(barcode = str_pad(as.character(position), 2, pad = "0"))
}

match_csv_column <- function(swab_id_raw, csv_colnames) {
  if (is.na(swab_id_raw)) return(NA_character_)
  if (str_detect(swab_id_raw, regex("control", ignore_case = TRUE))) {
    hit <- csv_colnames[str_detect(csv_colnames, regex("control", ignore_case = TRUE))]
    return(if (length(hit) >= 1) hit[1] else NA_character_)
  }
  digits <- str_extract(swab_id_raw, "\\d+")
  if (is.na(digits)) return(NA_character_)
  padded <- str_pad(digits, 9, pad = "0")
  if (padded %in% csv_colnames) return(padded)
  numeric_match <- csv_colnames[
    suppressWarnings(as.numeric(csv_colnames)) == suppressWarnings(as.numeric(digits))
  ]
  if (length(numeric_match) >= 1) return(numeric_match[1])
  NA_character_
}

parse_counts_csv <- function(path) {
  raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(raw)[1] <- "species"
  taxonomy_ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus","tax")
  non_sample_cols <- c("species", "total", taxonomy_ranks)
  sample_cols <- setdiff(names(raw), non_sample_cols)
  sample_cols <- sample_cols[sapply(raw[sample_cols], function(x) {
    suppressWarnings(!all(is.na(as.numeric(as.character(x)))))
  })]
  taxonomy_present <- intersect(taxonomy_ranks, names(raw))
  long <- raw %>%
    select(all_of(c("species", taxonomy_present, sample_cols))) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample_col", values_to = "reads") %>%
    mutate(reads = suppressWarnings(as.numeric(reads)), reads = ifelse(is.na(reads), 0, reads))
  control_col <- sample_cols[str_detect(sample_cols, regex("control", ignore_case = TRUE))]
  control_col <- if (length(control_col) >= 1) control_col[1] else NA_character_
  list(long = long, sample_cols = sample_cols, control_col = control_col, taxonomy_present = taxonomy_present)
}

pcr_comment <- function(pre_pcr, post_pcr, batch_pre = NULL, batch_post = NULL) {
  pre_val <- if (is.na(pre_pcr)) 0 else pre_pcr
  post_val <- if (is.na(post_pcr)) 0 else post_pcr
  pre_txt <- if (is.na(pre_pcr) || pre_pcr == 0) "undetected" else sprintf("%.3g ng/\u00b5L", pre_pcr)
  post_txt <- if (is.na(post_pcr) || post_pcr == 0) "undetected" else sprintf("%.3g ng/\u00b5L", post_pcr)
  pre_detectable <- pre_val > 0.2

  # Batch-relative checks: how does this sample compare to the REST of the
  # batch, rather than to one fixed number that may not fit every batch/kit
  # lot? Two distinct technical failure modes are distinguished here, both
  # from a low or underperforming yield, but with different implications:
  #   - PCR inhibition: DNA was clearly extracted (high pre-PCR), but it
  #     failed to amplify as well as similar/lower-input samples did.
  #   - Poor extraction: most of the batch yielded decent pre-PCR DNA, but
  #     THIS sample's is an outlier on the low end \u2014 the shortfall is
  #     upstream of PCR, in extraction, not a PCR chemistry problem.
  # Neither is the same as ordinary low biomass, where LOW pre-PCR is the
  # norm for most of the batch (a real feature of the sample type), not an
  # outlier for this one sample specifically.
  if (!is.null(batch_pre) && !is.null(batch_post)) {
    valid <- !is.na(batch_pre) & !is.na(batch_post)
    bp <- batch_pre[valid]; bq <- batch_post[valid]
    if (length(bp) >= 4) {
      pre_pct <- mean(bp <= pre_val)
      post_pct <- mean(bq <= post_val)
      frac_higher_pre <- mean(bp > pre_val)
      batch_pre_median <- median(bp)

      if (pre_detectable && pre_pct >= 0.75 && post_pct <= 0.5) {
        verdict <- sprintf(
          "this sample's pre-PCR concentration ranks among the highest in the batch (~%.0fth percentile), yet its post-PCR yield ranks only ~%.0fth percentile \u2014 well below what similarly- or lower-input samples achieved here. Strong evidence of PCR inhibition specific to this sample; consider repeating PCR on a diluted aliquot",
          100 * pre_pct, 100 * post_pct
        )
        return(sprintf("PCR: Pre-PCR %s \u2192 Post-PCR %s \u2014 %s.", pre_txt, post_txt, verdict))
      }

      if (frac_higher_pre >= 0.75 && batch_pre_median > 0.5) {
        verdict <- sprintf(
          "%.0f%% of this batch's other samples yielded more pre-PCR DNA than this one did (batch median %.3g ng/\u00b5L) \u2014 possible poor extraction specific to this sample, rather than a batch-wide low-biomass pattern; consider repeating extraction",
          100 * frac_higher_pre, batch_pre_median
        )
        return(sprintf("PCR: Pre-PCR %s \u2192 Post-PCR %s \u2014 %s.", pre_txt, post_txt, verdict))
      }
    }
  }

  # Fallback: single-sample thresholds, used when there isn't enough batch
  # context. The low/good post-PCR bands meet at the same boundary
  # (2 ng/uL) so no value can fall through uncategorised.
  pre_high <- pre_val >= 2
  post_low <- post_val < 2
  post_good <- post_val >= 2
  verdict <- if (!pre_detectable && post_good) {
    "efficient amplification from very low/undetected starting material \u2014 typical low-biomass pattern, no inhibition concern"
  } else if (!pre_detectable && post_low) {
    "low starting concentration and low PCR yield \u2014 consistent with genuinely low biomass rather than inhibition"
  } else if (pre_high && post_low) {
    "substantial starting DNA did not reach a strong yield \u2014 possible evidence of PCR inhibition; consider repeating PCR on a diluted aliquot"
  } else if (pre_detectable && post_low) {
    "DNA was present but PCR yield remained modest \u2014 possible evidence of inhibition; consider a diluted repeat"
  } else {
    "good starting material and efficient amplification \u2014 no inhibition concern"
  }
  sprintf("PCR: Pre-PCR %s \u2192 Post-PCR %s \u2014 %s.", pre_txt, post_txt, verdict)
}

add_markdown_paragraph <- function(doc, text, style = "Normal") {
  matches <- str_locate_all(text, "\\*\\*[^*]+\\*\\*")[[1]]
  fps <- list()
  pos <- 1
  if (nrow(matches) > 0) {
    for (i in seq_len(nrow(matches))) {
      if (matches[i, "start"] > pos) {
        fps[[length(fps) + 1]] <- ftext(str_sub(text, pos, matches[i, "start"] - 1), fp_text())
      }
      bold_text <- str_sub(text, matches[i, "start"] + 2, matches[i, "end"] - 2)
      fps[[length(fps) + 1]] <- ftext(bold_text, fp_text(bold = TRUE))
      pos <- matches[i, "end"] + 1
    }
  }
  if (pos <= nchar(text)) {
    fps[[length(fps) + 1]] <- ftext(str_sub(text, pos, nchar(text)), fp_text())
  }
  doc <- body_add_fpar(doc, fpar(fp_p = fp_par(), values = fps), style = style)
  doc
}

# -----------------------------------------------------------------------------
FLAGGED_VERDICTS <- c("Below control", "Possible index hopping", "Possible index hopping (weak evidence)")
TRUSTED_VERDICTS <- c("Trusted (recurring)", "Trusted (recurring, low biomass)")

methodology_paragraphs <- c(
  "This app reviews the **entire** wf-metagenomics counts table for a batch \u2014 every organism, not a curated reporting panel. \"Notable\" calls (used throughout search, the verdicts table, and the export) are those meeting BOTH a configurable minimum read count AND minimum % of that sample's classified reads, adjustable in the sidebar \u2014 requiring both avoids a large absolute count in a deep sample or a large percentage in a shallow one alone qualifying trivial noise, and keeps a table of ~2,000 species x 24 samples tractable.",
  "**Negative control check.** For every notable call, the app looks up that species' read count in the negative control. If a sample's own count for that species is at or below the control's count, the call is flagged **Below control** \u2014 the signal cannot be distinguished from background at that level, regardless of how deep the sample is overall.",
  "**Recurrence check.** For species that pass the control check, the app counts how many other samples in the batch also carry a notable amount of the same species. A species seen only once in the whole batch (a **Singleton**) isn't flagged as contamination, but it also doesn't get the confidence boost that comes from independent corroboration.",
  "**Index-hopping check, with corroboration.** The app also compares each call against the same species' count in every other sample. If a sample's count is a tiny fraction (under 2%) of another sample's count for the same species, and that other sample has far more (over 10x) total reads, that's an initial signature of possible cross-talk. Because a single ratio match is weak evidence on its own, the app then checks whether at least 2 of that carrier's other top 5 species also show up in the flagged sample \u2014 real bulk cross-talk should smear across several species, not just one. With that corroboration, the call is flagged **Possible index hopping**. Without it, if the species independently recurs elsewhere in the batch, it's instead labelled **Trusted (recurring, low biomass)** \u2014 a genuine low-abundance detection, not contamination. If corroboration is weak AND there's no recurrence, it becomes **Possible index hopping (weak evidence)** \u2014 inconclusive, flagged for a low-priority manual check.",
  "**PCR efficiency commentary** (only shown if a position-mapping file with Pre-/Post-PCR values is provided) is based first on two batch-relative comparisons, distinguishing three different failure modes rather than lumping them together: **PCR inhibition** \u2014 a sample's pre-PCR concentration ranks among the highest in the batch, but its post-PCR yield still ranks at or below the batch median (DNA was extracted, but didn't amplify as well as similar- or lower-input samples did); **possible poor extraction** \u2014 most of the batch achieved a decent pre-PCR reading, but this one sample's is an outlier on the low end (the shortfall is upstream of PCR, in extraction, not PCR chemistry); and ordinary **low biomass** \u2014 low pre-PCR is simply the norm for most of the batch, so a low reading here isn't an outlier needing explanation. This adapts to whatever a \"good\" range looks like in this particular batch/kit lot, rather than relying on one fixed number for every batch. When there isn't enough batch context (fewer than 4 samples with both values), it falls back to fixed thresholds: efficient low-biomass amplification, genuinely low biomass throughout, possible inhibition, or good material with efficient amplification.",
  "**GTX case grouping** (only available if the position-mapping spreadsheet includes an auxiliary row \u2014 blank Position, case code in the Swab ID column \u2014 immediately below each sample) rolls samples up into cases on the Case Summary tab and in the exported report, using the same confidence-tier logic, zero-call handling, and within-case species comparison as described above for recurrence and corroboration. A case with zero notable calls checks whether any of its samples shows PCR inhibition or poor-extraction evidence before describing the negative as reportable, exactly as the sample-level PCR commentary would.",
  "**What this does not do:** it cannot tell true low-abundance biology apart from extraction bias, cannot confirm species-level identity against closely related reference genomes, and is not a substitute for targeted PCR or repeat sequencing where a result needs to be reported with certainty. The index-hopping and corroboration checks are heuristics, not a formal statistical test \u2014 treat every verdict as a prompt to check the plate layout and processing order if certainty is required, not a definitive answer on its own.",
  "The **barcode-equals-position** assumption (used only when a position-mapping file is supplied) is not verified against the data \u2014 if a batch was pooled out of position order, every verdict here would still be computed correctly, but attached to the wrong swab ID. Without a position-mapping file, samples are identified directly by their CSV column (swab ID) and no barcode-level read/QC context or PCR commentary is available."
)

genus_notes <- c(
  "Aspergillus" = "conidia can be melanized/pigmented in some species, with a robust chitin-rich cell wall; generally amplifies well with standard fungal extraction, but incomplete mechanical lysis can under-represent it.",
  "Fusarium" = "thin-walled, non-melanized hyphae that lyse and amplify readily; some species produce trichothecene mycotoxins, but these are not typical PCR inhibitors.",
  "Malassezia" = "a lipid-dependent yeast with a thick, lipid-rich cell wall \u2014 standard extraction can under-represent it unless a lipid-disrupting lysis step is used.",
  "Candida" = "yeast with a chitin/glucan cell wall that is generally straightforward to lyse and amplify.",
  "Alternaria" = "dematiaceous (melanized) hyphae \u2014 melanin is a well-documented PCR inhibitor that can co-purify with DNA and suppress amplification.",
  "Penicillium" = "similar profile to Aspergillus \u2014 pigmented conidia in some species, generally amplifies well with standard fungal extraction.",
  "Cladosporium" = "dematiaceous (melanized) hyphae, a known source of co-purifying PCR-inhibitory pigment.",
  "Cryptococcus" = "yeast with a polysaccharide capsule that can resist standard lysis, potentially under-representing it without a capsule-disrupting step.",
  "Mycobacterium" = "waxy, lipid-rich (mycolic acid) cell envelope that resists standard lysis \u2014 often under-represented without a bead-beating or enzymatic step.",
  "Pseudomonas" = "gram-negative, thin cell wall; generally lyses and amplifies readily; some strains are common reagent/kit background organisms."
)

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

ui <- page_navbar(
  title = "Whole-Batch Metagenomics Confidence Explorer",
  theme = bs_theme(version = 5, primary = "#1F4E5F"),

  sidebar = sidebar(
    width = 330,
    h5("1. Upload batch files"),
    fileInput("html_file", "MinKNOW run report (.html)", accept = ".html"),
    fileInput("csv_file", "wf-metagenomics counts (.csv)", accept = ".csv"),
    fileInput("position_file", "Position-mapping spreadsheet (.xlsx) \u2014 optional", accept = c(".xlsx", ".xls")),
    fileInput("logo_file", "Company logo for report (optional, .png/.jpg)", accept = c(".png", ".jpg", ".jpeg")),
    hr(),
    numericInput("low_yield_threshold", "Low-yield reference (classified reads)", value = 5000, min = 100, step = 500),
    numericInput("min_reads", "Notable call: min reads", value = 100, min = 1, step = 10),
    numericInput("min_pct", "Notable call: min % of classified reads", value = 0.2, min = 0, step = 0.05),
    helpText(
      "Without a position-mapping file, samples are identified by their CSV column (swab ID) and depth is measured by classified reads. ",
      "With one, barcode-level total reads/QC and PCR commentary are added."
    ),
    hr(),
    downloadButton("download_report", "Download Full Report (.docx)", class = "btn-primary w-100")
  ),

  nav_panel(
    "Batch Overview",
    card(card_header("Classified reads per sample"),
         plotOutput("overview_plot", height = "420px", fill = FALSE)),
    card(card_header("Per-sample summary"), DTOutput("overview_table"))
  ),

  nav_panel(
    "Taxonomic Composition",
    card(
      card_header("Composition by taxonomic rank"),
      selectInput("tax_rank", "Rank", choices = c("superkingdom","kingdom","phylum","genus"), selected = "superkingdom"),
      plotOutput("composition_plot", height = "450px", fill = FALSE)
    ),
    card(card_header("Composition table"), DTOutput("composition_table"))
  ),

  nav_panel(
    "Negative Control",
    layout_columns(
      col_widths = c(5, 4, 3),
      card(card_header("Control composition (top taxa)"),
           plotOutput("control_plot", height = "380px", fill = FALSE)),
      card(card_header("Negative control read summary"),
           tableOutput("control_summary_table"),
           plotOutput("control_summary_plot", height = "220px", fill = FALSE)),
      card(card_header("About this tab"),
           markdown(
             "Every species listed here came out of the negative control. Any of these that is **also a notable call** elsewhere in the batch is a genuine contamination risk for that species specifically — check the flagged table below."
           ))
    ),
    card(card_header("Notable calls also found in the negative control"), DTOutput("control_flag_table"))
  ),

  nav_panel(
    "Organism Search & Confidence",
    layout_columns(
      col_widths = c(4, 8),
      card(
        card_header("Search"),
        textInput("organism_query", "Organism name (partial match ok)", placeholder = "e.g. aspergillus"),
        uiOutput("species_picker"),
        uiOutput("sample_picker_search")
      ),
      card(card_header("Prevalence across the batch"), plotOutput("prevalence_plot", height = "320px", fill = FALSE))
    ),
    card(card_header("Matches"), DTOutput("search_table")),
    card(
      card_header("Confidence assessment for selected organism + sample"),
      uiOutput("confidence_summary"),
      DTOutput("confidence_breakdown")
    )
  ),

  nav_panel(
    "All Notable Calls & Verdicts",
    card(
      card_header("Every notable (species, sample) call, with automatic verdict"),
      helpText(
        "\"Notable\" = meets BOTH the min-reads AND min-% thresholds set in the sidebar. Verdict logic (checked in order): \"Below control\" = at or below the ",
        "negative control's count for the same species. \"Possible index hopping\" = tiny fraction of another sample's count, that sample has far more ",
        "total reads, AND \u22652 of its other top species also appear here. Weak corroboration + recurrence elsewhere \u2192 \"Trusted (recurring, low biomass)\". ",
        "Weak corroboration + no recurrence \u2192 \"Possible index hopping (weak evidence)\". Otherwise \"Trusted\" (2+ samples) or \"Singleton\"."
      ),
      DTOutput("verdict_table"),
      downloadButton("download_verdicts", "Download full table (CSV)")
    ),
    card(
      card_header("Species vs negative control (select a species)"),
      uiOutput("species_picker_verdict"),
      plotOutput("species_vs_control_plot", height = "350px", fill = FALSE)
    )
  ),

  nav_panel(
    "Case Summary",
    card(
      card_header("GTX case overview"),
      helpText(
        "Only shown when the position-mapping spreadsheet includes GTX case codes (an auxiliary row, with blank Position, ",
        "carrying the case code in the Swab ID column). Confidence tier is computed from the mix of verdicts across all of a ",
        "case's samples: \"High\" = no flagged calls and most samples well-powered; \"Moderate\" = under a third of calls flagged ",
        "and at least half the samples well-powered; \"Low\" = otherwise, including cases where most samples are low-yield."
      ),
      DTOutput("case_overview_table")
    ),
    card(
      card_header("Case detail: sample-by-sample and overall"),
      uiOutput("case_picker"),
      uiOutput("case_detail_view")
    )
  ),

  nav_panel(
    "Sample Detail",
    card(
      card_header("Inspect a single sample"),
      uiOutput("sample_picker_detail"),
      uiOutput("sample_detail_view")
    )
  ),

  nav_panel(
    "Methodology",
    card(card_header("How this app assesses confidence"),
         markdown(paste(methodology_paragraphs, collapse = "\n\n")))
  )
)

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

server <- function(input, output, session) {

  barcode_stats <- reactive({ req(input$html_file); parse_html_barcode_stats(input$html_file$datapath) })
  counts <- reactive({ req(input$csv_file); parse_counts_csv(input$csv_file$datapath) })
  position_map <- reactive({
    if (is.null(input$position_file)) return(NULL)
    tryCatch(parse_position_map(input$position_file$datapath), error = function(e) NULL)
  })

  min_reads_v <- reactive({ v <- input$min_reads; if (is.null(v) || !is.numeric(v)) 100 else v })
  min_pct_v <- reactive({ v <- input$min_pct; if (is.null(v) || !is.numeric(v)) 0.2 else v })
  threshold_v <- reactive({ v <- input$low_yield_threshold; if (is.null(v) || !is.numeric(v)) 5000 else v })

  classified_lookup <- reactive({
    cnts <- counts()
    cnts$long %>%
      group_by(sample_col) %>%
      summarise(
        classified_reads = sum(reads, na.rm = TRUE),
        unknown_reads = sum(reads[species == "Unknown"], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(unknown_pct = ifelse(classified_reads > 0, 100 * unknown_reads / classified_reads, NA))
  })

  # Per-sample overview: always has sample_col/classified_reads/unknown info.
  # If a position-mapping file is supplied, also carries barcode, total_reads,
  # passed_reads_percent, pre_pcr, post_pcr, and a "label" for barcode; else
  # falls back to using the swab ID itself as the label.
  sample_overview <- reactive({
    cnts <- counts()
    cl <- classified_lookup()
    pm <- position_map()
    base <- tibble(sample_col = cnts$sample_cols) %>% left_join(cl, by = "sample_col")
    if (is.null(pm)) {
      base %>% mutate(barcode = NA_character_, swab_id_raw = sample_col, total_reads = NA_real_,
                       passed_reads_percent = NA_real_, pre_pcr = NA_real_, post_pcr = NA_real_,
                       case = NA_character_, label = sample_col) %>%
        arrange(sample_col)
    } else {
      bs <- barcode_stats()
      pm2 <- pm %>% rowwise() %>% mutate(sample_col = match_csv_column(swab_id_raw, cnts$sample_cols)) %>% ungroup()
      pm2 %>% left_join(bs, by = "barcode") %>% left_join(cl, by = "sample_col") %>%
        mutate(label = barcode) %>% arrange(barcode)
    }
  })

  # -------- All notable calls, with verdicts (species x sample, batch-wide) ----

  all_calls_long <- reactive({
    cnts <- counts()
    ctrl_col <- cnts$control_col
    cl <- classified_lookup()
    ov <- sample_overview()
    min_r <- min_reads_v(); min_p <- min_pct_v()

    base <- cnts$long %>%
      filter(reads > 0, species != "Unknown") %>%
      left_join(cl, by = "sample_col") %>%
      mutate(pct = ifelse(classified_reads > 0, 100 * reads / classified_reads, NA)) %>%
      filter(reads >= min_r & (!is.na(pct) & pct >= min_p))

    validate(need(nrow(base) > 0, "No calls meet the current notability thresholds \u2014 try lowering them in the sidebar."))

    base <- base %>%
      left_join(ov %>% select(sample_col, label, total_reads, case, swab_id_raw), by = "sample_col") %>%
      mutate(sample_total_reads = ifelse(is.na(total_reads), classified_reads, total_reads)) %>%
      select(-total_reads)

    ctrl_lookup <- if (!is.na(ctrl_col)) {
      cnts$long %>% filter(sample_col == ctrl_col) %>% select(species, control_reads = reads)
    } else {
      tibble(species = character(0), control_reads = numeric(0))
    }
    base <- base %>% left_join(ctrl_lookup, by = "species") %>%
      mutate(control_reads = ifelse(is.na(control_reads), 0, control_reads))

    calls_df <- base %>%
      group_by(species) %>%
      mutate(n_samples = n()) %>%
      group_modify(~ {
        df <- .x
        df$top_carrier_sample <- NA_character_
        df$top_carrier_reads <- NA_real_
        df$top_carrier_total_reads <- NA_real_
        for (i in seq_len(nrow(df))) {
          other <- df[-i, , drop = FALSE]
          if (nrow(other) == 0) next
          top <- other[which.max(other$reads), ]
          df$top_carrier_sample[i] <- top$sample_col[1]
          df$top_carrier_reads[i] <- top$reads[1]
          df$top_carrier_total_reads[i] <- top$sample_total_reads[1]
        }
        df
      }) %>%
      ungroup() %>%
      mutate(
        hop_ratio_hits = ifelse(!is.na(top_carrier_reads) & top_carrier_reads > 0, reads / top_carrier_reads, NA_real_),
        hop_ratio_depth = ifelse(!is.na(top_carrier_total_reads) & sample_total_reads > 0,
                                  top_carrier_total_reads / sample_total_reads, NA_real_),
        possible_hop = !is.na(hop_ratio_hits) & !is.na(hop_ratio_depth) & hop_ratio_hits < 0.02 & hop_ratio_depth > 10
      )

    calls_df$hop_n_checked <- NA_integer_
    calls_df$hop_n_matched <- NA_integer_
    hop_idx <- which(calls_df$possible_hop)
    for (i in hop_idx) {
      companions <- cnts$long %>%
        filter(sample_col == calls_df$top_carrier_sample[i],
               species != calls_df$species[i], species != "Unknown", reads > 0) %>%
        arrange(desc(reads)) %>% head(5)
      calls_df$hop_n_checked[i] <- nrow(companions)
      if (nrow(companions) > 0) {
        matched <- cnts$long %>% filter(sample_col == calls_df$sample_col[i], species %in% companions$species, reads > 0)
        calls_df$hop_n_matched[i] <- nrow(matched)
      } else {
        calls_df$hop_n_matched[i] <- 0
      }
    }

    calls_df %>%
      mutate(
        hop_weak_corroboration = possible_hop & !is.na(hop_n_matched) & hop_n_matched <= 1,
        verdict = case_when(
          control_reads > 0 & reads <= control_reads ~ "Below control",
          possible_hop & !hop_weak_corroboration ~ "Possible index hopping",
          possible_hop & hop_weak_corroboration & n_samples >= 2 ~ "Trusted (recurring, low biomass)",
          possible_hop & hop_weak_corroboration ~ "Possible index hopping (weak evidence)",
          n_samples >= 2 ~ "Trusted (recurring)",
          TRUE ~ "Singleton"
        ),
        pct_of_classified = round(pct, 3)
      ) %>%
      arrange(species, desc(reads))
  })

  # ---- Case Summary (roll-up of sample verdicts by GTX case) --------------

  sample_call_summary <- reactive({
    calls <- all_calls_long()
    if (nrow(calls) == 0) {
      return(tibble(sample_col = character(), n_calls = integer(), n_trusted = integer(),
                    n_singleton = integer(), n_below_control = integer(), n_hop = integer(),
                    species_trusted = character(), species_flagged = character()))
    }
    calls %>%
      group_by(sample_col) %>%
      summarise(
        n_calls = n(),
        n_trusted = sum(verdict %in% TRUSTED_VERDICTS),
        n_singleton = sum(verdict == "Singleton"),
        n_below_control = sum(verdict == "Below control"),
        n_hop = sum(verdict %in% c("Possible index hopping", "Possible index hopping (weak evidence)")),
        species_trusted = paste(species[verdict %in% c(TRUSTED_VERDICTS, "Singleton")], collapse = ", "),
        species_flagged = paste(species[verdict %in% FLAGGED_VERDICTS], collapse = ", "),
        .groups = "drop"
      )
  })

  sample_level <- reactive({
    thr <- threshold_v()
    ov <- sample_overview()
    scs <- sample_call_summary()
    ov %>%
      left_join(scs, by = "sample_col") %>%
      mutate(
        n_calls = ifelse(is.na(n_calls), 0L, n_calls),
        n_trusted = ifelse(is.na(n_trusted), 0L, n_trusted),
        n_singleton = ifelse(is.na(n_singleton), 0L, n_singleton),
        n_below_control = ifelse(is.na(n_below_control), 0L, n_below_control),
        n_hop = ifelse(is.na(n_hop), 0L, n_hop),
        species_trusted = ifelse(is.na(species_trusted), "", species_trusted),
        species_flagged = ifelse(is.na(species_flagged), "", species_flagged),
        well_powered = !is.na(classified_reads) & classified_reads >= thr
      )
  })

  case_summary <- reactive({
    sl <- sample_level() %>% filter(!is.na(case))
    validate(need(nrow(sl) > 0, "No GTX case information found \u2014 the position-mapping spreadsheet needs an auxiliary row (blank Position) carrying the case code."))
    sl %>%
      group_by(case) %>%
      summarise(
        n_samples = n(),
        n_well_powered = sum(well_powered, na.rm = TRUE),
        total_calls = sum(n_calls),
        n_trusted = sum(n_trusted),
        n_singleton = sum(n_singleton),
        n_below_control = sum(n_below_control),
        n_hop = sum(n_hop),
        .groups = "drop"
      ) %>%
      mutate(
        n_low_yield = n_samples - n_well_powered,
        flagged = n_below_control + n_hop,
        flagged_fraction = ifelse(total_calls > 0, flagged / total_calls, 0),
        power_fraction = n_well_powered / n_samples,
        tier = case_when(
          total_calls == 0 ~ "No positive calls",
          flagged_fraction == 0 & power_fraction >= 0.7 ~ "High",
          flagged_fraction < 0.34 & power_fraction >= 0.5 ~ "Moderate",
          TRUE ~ "Low"
        )
      ) %>%
      arrange(case)
  })

  # Plain-English narrative for one sample row.
  sample_narrative <- function(row, batch_pre = NULL, batch_post = NULL) {
    power_txt <- if (isTRUE(row$well_powered)) "well-powered" else "low-depth"
    depth_n <- if (is.na(row$classified_reads)) 0 else row$classified_reads
    pcr_txt <- if (!is.na(row$pre_pcr) || !is.na(row$post_pcr)) pcr_comment(row$pre_pcr, row$post_pcr, batch_pre, batch_post) else NULL
    base_txt <- if (row$n_calls == 0) {
      sprintf("%s (%s): %s (%s classified reads), no notable calls at the current thresholds.",
              row$label, row$swab_id_raw, power_txt, format(depth_n, big.mark = ","))
    } else {
      parts <- c()
      if (row$n_trusted > 0) parts <- c(parts, sprintf("%d trusted", row$n_trusted))
      if (row$n_singleton > 0) parts <- c(parts, sprintf("%d singleton", row$n_singleton))
      if (row$n_below_control > 0) parts <- c(parts, sprintf("%d below-control", row$n_below_control))
      if (row$n_hop > 0) parts <- c(parts, sprintf("%d possible-hopping", row$n_hop))
      call_txt <- paste(parts, collapse = ", ")
      n_flagged <- row$n_below_control + row$n_hop
      strength <- if (n_flagged == 0 & row$n_trusted > 0) "a solid, reliable contributor to this case's profile"
                  else if (n_flagged > 0 & row$n_trusted == 0 & row$n_singleton == 0) "an unreliable contributor \u2014 none of its calls hold up"
                  else if (n_flagged > 0) "a mixed contributor \u2014 some calls trustworthy, others do not hold up"
                  else "an unconfirmed contributor \u2014 its call(s) lack independent corroboration"
      sprintf("%s (%s): %s (%s classified reads), %d notable call(s) [%s] \u2014 %s.",
              row$label, row$swab_id_raw, power_txt, format(depth_n, big.mark = ","), row$n_calls, call_txt, strength)
    }
    if (is.null(pcr_txt)) base_txt else paste(base_txt, pcr_txt)
  }

  case_conclusion_text <- function(case_row, sl_case, batch_pre = NULL, batch_post = NULL) {
    trusted_species <- sl_case$species_trusted[nzchar(sl_case$species_trusted)]
    trusted_species <- unique(unlist(str_split(paste(trusted_species, collapse = ", "), ",\\s*")))
    trusted_species <- trusted_species[nzchar(trusted_species)]
    flagged_species <- sl_case$species_flagged[nzchar(sl_case$species_flagged)]
    flagged_species <- unique(unlist(str_split(paste(flagged_species, collapse = ", "), ",\\s*")))
    flagged_species <- flagged_species[nzchar(flagged_species)]

    reason_parts <- c()
    if (case_row$flagged > 0) {
      reason_parts <- c(reason_parts, sprintf(
        "%d call(s) were flagged and should be excluded from reporting (see the removal list below)", case_row$flagged))
    }
    if (case_row$power_fraction < 0.7) {
      reason_parts <- c(reason_parts, sprintf(
        "only %d of %d sample(s) were well-powered", case_row$n_well_powered, case_row$n_samples))
    }

    if (case_row$total_calls == 0) {
      inhibited <- FALSE
      if (!is.null(batch_pre)) {
        inhibited <- any(mapply(function(pr, po) {
          if (is.na(pr) && is.na(po)) return(FALSE)
          txt <- pcr_comment(pr, po, batch_pre, batch_post)
          str_detect(txt, "inhibition|poor extraction")
        }, sl_case$pre_pcr, sl_case$post_pcr))
      }
      tier_reason <- if (inhibited) {
        "no notable calls were made, and at least one sample shows batch-relative evidence of PCR inhibition or poor extraction (see the PCR commentary above) \u2014 this should NOT be treated as a confirmed absence of organisms; repeat extraction/PCR before concluding a clean negative."
      } else {
        "no notable calls were made in this case's sample(s), with no PCR inhibition or extraction concerns flagged \u2014 can be treated as genuinely low-signal, subject to normal review."
      }
    } else if (length(reason_parts) == 0) {
      tier_reason <- "no calls were flagged and sample depth was adequate throughout \u2014 this case's profile can be reported with good confidence."
    } else {
      tier_reason <- paste0(paste(reason_parts, collapse = "; "), ".")
    }

    lines <- c(
      sprintf("Case %s \u2014 %d sample(s), %d well-powered / %d low-yield. %d total notable call(s): %d trusted, %d singleton, %d below-control, %d possible index-hopping.",
              case_row$case, case_row$n_samples, case_row$n_well_powered, case_row$n_low_yield,
              case_row$total_calls, case_row$n_trusted, case_row$n_singleton, case_row$n_below_control, case_row$n_hop),
      sprintf("Overall confidence: %s \u2014 %s", case_row$tier, tier_reason)
    )
    if (length(trusted_species) > 0) {
      lines <- c(lines, sprintf("Species supported with reasonable confidence: %s.", paste(trusted_species, collapse = ", ")))
    }
    if (length(flagged_species) > 0) {
      lines <- c(lines, sprintf("Species NOT to be reported from this case (flagged): %s.", paste(flagged_species, collapse = ", ")))
    }

    if (nrow(sl_case) >= 2 && length(trusted_species) > 0) {
      sample_species_lists <- lapply(sl_case$species_trusted, function(s) {
        if (!nzchar(s)) return(character(0))
        trimws(str_split(s, ",")[[1]])
      })
      species_sample_count <- sapply(trusted_species, function(sp) {
        sum(sapply(sample_species_lists, function(lst) sp %in% lst))
      })
      consistent <- names(species_sample_count)[species_sample_count >= 2]
      case_specific <- names(species_sample_count)[species_sample_count == 1]
      if (length(consistent) > 0) {
        lines <- c(lines, sprintf(
          "Independently detected in 2+ of this case's own %d sample(s) \u2014 strong internal consistency: %s.",
          nrow(sl_case), paste(consistent, collapse = ", ")))
      }
      if (length(case_specific) > 0) {
        lines <- c(lines, sprintf(
          "Detected in only 1 of this case's own samples (no internal corroboration within this case, though still trusted via batch-wide recurrence): %s.",
          paste(case_specific, collapse = ", ")))
      }
    }
    lines
  }

  case_removal_items <- function(calls, case_name) {
    if (nrow(calls) == 0) return(character(0))
    flagged <- calls %>% filter(case == case_name, verdict %in% FLAGGED_VERDICTS)
    if (nrow(flagged) == 0) return(character(0))
    sapply(seq_len(nrow(flagged)), function(i) {
      r <- flagged[i, ]
      if (r$verdict == "Below control") {
        sprintf("Remove %s from %s (%s) \u2014 %d reads, at or below the negative control's %d reads for this species.",
                r$species, r$label, r$swab_id_raw, r$reads, r$control_reads)
      } else if (r$verdict == "Possible index hopping") {
        sprintf("Remove/verify %s from %s (%s) \u2014 possible index hopping from %s (%d reads), corroborated by %d/%d of that sample's other top species also appearing here.",
                r$species, r$label, r$swab_id_raw, r$top_carrier_sample, r$top_carrier_reads, r$hop_n_matched, r$hop_n_checked)
      } else {
        sprintf("Review (low priority) %s from %s (%s) \u2014 weak/uncorroborated hop signal, not independently recurring; verify against plate layout if certainty is required.",
                r$species, r$label, r$swab_id_raw)
      }
    })
  }

  case_organism_comment <- function(species_list) {
    if (length(species_list) == 0) {
      return("No species were confidently detected in this case, so no organism-specific extraction/inhibition comment applies.")
    }
    genera <- unique(word(species_list, 1))
    notes <- genus_notes[genera]
    has_note <- !is.na(notes)
    if (!any(has_note)) {
      return(sprintf("Genera detected (%s) carry no specific extraction or PCR-inhibition concerns in this app's reference table.",
                      paste(genera, collapse = ", ")))
    }
    flagged_notes <- sprintf("%s \u2014 %s", genera[has_note], str_remove(notes[has_note], "\\.$"))
    other_genera <- genera[!has_note]
    txt <- paste0("Organisms detected include: ", paste(flagged_notes, collapse = "; "), ".")
    if (length(other_genera) > 0) {
      txt <- paste0(txt, " Other genera detected (", paste(other_genera, collapse = ", "),
                    ") carry no specific extraction/inhibition notes in this app's reference table.")
    }
    txt
  }

  output$case_overview_table <- renderDT({
    cs <- case_summary()
    out <- cs %>%
      transmute(
        Case = case, Samples = n_samples, `Well-powered` = n_well_powered,
        `Low-yield` = n_low_yield, `Total Calls` = total_calls, Trusted = n_trusted,
        Singleton = n_singleton, `Below Control` = n_below_control, `Possible Hopping` = n_hop,
        `Confidence Tier` = tier
      )
    datatable(out, options = list(pageLength = 20), rownames = FALSE)
  })

  output$case_picker <- renderUI({
    cs <- case_summary()
    validate(need(nrow(cs) > 0, "Upload files with a position-mapping spreadsheet that includes case codes."))
    selectInput("selected_case", "GTX Case", choices = sort(cs$case))
  })

  output$case_detail_view <- renderUI({
    req(input$selected_case)
    cs <- case_summary() %>% filter(case == input$selected_case)
    validate(need(nrow(cs) == 1, "Case not found."))
    sl_full <- sample_level()
    sl <- sl_full %>% filter(case == input$selected_case) %>% arrange(barcode)
    calls <- all_calls_long()

    sample_items <- lapply(seq_len(nrow(sl)), function(i) tags$li(sample_narrative(sl[i, ], sl_full$pre_pcr, sl_full$post_pcr)))
    conclusion_lines <- case_conclusion_text(cs, sl, sl_full$pre_pcr, sl_full$post_pcr)
    conclusion_items <- lapply(conclusion_lines, function(l) tags$p(l))

    removal_items <- case_removal_items(calls, input$selected_case)
    removal_block <- if (length(removal_items) == 0) {
      tags$p("No calls need to be removed from this case's report.")
    } else {
      tags$ul(lapply(removal_items, tags$li))
    }

    trusted_species <- sl$species_trusted[nzchar(sl$species_trusted)]
    trusted_species <- unique(unlist(str_split(paste(trusted_species, collapse = ", "), ",\\s*")))
    trusted_species <- trusted_species[nzchar(trusted_species)]
    organism_text <- case_organism_comment(trusted_species)

    tagList(
      tags$h5("Sample-by-sample"),
      tags$ul(sample_items),
      tags$h5("Case as a whole"),
      conclusion_items,
      tags$h5("To be removed from this case's report"),
      removal_block,
      tags$h5("Organisms present \u2014 extraction / PCR-inhibition notes"),
      tags$p(organism_text)
    )
  })

  # ---- Batch Overview -------------------------------------------------------

  output$overview_plot <- renderPlot({
    df <- sample_overview()
    validate(need(nrow(df) > 0, "Upload files to see the overview."))
    thr <- threshold_v()
    ggplot(df, aes(x = reorder(label, classified_reads), y = pmax(classified_reads, 1))) +
      geom_col(fill = "#1F4E5F") +
      geom_text(aes(label = label), vjust = -0.4, size = 3, color = "black") +
      scale_y_log10(labels = label_comma()) +
      geom_hline(yintercept = thr, linetype = "dashed", color = "red") +
      labs(x = NULL, y = "Classified reads (log scale)", title = "Classified reads per sample") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$overview_table <- renderDT({
    df <- sample_overview() %>%
      transmute(
        Sample = label, `Swab ID` = swab_id_raw,
        `Total Reads` = total_reads, `Passed Reads %` = passed_reads_percent,
        `Classified Reads` = classified_reads, `Unknown %` = round(unknown_pct, 1),
        `Pre-PCR` = pre_pcr, `Post-PCR` = post_pcr
      )
    datatable(df, options = list(pageLength = 24), rownames = FALSE)
  })

  # ---- Taxonomic Composition -------------------------------------------------

  composition_data <- reactive({
    cnts <- counts()
    rank <- input$tax_rank
    validate(need(rank %in% cnts$taxonomy_present, "This rank is not present in the counts CSV."))
    ov <- sample_overview()
    cnts$long %>%
      filter(sample_col != cnts$control_col | is.na(cnts$control_col)) %>%
      rename(rank_val = !!rank) %>%
      group_by(sample_col, rank_val) %>%
      summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop") %>%
      group_by(sample_col) %>%
      mutate(pct = 100 * reads / sum(reads)) %>%
      ungroup() %>%
      left_join(ov %>% select(sample_col, label), by = "sample_col")
  })

  output$composition_plot <- renderPlot({
    df <- composition_data()
    validate(need(nrow(df) > 0, "No composition data."))
    # collapse rare categories per-plot to keep the legend readable
    top_cats <- df %>% group_by(rank_val) %>% summarise(tot = sum(reads)) %>% arrange(desc(tot)) %>% head(12) %>% pull(rank_val)
    df <- df %>% mutate(rank_val2 = ifelse(rank_val %in% top_cats, rank_val, "Other"))
    ggplot(df, aes(x = label, y = pct, fill = rank_val2)) +
      geom_col() +
      labs(x = NULL, y = "% of classified reads", fill = str_to_title(input$tax_rank),
           title = paste0("Composition by ", input$tax_rank)) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$composition_table <- renderDT({
    df <- composition_data() %>%
      transmute(Sample = label, !!str_to_title(input$tax_rank) := rank_val, Reads = reads, `% of Sample` = round(pct, 2)) %>%
      arrange(Sample, desc(Reads))
    datatable(df, options = list(pageLength = 25), rownames = FALSE)
  })

  # ---- Negative Control -------------------------------------------------------

  output$control_plot <- renderPlot({
    cnts <- counts()
    validate(need(!is.na(cnts$control_col), "No control column detected in the counts CSV."))
    ctrl <- cnts$long %>% filter(sample_col == cnts$control_col, reads > 0)
    total <- sum(ctrl$reads)
    validate(need(total > 0, "Negative control has zero classified reads."))
    top <- ctrl %>% mutate(pct = 100 * reads / total) %>% arrange(desc(pct)) %>% head(10)
    ggplot(top, aes(x = reorder(species, pct), y = pct)) +
      geom_col(fill = "#1F4E5F") +
      coord_flip() +
      labs(x = NULL, y = "% of classified reads", title = paste0("Control composition (", total, " classified reads)")) +
      theme_minimal(base_size = 12)
  })

  control_summary <- reactive({
    cnts <- counts()
    ctrl_col <- cnts$control_col
    validate(need(!is.na(ctrl_col), "No negative control column found."))
    ov <- sample_overview()
    row <- ov %>% filter(sample_col == ctrl_col)
    validate(need(nrow(row) == 1, "Could not locate the negative control's read/QC data."))
    total_reads <- ifelse(is.na(row$total_reads[1]), row$classified_reads[1], row$total_reads[1])
    classified_reads <- ifelse(is.na(row$classified_reads[1]), 0, row$classified_reads[1])
    unknown_reads <- ifelse(is.na(row$unknown_reads[1]), 0, row$unknown_reads[1])
    mapped_reads <- classified_reads - unknown_reads
    tibble(
      Metric = factor(c("Total Reads", "Classified Reads", "Mapped Reads", "Unknown Reads"),
                       levels = c("Total Reads", "Classified Reads", "Mapped Reads", "Unknown Reads")),
      Value = c(total_reads, classified_reads, mapped_reads, unknown_reads)
    )
  })

  output$control_summary_table <- renderTable({
    control_summary() %>% mutate(Value = format(Value, big.mark = ","))
  }, striped = TRUE, bordered = TRUE, width = "100%", colnames = TRUE)

  output$control_summary_plot <- renderPlot({
    df <- control_summary()
    ggplot(df, aes(x = Metric, y = Value)) +
      geom_col(fill = "#1F4E5F") +
      geom_text(aes(label = format(Value, big.mark = ",")), vjust = -0.4, size = 3.3) +
      labs(x = NULL, y = "Reads") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
      expand_limits(y = max(df$Value) * 1.15)
  })

  output$control_flag_table <- renderDT({
    calls <- all_calls_long()
    flagged_species <- calls %>% filter(control_reads > 0) %>% distinct(species, control_reads)
    validate(need(nrow(flagged_species) > 0, "No notable calls were found in the negative control \u2014 no contamination flags."))
    out <- calls %>%
      filter(species %in% flagged_species$species) %>%
      transmute(Species = species, Sample = label, Reads = reads, `Control Reads` = control_reads,
                Verdict = ifelse(verdict %in% FLAGGED_VERDICTS, paste0("\u26a0 ", verdict), verdict)) %>%
      arrange(Species, desc(Reads))
    datatable(out, options = list(pageLength = 20), rownames = FALSE)
  })

  # ---- Organism Search --------------------------------------------------------

  matching_species <- reactive({
    req(input$organism_query)
    cnts <- counts()
    unique(cnts$long$species[str_detect(cnts$long$species, regex(input$organism_query, ignore_case = TRUE))])
  })

  output$species_picker <- renderUI({
    req(input$organism_query)
    ms <- matching_species()
    validate(need(length(ms) > 0, "No species match this search."))
    selectInput("selected_species", "Select species", choices = sort(ms))
  })

  species_hits <- reactive({
    req(input$selected_species)
    cnts <- counts()
    ov <- sample_overview()
    cl <- classified_lookup()
    cnts$long %>%
      filter(species == input$selected_species, reads > 0) %>%
      left_join(ov %>% select(sample_col, label, total_reads), by = "sample_col") %>%
      left_join(cl, by = "sample_col") %>%
      mutate(pct_of_classified = ifelse(classified_reads > 0, 100 * reads / classified_reads, NA)) %>%
      arrange(desc(reads))
  })

  output$sample_picker_search <- renderUI({
    req(input$selected_species)
    hits <- species_hits()
    validate(need(nrow(hits) > 0, "This species has no positive read counts."))
    choices <- setNames(hits$sample_col, paste0(hits$label, " (", hits$reads, " reads)"))
    selectInput("selected_sample_search", "Select sample to assess", choices = choices)
  })

  output$prevalence_plot <- renderPlot({
    req(input$selected_species)
    hits <- species_hits()
    validate(need(nrow(hits) > 0, "No data to plot."))
    hits <- hits %>% mutate(is_selected = if (is.null(input$selected_sample_search)) FALSE else sample_col == input$selected_sample_search)
    ggplot(hits, aes(x = reorder(label, reads), y = reads, fill = is_selected)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c(`TRUE` = "#c0504d", `FALSE` = "#1F4E5F"), guide = "none") +
      labs(x = NULL, y = "Reads", title = paste0("Detections of \"", input$selected_species, "\" across the batch")) +
      theme_minimal(base_size = 13)
  })

  output$search_table <- renderDT({
    req(input$selected_species)
    hits <- species_hits() %>%
      transmute(Sample = label, `Swab ID` = sample_col, Reads = reads,
                `% of Classified` = round(pct_of_classified, 2),
                `Sample Classified Reads` = classified_reads, `Sample Total Reads` = total_reads)
    datatable(hits, options = list(pageLength = 24), rownames = FALSE)
  })

  confidence_result <- reactive({
    req(input$selected_species, input$selected_sample_search)
    calls <- all_calls_long()
    row <- calls %>% filter(species == input$selected_species, sample_col == input$selected_sample_search)
    if (nrow(row) == 0) {
      hits <- species_hits() %>% filter(sample_col == input$selected_sample_search)
      validate(need(nrow(hits) > 0, "No data for this selection."))
      return(list(below_threshold = TRUE, reads = hits$reads[1], pct = hits$pct_of_classified[1]))
    }
    list(below_threshold = FALSE, row = row[1, ])
  })

  output$confidence_summary <- renderUI({
    res <- confidence_result()
    if (isTRUE(res$below_threshold)) {
      return(div(
        style = "padding: 14px 18px; border-radius: 8px; background:#f0f0f022; border:1px solid #888;",
        h4("Below notability threshold", style = "margin:0;"),
        p(sprintf("%d reads (%.2f%% of classified) \u2014 below the sidebar's min-reads/min-%% thresholds, so no automatic verdict is computed. Lower the thresholds to include it.",
                   res$reads, ifelse(is.na(res$pct), 0, res$pct)), style = "margin-top:6px; margin-bottom:0;")
      ))
    }
    r <- res$row
    colour <- if (r$verdict %in% FLAGGED_VERDICTS) "#c62828" else if (r$verdict %in% TRUSTED_VERDICTS) "#2e7d32" else "#f9a825"
    div(
      style = sprintf("padding: 14px 18px; border-radius: 8px; background:%s22; border:1px solid %s;", colour, colour),
      h4(r$verdict, style = sprintf("color:%s; margin:0;", colour)),
      p(sprintf("%d reads, %.2f%% of classified reads in this sample.", r$reads, ifelse(is.na(r$pct_of_classified), 0, r$pct_of_classified)),
        style = "margin-top:6px; margin-bottom:0;")
    )
  })

  output$confidence_breakdown <- renderDT({
    res <- confidence_result()
    validate(need(!isTRUE(res$below_threshold), "No breakdown \u2014 this call is below the notability threshold."))
    r <- res$row
    breakdown <- tibble(
      Factor = c("Reads", "Control reads (same species)", "Recurrence", "Top carrier (if hop-tested)", "Hop corroboration"),
      Detail = c(
        as.character(r$reads),
        as.character(r$control_reads),
        sprintf("present in %d sample(s) total (incl. this one)", r$n_samples),
        ifelse(!is.na(r$top_carrier_sample) & r$possible_hop, paste0(r$top_carrier_sample, " (", r$top_carrier_reads, " reads)"), "n/a"),
        ifelse(!is.na(r$hop_n_checked), paste0(r$hop_n_matched, "/", r$hop_n_checked, " companion species also present"), "n/a")
      )
    )
    datatable(breakdown, options = list(dom = "t"), rownames = FALSE)
  })

  # ---- All Notable Calls & Verdicts -----------------------------------------

  output$verdict_table <- renderDT({
    calls <- all_calls_long()
    out <- calls %>%
      transmute(Species = species, Sample = label, Reads = reads, `% Classified` = pct_of_classified,
                `Sample Classified Reads` = classified_reads, `Control Reads` = control_reads,
                `# Samples` = n_samples,
                `Top Carrier` = ifelse(!is.na(top_carrier_sample) & possible_hop, paste0(top_carrier_sample, " (", top_carrier_reads, " reads)"), ""),
                `Hop Corroboration` = ifelse(!is.na(hop_n_checked), paste0(hop_n_matched, "/", hop_n_checked), ""),
                Verdict = ifelse(verdict %in% FLAGGED_VERDICTS, paste0("\u26a0 ", verdict), verdict))
    datatable(out, filter = "top", options = list(pageLength = 25), rownames = FALSE)
  })

  output$download_verdicts <- downloadHandler(
    filename = function() "all_notable_calls.csv",
    content = function(file) write.csv(all_calls_long(), file, row.names = FALSE)
  )

  output$species_picker_verdict <- renderUI({
    calls <- all_calls_long()
    validate(need(nrow(calls) > 0, "No species to select."))
    selectInput("selected_species_verdict", "Species", choices = sort(unique(calls$species)))
  })

  output$species_vs_control_plot <- renderPlot({
    req(input$selected_species_verdict)
    calls <- all_calls_long() %>% filter(species == input$selected_species_verdict)
    validate(need(nrow(calls) > 0, "No data for this species."))
    ctrl_reads <- calls$control_reads[1]
    df <- calls %>%
      transmute(lab = label, reads, group = ifelse(reads <= ctrl_reads & ctrl_reads > 0, "below", "above")) %>%
      arrange(reads)
    if (ctrl_reads > 0) df <- bind_rows(tibble(lab = "Control", reads = ctrl_reads, group = "ctrl"), df)
    ggplot(df, aes(x = reorder(lab, reads), y = reads, fill = group)) +
      geom_col() +
      scale_fill_manual(values = c(ctrl = "#333333", below = "#c0504d", above = "#2e7d32"), guide = "none") +
      { if (ctrl_reads > 0) geom_hline(yintercept = ctrl_reads, linetype = "dashed") } +
      labs(x = NULL, y = "Reads", title = paste0(input$selected_species_verdict, ": all notable calls vs control")) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
  })

  # ---- Sample Detail -----------------------------------------------------------

  output$sample_picker_detail <- renderUI({
    ov <- sample_overview()
    validate(need(nrow(ov) > 0, "Upload files first."))
    choices <- setNames(ov$sample_col, paste0(ov$label, " (", ov$swab_id_raw, ")"))
    selectInput("selected_sample_detail", "Sample", choices = choices)
  })

  output$sample_detail_view <- renderUI({
    req(input$selected_sample_detail)
    calls <- all_calls_long() %>% filter(sample_col == input$selected_sample_detail) %>% arrange(desc(reads))
    ov_full <- sample_overview()
    ov <- ov_full %>% filter(sample_col == input$selected_sample_detail)
    pcr_block <- if (nrow(ov) == 1 && !is.na(ov$pre_pcr[1]) || (nrow(ov) == 1 && !is.na(ov$post_pcr[1]))) {
      tags$p(pcr_comment(ov$pre_pcr[1], ov$post_pcr[1], ov_full$pre_pcr, ov_full$post_pcr))
    } else NULL
    if (nrow(calls) == 0) {
      return(tagList(pcr_block, tags$p("No notable calls for this sample at the current thresholds.")))
    }
    tbl <- calls %>%
      transmute(Species = species, Reads = reads, `% Classified` = pct_of_classified,
                `Control Reads` = control_reads, `# Samples` = n_samples,
                Verdict = ifelse(verdict %in% FLAGGED_VERDICTS, paste0("\u26a0 ", verdict), verdict))
    tagList(
      pcr_block,
      renderDT(datatable(tbl, options = list(pageLength = 25), rownames = FALSE))
    )
  })

  # ---- Export ------------------------------------------------------------------

  output$download_report <- downloadHandler(
    filename = function() paste0("WholeBatch_Confidence_Report_", format(Sys.Date(), "%Y%m%d"), ".docx"),
    content = function(file) {
      req(input$html_file, input$csv_file)

      thr <- threshold_v()
      ov <- sample_overview()
      calls <- all_calls_long()
      cnts <- counts()

      chart_path <- tempfile(fileext = ".png")
      p <- ggplot(ov, aes(x = reorder(label, classified_reads), y = pmax(classified_reads, 1))) +
        geom_col(fill = "#1F4E5F") +
        geom_text(aes(label = label), vjust = -0.4, size = 3, color = "black") +
        scale_y_log10(labels = label_comma()) +
        geom_hline(yintercept = thr, linetype = "dashed", color = "red") +
        labs(x = NULL, y = "Classified reads (log scale)", title = "Classified reads per sample") +
        theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(chart_path, p, width = 9, height = 5, dpi = 200)

      comp_path <- tempfile(fileext = ".png")
      comp_df <- composition_data()
      top_cats <- comp_df %>% group_by(rank_val) %>% summarise(tot = sum(reads)) %>% arrange(desc(tot)) %>% head(12) %>% pull(rank_val)
      comp_df <- comp_df %>% mutate(rank_val2 = ifelse(rank_val %in% top_cats, rank_val, "Other"))
      pc <- ggplot(comp_df, aes(x = label, y = pct, fill = rank_val2)) +
        geom_col() +
        labs(x = NULL, y = "% of classified reads", fill = str_to_title(input$tax_rank),
             title = paste0("Composition by ", input$tax_rank)) +
        theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(comp_path, pc, width = 9, height = 5, dpi = 200)

      doc <- read_docx()
      if (!is.null(input$logo_file)) {
        tryCatch({ doc <- doc %>% body_add_img(input$logo_file$datapath, width = 2, height = 2, style = "centered") },
                 error = function(e) NULL)
      }
      doc <- doc %>%
        body_add_par("Whole-Batch Metagenomics Confidence Report", style = "heading 1") %>%
        body_add_par(paste0("Generated: ", format(Sys.time(), "%d %B %Y, %H:%M")), style = "Normal") %>%
        body_add_par("", style = "Normal")

      run_info <- tryCatch(parse_html_run_info(input$html_file$datapath), error = function(e) NULL)
      if (!is.null(run_info)) {
        doc <- doc %>% body_add_par("Run Information", style = "heading 1")
        hdr <- run_info$header
        run_summary_df <- data.frame(
          Field = c("Experiment name", "Sample ID", "Device", "Position", "Protocol run ID",
                    "Flow cell ID", "Run status", "Run start", "Run end", "Estimated N50"),
          Value = c(
            ifelse(is.null(hdr$experiment_name), "\u2014", hdr$experiment_name),
            ifelse(is.null(hdr$sample_id), "\u2014", hdr$sample_id),
            ifelse(is.null(hdr$device_type), "\u2014", paste0(hdr$device_type, " (", ifelse(is.null(hdr$serial), "?", hdr$serial), ")")),
            ifelse(is.null(hdr$position), "\u2014", hdr$position),
            ifelse(is.null(hdr$protocol_run_id), "\u2014", hdr$protocol_run_id),
            ifelse(is.null(run_info$flow_cell_id), "\u2014", run_info$flow_cell_id),
            ifelse(is.null(run_info$run_status), "\u2014",
                   paste0(run_info$run_status, ifelse(is.null(run_info$run_status_additional_context), "",
                                                       paste0(" (", run_info$run_status_additional_context, ")")))),
            format_iso_time(run_info$run_start_time), format_iso_time(run_info$run_end_time),
            ifelse(is.null(run_info$estimated_n50), "\u2014", as.character(run_info$estimated_n50))
          ), stringsAsFactors = FALSE
        )
        doc <- doc %>% body_add_table(run_summary_df, style = "table_template")
        doc <- doc %>% body_add_par("Run Configuration", style = "heading 1")
        config_df <- rbind(titlevalue_to_df(run_info$run_setup), titlevalue_to_df(run_info$run_settings))
        if (nrow(config_df) > 0) doc <- doc %>% body_add_table(config_df, style = "table_template")
      }

      n_flagged <- sum(calls$verdict %in% FLAGGED_VERDICTS)
      doc <- doc %>%
        body_add_par("Executive Summary", style = "heading 1") %>%
        body_add_par(sprintf(
          "This batch has %d sample(s). Of %d notable (species, sample) call(s) across the whole dataset (min %d reads AND %.2f%% of classified reads), %d are flagged and should not be reported as-is \u2014 see the Actions Required section below.",
          nrow(ov), nrow(calls), min_reads_v(), min_pct_v(), n_flagged
        ), style = "Normal")

      doc <- doc %>%
        body_add_par("Read Yield Overview", style = "heading 1") %>%
        body_add_par(paste0("Classified reads per sample. The dashed line marks the low-yield reference threshold (",
                             format(thr, big.mark = ","), " reads)."), style = "Normal") %>%
        body_add_img(chart_path, width = 6.5, height = 6.5 * 5 / 9)

      doc <- doc %>%
        body_add_par("Taxonomic Composition", style = "heading 1") %>%
        body_add_par(paste0("Composition by ", input$tax_rank, ", across all samples:"), style = "Normal") %>%
        body_add_img(comp_path, width = 6.5, height = 6.5 * 5 / 9)

      doc <- doc %>% body_add_par("Negative Control", style = "heading 1")
      ctrl_col <- cnts$control_col
      if (!is.na(ctrl_col)) {
        cs <- control_summary()
        cs_tab <- cs %>% mutate(Value = format(Value, big.mark = ",")) %>% as.data.frame()
        doc <- doc %>% body_add_table(cs_tab, style = "table_template")
        flagged_species <- calls %>% filter(control_reads > 0) %>% distinct(species, control_reads)
        doc <- doc %>% body_add_par(
          if (nrow(flagged_species) == 0) "No notable calls were found in the negative control."
          else paste0("Notable calls also present in the negative control: ",
                       paste(sprintf("%s (%d reads)", flagged_species$species, flagged_species$control_reads), collapse = "; "), "."),
          style = "Normal")
      }

      cs_all <- tryCatch(case_summary(), error = function(e) NULL)
      if (!is.null(cs_all) && nrow(cs_all) > 0) {
        doc <- doc %>%
          body_add_par("Case Summaries", style = "heading 1") %>%
          body_add_par("Overview of all GTX cases in this batch:", style = "Normal")
        case_tab <- cs_all %>%
          transmute(Case = case, N = n_samples, `Well-Pwr` = n_well_powered,
                    `Low-Yld` = n_low_yield, Calls = total_calls, Trust = n_trusted,
                    Single = n_singleton, `Blw Ctrl` = n_below_control, Hop = n_hop, Tier = tier) %>%
          as.data.frame()
        doc <- doc %>% body_add_table(case_tab, style = "table_template") %>%
          body_add_par(
            "N = samples; Well-Pwr / Low-Yld = well-powered / low-yield sample counts; Calls = total notable calls; Trust / Single / Blw Ctrl / Hop = Trusted, Singleton, Below Control, and Possible Index Hopping call counts; Tier = overall case confidence.",
            style = "Normal")

        sl_full_export <- sample_level()
        for (cs_name in cs_all$case) {
          crow <- cs_all %>% filter(case == cs_name)
          slc <- sl_full_export %>% filter(case == cs_name) %>% arrange(barcode)
          doc <- doc %>% body_add_par(sprintf("%s \u2014 Confidence: %s", cs_name, crow$tier), style = "heading 2")
          for (i in seq_len(nrow(slc))) {
            doc <- doc %>% body_add_par(paste0("\u2022  ", sample_narrative(slc[i, ], sl_full_export$pre_pcr, sl_full_export$post_pcr)), style = "Normal")
          }
          for (l in case_conclusion_text(crow, slc, sl_full_export$pre_pcr, sl_full_export$post_pcr)) {
            doc <- doc %>% body_add_par(l, style = "Normal")
          }
          removal_items_export <- case_removal_items(calls, cs_name)
          if (length(removal_items_export) > 0) {
            doc <- doc %>% body_add_par("To be removed from this case's report:", style = "Normal")
            for (r in removal_items_export) doc <- doc %>% body_add_par(paste0("\u2022  ", r), style = "Normal")
          }
          trusted_species_c <- slc$species_trusted[nzchar(slc$species_trusted)]
          trusted_species_c <- unique(unlist(str_split(paste(trusted_species_c, collapse = ", "), ",\\s*")))
          trusted_species_c <- trusted_species_c[nzchar(trusted_species_c)]
          doc <- doc %>% body_add_par(case_organism_comment(trusted_species_c), style = "Normal")
        }
      }

      doc <- doc %>% body_add_par("Actions Required Before Results Can Be Reliably Reported", style = "heading 1")
      flagged <- calls %>% filter(verdict %in% FLAGGED_VERDICTS)
      if (nrow(flagged) == 0) {
        doc <- doc %>% body_add_par("No blocking issues identified at the current notability thresholds.", style = "Normal")
      } else {
        for (i in seq_len(nrow(flagged))) {
          r <- flagged[i, ]
          txt <- if (r$verdict == "Below control") {
            sprintf("Exclude %s from %s \u2014 %d reads, at or below the negative control's %d reads for this species.",
                    r$species, r$label, r$reads, r$control_reads)
          } else if (r$verdict == "Possible index hopping") {
            sprintf("Verify %s in %s against the plate layout \u2014 possible index hopping from %s (%d reads), corroborated by %d/%d companion species.",
                    r$species, r$label, r$top_carrier_sample, r$top_carrier_reads, r$hop_n_matched, r$hop_n_checked)
          } else {
            sprintf("Low-priority review: %s in %s \u2014 weak/uncorroborated hop signal, not independently recurring; verify if certainty is required.",
                    r$species, r$label)
          }
          doc <- doc %>% body_add_par(paste0("\u2022  ", txt), style = "Normal")
        }
      }
      cs_low <- tryCatch(case_summary(), error = function(e) NULL)
      if (!is.null(cs_low)) {
        low_cases <- cs_low %>% filter(tier == "Low")
        for (i in seq_len(nrow(low_cases))) {
          r <- low_cases[i, ]
          doc <- doc %>% body_add_par(sprintf(
            "\u2022  Case %s: overall confidence Low (%d/%d samples well-powered) \u2014 review before reporting the case as a whole; consider repeating its weaker samples.",
            r$case, r$n_well_powered, r$n_samples), style = "Normal")
        }
      }

      doc <- doc %>% body_add_break() %>% body_add_par("Annexe: Methodology", style = "heading 1")
      for (para in methodology_paragraphs) doc <- add_markdown_paragraph(doc, para)

      print(doc, target = file)
    }
  )

}

shinyApp(ui, server)

#   shiny::runApp("app.R")
