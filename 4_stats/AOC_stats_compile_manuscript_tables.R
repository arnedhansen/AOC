#!/usr/bin/env Rscript
# Compile manuscript in-text and supplementary full-model table documents.

suppressPackageStartupMessages({
  library(officer)
  library(xml2)
})

script_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
}
source(file.path(script_dir, "AOC_stats_glmm_helpers.R"))

stats_dir <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/stats"
combined_dir <- file.path(stats_dir, "combinedTables")
full_dir <- file.path(stats_dir, "fullModelTables")
full_trials_dir <- file.path(stats_dir, "fullModelTables_trials")
manuscript_dir <- file.path(stats_dir, "manuscript")

dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

resolve_src <- function(folder, filename) {
  path <- file.path(folder, filename)
  if (!file.exists(path)) {
    stop("Missing table source: ", path)
  }
  path
}

strip_docx_title_paragraph <- function(src_path) {
  work <- tempfile()
  dir.create(work)
  out <- tempfile(fileext = ".docx")
  utils::unzip(src_path, exdir = work)
  doc_xml_path <- file.path(work, "word", "document.xml")
  if (!file.exists(doc_xml_path)) {
    unlink(work, recursive = TRUE)
    return(src_path)
  }
  xml <- read_xml(doc_xml_path)
  body <- xml_find_first(xml, "/*[local-name()='document']/*[local-name()='body']")
  first_p <- xml_find_first(body, "*[local-name()='p'][1]")
  if (!is.na(first_p)) {
    xml_remove(first_p)
    write_xml(xml, doc_xml_path)
  }
  owd <- getwd()
  on.exit({
    setwd(owd)
    unlink(work, recursive = TRUE)
  }, add = TRUE)
  setwd(work)
  status <- system2("zip", c("-qr", out, "."), stdout = FALSE, stderr = FALSE)
  if (!identical(status, 0L) || !file.exists(out)) {
    return(src_path)
  }
  out
}

compile_manuscript_docx <- function(specs, out_path, doc_title, strip_embedded_titles = FALSE) {
  master <- read_docx()
  master <- body_add_fpar(
    master,
    fpar(ftext(doc_title, fp_text(bold = TRUE, font.size = 12)))
  )
  temp_files <- character()
  for (spec in specs) {
    src <- resolve_src(spec$folder, spec$file)
    embed_src <- if (isTRUE(strip_embedded_titles)) {
      stripped <- strip_docx_title_paragraph(src)
      if (!identical(stripped, src)) {
        temp_files <- c(temp_files, stripped)
      }
      stripped
    } else {
      src
    }
    master <- body_add_par(master, "")
    master <- body_add_fpar(
      master,
      fpar(ftext(spec$label, fp_text(bold = TRUE, font.size = 11)))
    )
    master <- body_add_docx(master, src = embed_src)
  }
  print(master, target = out_path)
  if (length(temp_files) > 0) {
    unlink(temp_files)
  }
  message("Wrote ", out_path)
}

build_numbered_specs <- function(entries, label_fmt, index_offset = 0L) {
  included <- Filter(function(entry) {
    file.exists(file.path(entry$folder, entry$file))
  }, entries)
  Map(function(i, entry) {
    list(
      label = sprintf(label_fmt, i + index_offset, entry$title),
      folder = entry$folder,
      file = entry$file
    )
  }, seq_along(included), included)
}

main_text_specs <- build_numbered_specs(
  list(
    list(title = "Reaction Time", folder = combined_dir, file = "reaction_time_results_table.docx"),
    list(title = "Accuracy", folder = combined_dir, file = "accuracy_results_table.docx"),
    list(title = "Gaze Deviation", folder = combined_dir, file = "gaze_deviation_results_table.docx"),
    list(title = "Microsaccade Rate", folder = combined_dir, file = "microsaccade_rate_results_table.docx"),
    list(title = "ERS/ERD", folder = combined_dir, file = "ers_erd_results_table.docx"),
    list(
      title = "ERS/ERD Co-Variation with Gaze Deviation",
      folder = combined_dir,
      file = "ersd_gaze_deviation_reduced_results_table.docx"
    ),
    list(
      title = "ERS/ERD Co-Variation with Microsaccade Rate",
      folder = combined_dir,
      file = "ersd_microsaccade_rate_reduced_results_table.docx"
    )
  ),
  "Table %d: %s",
  index_offset = 1L
)

supplement_entries <- list(
  list(title = "Reaction Time", folder = full_dir, file = "reaction_time_full_model_table.docx"),
  list(title = "Accuracy", folder = full_dir, file = "accuracy_full_model_table.docx"),
  list(title = "Gaze Deviation", folder = full_dir, file = "gaze_deviation_full_model_table.docx"),
  list(title = "Microsaccade Rate", folder = full_dir, file = "microsaccade_rate_full_model_table.docx"),
  list(title = "ERS/ERD", folder = full_dir, file = "ers_erd_full_model_table.docx"),
  list(
    title = "ERS/ERD Co-Variation with Gaze Deviation",
    folder = full_dir,
    file = "ersd_gaze_deviation_reduced_full_model_table.docx"
  ),
  list(
    title = "ERS/ERD Co-Variation with Microsaccade Rate",
    folder = full_dir,
    file = "ersd_microsaccade_rate_reduced_full_model_table.docx"
  ),
  list(
    title = "Baselined Gaze Deviation",
    folder = full_dir,
    file = "gaze_deviation_bl_full_model_table.docx"
  ),
  list(
    title = "Baselined Microsaccade Rate",
    folder = full_dir,
    file = "microsaccade_rate_bl_full_model_table.docx"
  ),
  list(
    title = "ERS/ERD Co-Variation with Baselined Gaze Deviation",
    folder = full_dir,
    file = "ersd_gaze_deviation_bl_reduced_full_model_table.docx"
  ),
  list(
    title = "ERS/ERD Co-Variation with Baselined Microsaccade Rate",
    folder = full_dir,
    file = "ersd_microsaccade_rate_bl_reduced_full_model_table.docx"
  ),
  list(
    title = "Baseline-Window Alpha Power",
    folder = full_dir,
    file = "alpha_power_baselineWindow_full_model_table.docx"
  ),
  list(
    title = "ERS/ERD Co-Variation with Baselined Gaze Deviation",
    folder = full_trials_dir,
    file = "ersd_gaze_deviation_bl_reduced_full_model_table.docx"
  ),
  list(
    title = "ERS/ERD Co-Variation with Baselined Microsaccade Rate",
    folder = full_trials_dir,
    file = "ersd_microsaccade_rate_bl_full_full_model_table.docx"
  )
)

supplement_specs <- build_numbered_specs(supplement_entries, "Table S3.%d: %s")

compile_manuscript_docx(
  main_text_specs,
  file.path(manuscript_dir, "manuscript_intext_tables_TEMPLATE.docx"),
  "Manuscript In-Text Tables (Tables 2–8)",
  strip_embedded_titles = TRUE
)

compile_manuscript_docx(
  supplement_specs,
  file.path(manuscript_dir, "manuscript_supplement_tables_TEMPLATE.docx"),
  "Supplementary Full Model Tables (S3)",
  strip_embedded_titles = TRUE
)
