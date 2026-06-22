#!/usr/bin/env Rscript
# Builds Word / CSV full model tables from exported GLMM stats CSVs.
# Run after AOC_stats_glmm_export.R (or run export script, which sources this file).

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

stats_dir <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/stats"
out_dir <- file.path(stats_dir, "FullModelTables")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

fmt_num <- function(x) {
  ifelse(is.na(x), "", sprintf("%.3f", as.numeric(x)))
}

fmt_p <- function(p) {
  p <- as.numeric(p)
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.001] <- "< .001"
  idx <- !is.na(p) & p >= 0.001
  out[idx] <- sub("^0", "", sprintf("%.3f", p[idx]))
  out
}

fmt_ci <- function(lo, hi) {
  paste0("[", sprintf("%.3f", as.numeric(lo)), ", ", sprintf("%.3f", as.numeric(hi)), "]")
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
}

merge_pairwise_effectsizes <- function(pair_df, eff_df, task_name, dv_name) {
  if (is.null(pair_df) || nrow(pair_df) == 0 || is.null(eff_df) || nrow(eff_df) == 0) {
    return(pair_df)
  }
  req_cols <- c("Task", "Variable", "Group1", "Group2", "Cohens_dz")
  if (!all(req_cols %in% names(eff_df))) return(pair_df)

  p <- pair_df
  p$Task <- if ("Task" %in% names(p)) as.character(p$Task) else task_name
  p$Variable <- if ("Variable" %in% names(p)) as.character(p$Variable) else if ("AlphaDV" %in% names(p)) as.character(p$AlphaDV) else dv_name
  p$GroupA <- pmin(as.character(p$Group1), as.character(p$Group2))
  p$GroupB <- pmax(as.character(p$Group1), as.character(p$Group2))
  p$.pair_order <- seq_len(nrow(p))

  e <- eff_df
  e <- e[e$Task == task_name & e$Variable == dv_name, , drop = FALSE]
  if (nrow(e) == 0) return(pair_df)
  e$GroupA <- pmin(as.character(e$Group1), as.character(e$Group2))
  e$GroupB <- pmax(as.character(e$Group1), as.character(e$Group2))
  e <- e[, c("GroupA", "GroupB", "Cohens_dz"), drop = FALSE]

  merged <- merge(p, e, by = c("GroupA", "GroupB"), all.x = TRUE, sort = FALSE)
  merged <- merged[order(merged$.pair_order), , drop = FALSE]
  merged$GroupA <- NULL
  merged$GroupB <- NULL
  merged$.pair_order <- NULL
  merged
}

display_dv <- function(dv) {
  out <- as.character(dv)
  out[out == "ERSD"] <- "ERS/ERD"
  out
}

prettify_model_label <- function(model_label, dv_label = NULL, predictor_label = NULL) {
  x <- as.character(model_label)
  x <- gsub("^DV\\s*~", paste0(dv_label, " ~"), x)
  x <- gsub("\\(1\\|ID\\)", "(1|Subject)", x)
  x <- gsub("\\(1\\|SubjectID\\)", "(1|Subject)", x)
  x <- gsub("_c\\b", "", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  if (!is.null(predictor_label) && nzchar(predictor_label)) {
    x <- gsub(predictor_label, predictor_label, x, fixed = TRUE)
  }
  x
}

clean_term <- function(term, load_ref = NULL) {
  term <- as.character(term)
  out <- term
  out[out == "Intercept"] <- "Intercept"
  out[out == "Group Var"] <- "Random intercept variance (Subject)"
  out[out == "Load"] <- "Load"
  out[out == "Interaction"] <- "Interaction term"
  out[out == "GazeDeviation_c"] <- "Gaze deviation (centered)"
  out[out == "MSRate_c"] <- "Microsaccade rate (centered)"
  out[out == "GazeDeviationBL_c"] <- "Gaze deviation baseline (centered)"
  out[out == "MSRateBL_c"] <- "Microsaccade rate baseline (centered)"

  # lme4 Load main effects (treatment coding)
  load_main <- grepl("^Load", out) & !grepl(":", out)
  if (any(load_main)) {
    lvls <- sub("^Load", "", out[load_main])
    ref_txt <- if (!is.null(load_ref) && nzchar(load_ref)) load_ref else "reference"
    out[load_main] <- paste0("Load: ", lvls, " vs ", ref_txt)
  }

  # lme4 gaze:Load interactions
  int_lme4 <- grepl("_c:Load|^Load.*_c:", out)
  if (any(int_lme4)) {
    for (i in which(int_lme4)) {
      parts <- strsplit(out[i], ":", fixed = TRUE)[[1]]
      pred <- parts[1]
      load_part <- parts[2]
      lvl <- sub("^Load", "", load_part)
      pred <- gsub("_c$", "", pred)
      pred <- switch(
        pred,
        GazeDeviation = "Gaze deviation (centered)",
        MSRate = "Microsaccade rate (centered)",
        GazeDeviationBL = "Gaze deviation baseline (centered)",
        MSRateBL = "Microsaccade rate baseline (centered)",
        pred
      )
      ref_txt <- if (!is.null(load_ref) && nzchar(load_ref)) load_ref else "reference"
      out[i] <- paste0("Interaction: ", pred, " x Load (", lvl, " vs ", ref_txt, ")")
    }
  }

  # statsmodels legacy Condition terms (fallback)
  cond_pattern <- "^C\\(Condition, Treatment\\(reference=\"([^\"]+)\"\\)\\)\\[T\\.(.+)\\]$"
  cond_hits <- grepl(cond_pattern, out)
  if (any(cond_hits)) {
    refs <- sub(cond_pattern, "\\1", out[cond_hits])
    lvls <- sub(cond_pattern, "\\2", out[cond_hits])
    out[cond_hits] <- paste0("Condition: ", lvls, " vs ", refs)
  }

  int_pattern <- "^([A-Za-z0-9_]+):C\\(Condition, Treatment\\(reference=\"([^\"]+)\"\\)\\)\\[T\\.(.+)\\]$"
  int_hits <- grepl(int_pattern, out)
  if (any(int_hits)) {
    pred <- sub(int_pattern, "\\1", out[int_hits])
    ref <- sub(int_pattern, "\\2", out[int_hits])
    lvl <- sub(int_pattern, "\\3", out[int_hits])
    pred <- ifelse(pred == "GazeDeviation_c", "Gaze deviation (centered)", pred)
    pred <- ifelse(pred == "MSRate_c", "Microsaccade rate (centered)", pred)
    out[int_hits] <- paste0("Interaction: ", pred, " x Condition (", lvl, " vs ", ref, ")")
  }
  out
}

clean_fixed <- function(df, load_ref = NULL) {
  out <- df
  out$Term <- clean_term(out$Term, load_ref = load_ref)
  out$`β` <- fmt_num(out$beta)
  out$SE <- fmt_num(out$SE)
  out$Statistic <- fmt_num(out$stat)
  out$`p-value` <- fmt_p(out$p)
  out$CI <- fmt_ci(out$CI_low, out$CI_high)
  out <- out[, c("Term", "β", "SE", "Statistic", "p-value", "CI")]
  names(out)[1] <- "Term"
  out
}

clean_random <- function(df, load_ref = NULL) {
  out <- df
  out$Term <- clean_term(out$Term, load_ref = load_ref)
  out$Variance <- fmt_num(out$beta)
  out$SE <- fmt_num(out$SE)
  out$Statistic <- fmt_num(out$stat)
  out$`p-value` <- fmt_p(out$p)
  out$CI <- fmt_ci(out$CI_low, out$CI_high)
  out <- out[, c("Term", "Variance", "SE", "Statistic", "p-value", "CI")]
  names(out)[1] <- "Term"
  out
}

clean_lrt <- function(df, load_ref = NULL) {
  out <- df
  out$Term <- clean_term(out$Term, load_ref = load_ref)
  out$df <- fmt_num(out$df)
  out$`χ²` <- fmt_num(out$LR_Chi2)
  out$`p-value` <- fmt_p(out$p)
  out$`R2 (LR)` <- fmt_num(out$R2_LR)
  out$`f2 (LR)` <- fmt_num(out$f2_LR)
  out[, c("Term", "df", "χ²", "p-value", "R2 (LR)", "f2 (LR)")]
}

clean_pairwise <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame(
      Contrast = character(),
      `β` = character(),
      SE = character(),
      Statistic = character(),
      padj = character(),
      CI = character(),
      `Cohen's d` = character(),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  out <- df
  # Estimate is Group2 - Group1; label as Group2 vs Group1 (matches fixed effects)
  out$Contrast <- paste(out$Group2, "vs", out$Group1)
  out$`β` <- fmt_num(out$Estimate)
  out$SE <- fmt_num(out$SE)
  out$Statistic <- fmt_num(out$z)
  if ("p_adj" %in% names(out)) {
    out$padj <- fmt_p(out$p_adj)
  } else {
    out$padj <- fmt_p(out$p)
  }
  if ("CI95_low" %in% names(out)) {
    out$CI <- fmt_ci(out$CI95_low, out$CI95_high)
  } else {
    out$CI <- fmt_ci(out$CI_low, out$CI_high)
  }
  if (!("Cohens_dz" %in% names(out))) out$Cohens_dz <- NA_real_
  out$`Cohen's d` <- fmt_num(out$Cohens_dz)
  out[, c("Contrast", "β", "SE", "Statistic", "padj", "CI", "Cohen's d")]
}

style_ft <- function(ft, col_widths = NULL) {
  ft <- theme_booktabs(ft)
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- bold(ft, part = "header")
  ft <- padding(ft, padding = 1, part = "all")
  ft <- align(ft, align = "left", j = 1, part = "all")
  if (ncol(ft$body$dataset) > 1) {
    ft <- align(ft, align = "center", j = 2:ncol(ft$body$dataset), part = "all")
  }
  ft <- set_table_properties(ft, layout = "fixed")
  ft <- line_spacing(ft, space = 0.9, part = "all")
  ft <- hrule(ft, rule = "exact", part = "all")
  ft <- height_all(ft, height = 0.18, part = "all")
  if (!is.null(col_widths) && length(col_widths) == ncol(ft$body$dataset)) {
    ft <- width(ft, j = seq_along(col_widths), width = col_widths)
  } else {
    ft <- autofit(ft)
  }
  ft
}

doc_add_title <- function(doc, text, size = 10) {
  body_add_fpar(
    doc,
    fpar(
      ftext(text, fp_text(bold = TRUE, font.size = size))
    )
  )
}

doc_add_table <- function(doc, title, df, table_type = c("fixed", "random", "lrt", "pairwise")) {
  table_type <- match.arg(table_type)
  ft <- flextable(df)
  if (table_type %in% c("fixed", "random")) {
    ft <- style_ft(ft, col_widths = c(3.6, 0.9, 0.8, 1.0, 1.0, 1.6))
  } else if (table_type == "lrt") {
    ft <- style_ft(ft, col_widths = c(2.5, 0.7, 1.0, 1.0, 0.9, 0.8))
  } else {
    ft <- style_ft(ft, col_widths = c(2.2, 0.9, 0.8, 1.0, 0.95, 1.2, 0.95))
  }
  doc <- doc_add_title(doc, title, size = 10)
  body_add_flextable(doc, ft)
}

pair_nback <- safe_read_csv(file.path(stats_dir, "AOC_pairwise_mixedlm_nback.csv"))
pair_stern <- safe_read_csv(file.path(stats_dir, "AOC_pairwise_mixedlm_sternberg.csv"))
if (is.null(pair_nback)) pair_nback <- data.frame()
if (is.null(pair_stern)) pair_stern <- data.frame()
eff_nback <- safe_read_csv(file.path(stats_dir, "AOC_pairwise_effectsizes_nback.csv"))
eff_stern <- safe_read_csv(file.path(stats_dir, "AOC_pairwise_effectsizes_sternberg.csv"))

load_ref_for_task <- function(task_name) {
  if (tolower(task_name) == "nback") "1-back" else "WM load 2"
}

model_specs <- list(
  list(
    id = "nback_accuracy",
    title = "N-back Accuracy",
    fixed = "AOC_mixedlm_fixed_acc_nback.csv",
    pair_src = if (nrow(pair_nback) > 0) pair_nback[pair_nback$Variable == "Accuracy", ] else NULL
  ),
  list(
    id = "sternberg_accuracy",
    title = "Sternberg Accuracy",
    fixed = "AOC_mixedlm_fixed_acc_sternberg.csv",
    pair_src = if (nrow(pair_stern) > 0) pair_stern[pair_stern$Variable == "Accuracy", ] else NULL
  ),
  list(
    id = "nback_reactiontime",
    title = "N-back Reaction Time",
    fixed = "AOC_mixedlm_fixed_rt_nback.csv",
    pair_src = if (nrow(pair_nback) > 0) pair_nback[pair_nback$Variable == "ReactionTime", ] else NULL
  ),
  list(
    id = "sternberg_reactiontime",
    title = "Sternberg Reaction Time",
    fixed = "AOC_mixedlm_fixed_rt_sternberg.csv",
    pair_src = if (nrow(pair_stern) > 0) pair_stern[pair_stern$Variable == "ReactionTime", ] else NULL
  ),
  list(
    id = "nback_gazedeviation",
    title = "N-back Gaze Deviation",
    fixed = "AOC_mixedlm_fixed_gazedev_nback.csv",
    pair_src = if (nrow(pair_nback) > 0) pair_nback[pair_nback$Variable == "GazeDeviation", ] else NULL
  ),
  list(
    id = "sternberg_gazedeviation",
    title = "Sternberg Gaze Deviation",
    fixed = "AOC_mixedlm_fixed_gazedev_sternberg.csv",
    pair_src = if (nrow(pair_stern) > 0) pair_stern[pair_stern$Variable == "GazeDeviation", ] else NULL
  ),
  list(
    id = "nback_msrate",
    title = "N-back Microsaccade Rate",
    fixed = "AOC_mixedlm_fixed_ms_nback.csv",
    pair_src = if (nrow(pair_nback) > 0) pair_nback[pair_nback$Variable == "MSRate", ] else NULL
  ),
  list(
    id = "sternberg_msrate",
    title = "Sternberg Microsaccade Rate",
    fixed = "AOC_mixedlm_fixed_ms_sternberg.csv",
    pair_src = if (nrow(pair_stern) > 0) pair_stern[pair_stern$Variable == "MSRate", ] else NULL
  ),
  list(
    id = "nback_ersd",
    title = "N-back ERS/ERD",
    fixed = "AOC_mixedlm_fixed_ersd_nback.csv",
    pair_src = if (nrow(pair_nback) > 0) pair_nback[pair_nback$Variable == "ERSD", ] else NULL
  ),
  list(
    id = "sternberg_ersd",
    title = "Sternberg ERS/ERD",
    fixed = "AOC_mixedlm_fixed_ersd_sternberg.csv",
    pair_src = if (nrow(pair_stern) > 0) pair_stern[pair_stern$Variable == "ERSD", ] else NULL
  ),
  list(
    id = "nback_ersd_by_gazedeviation",
    title = "N-back ERS/ERD by GazeDeviation",
    fixed = "AOC_ersd_fixed_gazedeviation_nback.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_nback.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_nback.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedeviation",
    title = "Sternberg ERS/ERD by GazeDeviation",
    fixed = "AOC_ersd_fixed_gazedeviation_sternberg.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_sternberg.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_sternberg.csv"
  ),
  list(
    id = "nback_ersd_by_msrate",
    title = "N-back ERS/ERD by MSRate",
    fixed = "AOC_ersd_fixed_msrate_nback.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_nback.csv",
    pair = "AOC_ersd_contrasts_msrate_nback.csv"
  ),
  list(
    id = "sternberg_ersd_by_msrate",
    title = "Sternberg ERS/ERD by MSRate",
    fixed = "AOC_ersd_fixed_msrate_sternberg.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_sternberg.csv",
    pair = "AOC_ersd_contrasts_msrate_sternberg.csv"
  ),
  list(
    id = "nback_ersd_by_gazedev_bl",
    title = "N-back ERS/ERD by GazeDeviationBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_nback.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_nback.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedev_bl",
    title = "Sternberg ERS/ERD by GazeDeviationBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_sternberg.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_sternberg.csv"
  ),
  list(
    id = "nback_ersd_by_ms_bl",
    title = "N-back ERS/ERD by MSRateBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_nback.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_nback.csv"
  ),
  list(
    id = "sternberg_ersd_by_ms_bl",
    title = "Sternberg ERS/ERD by MSRateBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_sternberg.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_sternberg.csv"
  )
)

for (spec in model_specs) {
  fixed_path <- file.path(stats_dir, spec$fixed)
  if (!file.exists(fixed_path)) {
    message("Skipping ", spec$id, " (missing ", spec$fixed, ")")
    next
  }
  fixed_raw <- read.csv(fixed_path, stringsAsFactors = FALSE)
  task_name <- tolower(as.character(unique(fixed_raw$Task)[1]))
  load_ref <- load_ref_for_task(task_name)
  random_mask <- grepl("Var$|Cov$|Group Var", fixed_raw$Term)
  fixed_tbl <- clean_fixed(fixed_raw[!random_mask, , drop = FALSE], load_ref = load_ref)
  random_tbl <- NULL
  if (any(random_mask)) {
    random_tbl <- clean_random(fixed_raw[random_mask, , drop = FALSE], load_ref = load_ref)
  }

  fixed_csv_path <- file.path(out_dir, paste0(spec$id, "_fixed.csv"))
  write.csv(fixed_tbl, fixed_csv_path, row.names = FALSE)
  if (!is.null(random_tbl)) {
    write.csv(random_tbl, file.path(out_dir, paste0(spec$id, "_random.csv")), row.names = FALSE)
  }

  pair_tbl <- NULL
  pair_raw <- NULL
  dv_name <- as.character(unique(fixed_raw$DV)[1])
  eff_df <- if (!is.na(task_name) && task_name == "nback") eff_nback else eff_stern
  if (!is.null(spec$pair_src)) {
    pair_raw <- spec$pair_src
  } else if (!is.null(spec$pair)) {
    pair_raw <- read.csv(file.path(stats_dir, spec$pair), stringsAsFactors = FALSE)
  }
  if (!is.null(pair_raw) && nrow(pair_raw) > 0) {
    pair_raw <- merge_pairwise_effectsizes(
      pair_raw,
      eff_df = eff_df,
      task_name = task_name,
      dv_name = dv_name
    )
    pair_tbl <- clean_pairwise(pair_raw)
  } else if (!is.null(pair_raw) && nrow(pair_raw) == 0) {
    pair_tbl <- NULL
    message(
      "Skipping pairwise for ", spec$id,
      " (no rows in pairwise source; regenerate AOC_pairwise_* when Variable names match)."
    )
  }
  if (!is.null(pair_tbl) && nrow(pair_tbl) > 0) {
    write.csv(pair_tbl, file.path(out_dir, paste0(spec$id, "_pairwise.csv")), row.names = FALSE)
  }

  lrt_tbl <- NULL
  if (!is.null(spec$lrt)) {
    lrt_path <- file.path(stats_dir, spec$lrt)
    if (file.exists(lrt_path)) {
      lrt_raw <- read.csv(lrt_path, stringsAsFactors = FALSE)
      lrt_tbl <- clean_lrt(lrt_raw, load_ref = load_ref)
      write.csv(lrt_tbl, file.path(out_dir, paste0(spec$id, "_lrt.csv")), row.names = FALSE)
    }
  }

  doc <- read_docx()
  doc <- doc_add_title(doc, spec$title, size = 10)
  formula_text <- unique(fixed_raw$ModelLabel)
  dv_text <- unique(fixed_raw$DV)
  dv_text <- if (length(dv_text) > 0 && !is.na(dv_text[1])) display_dv(dv_text[1]) else "Outcome"
  if (length(formula_text) > 0 && !is.na(formula_text[1])) {
    model_text <- prettify_model_label(formula_text[1], dv_label = dv_text)
    doc <- body_add_fpar(
      doc,
      fpar(ftext(paste0("Model: ", model_text), fp_text(font.size = 10)))
    )
  }
  doc <- body_add_par(doc, "")
  doc <- doc_add_table(doc, "Fixed effects", fixed_tbl, table_type = "fixed")
  if (!is.null(random_tbl)) {
    doc <- body_add_par(doc, "")
    doc <- doc_add_table(doc, "Random effects", random_tbl, table_type = "random")
  }

  if (!is.null(lrt_tbl)) {
    doc <- body_add_par(doc, "")
    doc <- doc_add_table(doc, "Likelihood-ratio tests", lrt_tbl, table_type = "lrt")
  }

  if (!is.null(pair_tbl) && nrow(pair_tbl) > 0) {
    doc <- body_add_par(doc, "")
    doc <- doc_add_table(doc, "Pairwise contrasts", pair_tbl, table_type = "pairwise")
  }

  print(doc, target = file.path(out_dir, paste0(spec$id, "_full_model_table.docx")))
}

message("Full model CSVs and DOCX tables generated in: ", out_dir)
