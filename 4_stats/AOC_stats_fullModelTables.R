#!/usr/bin/env Rscript

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
  out[out == "AlphaPower_FOOOF_bl"] <- "AlphaPowerSpecParamBaseline"
  out[out == "AlphaPower"] <- "AlphaPower"
  out
}

prettify_model_label <- function(model_label, dv_label = NULL, predictor_label = NULL) {
  x <- as.character(model_label)
  x <- gsub("^DV\\s*~", paste0(dv_label, " ~"), x)
  x <- gsub("\\(1\\|ID\\)", "(1|SubjectID)", x)
  x <- gsub("_c\\b", "", x)
  x <- gsub("AlphaPower_FOOOF_bl", "AlphaPowerSpecParamBaseline", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  if (!is.null(predictor_label) && nzchar(predictor_label)) {
    x <- gsub(predictor_label, predictor_label, x, fixed = TRUE)
  }
  x
}

clean_term <- function(term) {
  term <- as.character(term)
  out <- term
  out[out == "Intercept"] <- "Intercept"
  out[out == "Group Var"] <- "Random intercept variance"
  out[out == "Condition"] <- "Condition"
  out[out == "Interaction"] <- "Interaction term"
  out[out == "GazeDeviation_c"] <- "Gaze deviation (centered)"
  out[out == "MSRate_c"] <- "Microsaccade rate (centered)"

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

clean_fixed <- function(df) {
  out <- df
  out$Term <- clean_term(out$Term)
  out$`β` <- fmt_num(out$beta)
  out$SE <- fmt_num(out$SE)
  out$`Statistic (z/t)` <- fmt_num(out$stat)
  out$`p-value` <- fmt_p(out$p)
  out$CI <- fmt_ci(out$CI_low, out$CI_high)
  out <- out[, c("Term", "β", "SE", "Statistic (z/t)", "p-value", "CI")]
  names(out)[1] <- "Term"
  out
}

clean_lrt <- function(df) {
  out <- df
  out$Term <- clean_term(out$Term)
  out$df <- fmt_num(out$df)
  out$`χ²` <- fmt_num(out$LR_Chi2)
  out$`p-value` <- fmt_p(out$p)
  out$`R2 (LR)` <- fmt_num(out$R2_LR)
  out$`f2 (LR)` <- fmt_num(out$f2_LR)
  out[, c("Term", "df", "χ²", "p-value", "R2 (LR)", "f2 (LR)")]
}

clean_pairwise <- function(df) {
  out <- df
  out$Contrast <- paste(out$Group1, "vs", out$Group2)
  out$`β` <- fmt_num(out$Estimate)
  out$SE <- fmt_num(out$SE)
  out$`Statistic (z)` <- fmt_num(out$z)
  out$`p-value` <- fmt_p(out$p)
  out$`p_adj` <- fmt_p(out$p_adj)
  if ("CI95_low" %in% names(out)) {
    out$CI <- fmt_ci(out$CI95_low, out$CI95_high)
  } else {
    out$CI <- fmt_ci(out$CI_low, out$CI_high)
  }
  if (!("Cohens_dz" %in% names(out))) out$Cohens_dz <- NA_real_
  out$`Cohen's d` <- fmt_num(out$Cohens_dz)
  out[, c("Contrast", "β", "SE", "Statistic (z)", "p-value", "p_adj", "CI", "Cohen's d")]
}

style_ft <- function(ft, col_widths = NULL) {
  ft <- theme_booktabs(ft)
  ft <- fontsize(ft, size = 10, part = "all")
  ft <- bold(ft, part = "header")
  ft <- padding(ft, padding = 4, part = "all")
  ft <- align(ft, align = "left", j = 1, part = "all")
  if (ncol(ft$body$dataset) > 1) {
    ft <- align(ft, align = "center", j = 2:ncol(ft$body$dataset), part = "all")
  }
  ft <- set_table_properties(ft, layout = "fixed")
  ft <- line_spacing(ft, space = 1, part = "all")
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

doc_add_table <- function(doc, title, df, table_type = c("fixed", "lrt", "pairwise")) {
  table_type <- match.arg(table_type)
  ft <- flextable(df)
  if (table_type == "fixed") {
    ft <- style_ft(ft, col_widths = c(3.6, 0.9, 0.8, 1.0, 1.0, 1.6))
  } else if (table_type == "lrt") {
    ft <- style_ft(ft, col_widths = c(2.5, 0.7, 1.0, 1.0, 0.9, 0.8))
  } else {
    ft <- style_ft(ft, col_widths = c(2.0, 0.85, 0.75, 0.95, 0.85, 0.85, 1.2, 0.9))
    if ("p_adj" %in% names(df)) {
      ft <- compose(
        ft,
        part = "header",
        j = "p_adj",
        value = as_paragraph(
          as_chunk("p"),
          as_chunk("adj", props = fp_text(vertical.align = "subscript"))
        )
      )
    }
  }
  doc <- doc_add_title(doc, title, size = 10)
  body_add_flextable(doc, ft)
}

pair_nback <- read.csv(file.path(stats_dir, "AOC_pairwise_mixedlm_nback.csv"), stringsAsFactors = FALSE)
pair_stern <- read.csv(file.path(stats_dir, "AOC_pairwise_mixedlm_sternberg.csv"), stringsAsFactors = FALSE)
eff_nback <- safe_read_csv(file.path(stats_dir, "AOC_pairwise_effectsizes_nback.csv"))
eff_stern <- safe_read_csv(file.path(stats_dir, "AOC_pairwise_effectsizes_sternberg.csv"))

model_specs <- list(
  list(
    id = "nback_accuracy",
    title = "N-back Accuracy",
    fixed = "AOC_mixedlm_fixed_acc_nback.csv",
    pair_src = pair_nback[pair_nback$Variable == "Accuracy", ]
  ),
  list(
    id = "sternberg_accuracy",
    title = "Sternberg Accuracy",
    fixed = "AOC_mixedlm_fixed_acc_sternberg.csv",
    pair_src = pair_stern[pair_stern$Variable == "Accuracy", ]
  ),
  list(
    id = "nback_reactiontime",
    title = "N-back Reaction Time",
    fixed = "AOC_mixedlm_fixed_rt_nback.csv",
    pair_src = pair_nback[pair_nback$Variable == "ReactionTime", ]
  ),
  list(
    id = "sternberg_reactiontime",
    title = "Sternberg Reaction Time",
    fixed = "AOC_mixedlm_fixed_rt_sternberg.csv",
    pair_src = pair_stern[pair_stern$Variable == "ReactionTime", ]
  ),
  list(
    id = "nback_gazedeviation",
    title = "N-back Gaze Deviation",
    fixed = "AOC_mixedlm_fixed_gazedev_nback.csv",
    pair_src = pair_nback[pair_nback$Variable == "GazeDeviation", ]
  ),
  list(
    id = "sternberg_gazedeviation",
    title = "Sternberg Gaze Deviation",
    fixed = "AOC_mixedlm_fixed_gazedev_sternberg.csv",
    pair_src = pair_stern[pair_stern$Variable == "GazeDeviation", ]
  ),
  list(
    id = "nback_msrate",
    title = "N-back Microsaccade Rate",
    fixed = "AOC_mixedlm_fixed_ms_nback.csv",
    pair_src = pair_nback[pair_nback$Variable == "MSRate", ]
  ),
  list(
    id = "sternberg_msrate",
    title = "Sternberg Microsaccade Rate",
    fixed = "AOC_mixedlm_fixed_ms_sternberg.csv",
    pair_src = pair_stern[pair_stern$Variable == "MSRate", ]
  ),
  list(
    id = "nback_alphapower",
    title = "N-back Alpha Power",
    fixed = "AOC_mixedlm_fixed_pow_nback.csv",
    pair_src = pair_nback[pair_nback$Variable == "AlphaPower", ]
  ),
  list(
    id = "sternberg_alphapower",
    title = "Sternberg Alpha Power",
    fixed = "AOC_mixedlm_fixed_pow_sternberg.csv",
    pair_src = pair_stern[pair_stern$Variable == "AlphaPower", ]
  ),
  list(
    id = "nback_alphapower_fooof_bl",
    title = "N-back AlphaPower specParam Baseline",
    fixed = "AOC_mixedlm_fixed_pow_specparam_bl_nback.csv",
    pair_src = pair_nback[pair_nback$Variable == "AlphaPower_FOOOF_bl", ]
  ),
  list(
    id = "sternberg_alphapower_fooof_bl",
    title = "Sternberg AlphaPower specParam Baseline",
    fixed = "AOC_mixedlm_fixed_pow_specparam_bl_sternberg.csv",
    pair_src = pair_stern[pair_stern$Variable == "AlphaPower_FOOOF_bl", ]
  ),
  list(
    id = "nback_alpha_by_gazedeviation",
    title = "N-back AlphaPower by GazeDeviation",
    fixed = "AOC_alpha_fixed_alphapower_gazedeviation_nback.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_gazedeviation_nback.csv",
    pair = "AOC_alpha_contrasts_alphapower_gazedeviation_nback.csv"
  ),
  list(
    id = "sternberg_alpha_by_gazedeviation",
    title = "Sternberg AlphaPower by GazeDeviation",
    fixed = "AOC_alpha_fixed_alphapower_gazedeviation_sternberg.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_gazedeviation_sternberg.csv",
    pair = "AOC_alpha_contrasts_alphapower_gazedeviation_sternberg.csv"
  ),
  list(
    id = "nback_alpha_by_msrate",
    title = "N-back AlphaPower by MSRate",
    fixed = "AOC_alpha_fixed_alphapower_msrate_nback.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_msrate_nback.csv",
    pair = "AOC_alpha_contrasts_alphapower_msrate_nback.csv"
  ),
  list(
    id = "sternberg_alpha_by_msrate",
    title = "Sternberg AlphaPower by MSRate",
    fixed = "AOC_alpha_fixed_alphapower_msrate_sternberg.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_msrate_sternberg.csv",
    pair = "AOC_alpha_contrasts_alphapower_msrate_sternberg.csv"
  )
)

for (spec in model_specs) {
  fixed_raw <- read.csv(file.path(stats_dir, spec$fixed), stringsAsFactors = FALSE)
  fixed_tbl <- clean_fixed(fixed_raw)

  fixed_csv_path <- file.path(out_dir, paste0(spec$id, "_fixed.csv"))
  write.csv(fixed_tbl, fixed_csv_path, row.names = FALSE)

  pair_tbl <- NULL
  pair_raw <- NULL
  task_name <- tolower(as.character(unique(fixed_raw$Task)[1]))
  dv_name <- as.character(unique(fixed_raw$DV)[1])
  eff_df <- if (!is.na(task_name) && task_name == "nback") eff_nback else eff_stern
  if (!is.null(spec$pair_src)) {
    pair_raw <- spec$pair_src
  } else if (!is.null(spec$pair)) {
    pair_raw <- read.csv(file.path(stats_dir, spec$pair), stringsAsFactors = FALSE)
  }
  if (!is.null(pair_raw)) {
    pair_raw <- merge_pairwise_effectsizes(
      pair_raw,
      eff_df = eff_df,
      task_name = task_name,
      dv_name = dv_name
    )
    pair_tbl <- clean_pairwise(pair_raw)
  }
  if (!is.null(pair_tbl)) {
    write.csv(pair_tbl, file.path(out_dir, paste0(spec$id, "_pairwise.csv")), row.names = FALSE)
  }

  lrt_tbl <- NULL
  if (!is.null(spec$lrt)) {
    lrt_raw <- read.csv(file.path(stats_dir, spec$lrt), stringsAsFactors = FALSE)
    lrt_tbl <- clean_lrt(lrt_raw)
    write.csv(lrt_tbl, file.path(out_dir, paste0(spec$id, "_lrt.csv")), row.names = FALSE)
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

  if (!is.null(lrt_tbl)) {
    doc <- body_add_par(doc, "")
    doc <- doc_add_table(doc, "Likelihood-ratio tests", lrt_tbl, table_type = "lrt")
  }

  if (!is.null(pair_tbl)) {
    doc <- body_add_par(doc, "")
    doc <- doc_add_table(doc, "Pairwise contrasts", pair_tbl, table_type = "pairwise")
  }

  print(doc, target = file.path(out_dir, paste0(spec$id, "_full_model_table.docx")))
}

message("Full model CSVs and DOCX tables generated in: ", out_dir)
