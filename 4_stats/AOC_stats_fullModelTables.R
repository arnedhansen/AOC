#!/usr/bin/env Rscript
# Builds Word / CSV full model tables from exported GLMM stats CSVs.
# Run after AOC_stats_glmm_export.R (or run export script, which sources this file).

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
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
source(file.path(script_dir, "AOC_stats_glmm_preprocess.R"))

stats_dir <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/data/stats"
out_dir <- file.path(stats_dir, "FullModelTables")
combined_dir <- file.path(out_dir, "combined_tables")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
}

# Usable text width for a standard Word page with 1-inch margins (letter).
SUMMARY_PAIRWISE_TABLE_WIDTH_IN <- 6.5

# Font size (pt) for the combined summary results tables.
SUMMARY_FONT_SIZE <- 7

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
  out[out == "AlphaPower_baselineWindow"] <- "Baseline-window alpha power"
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

format_load_level <- function(lvl, load_ref = NULL) {
  lvl <- as.character(lvl)
  if (lvl %in% c(".L", "Load.L")) return("linear trend")
  if (lvl %in% c(".Q", "Load.Q")) return("quadratic trend")
  ref_txt <- if (!is.null(load_ref) && nzchar(load_ref)) load_ref else "reference"
  paste0(lvl, " vs ", ref_txt)
}

clean_term <- function(term, load_ref = NULL) {
  term <- as.character(term)
  out <- term
  out[out == "Intercept"] <- "Intercept"
  out[out == "Group Var"] <- "Random intercept variance (Subject)"
  out[out == "Load"] <- "Load"
  out[out == "Interaction"] <- "Interaction term"
  out[out == "GazeDeviation_c"] <- "Gaze deviation"
  out[out == "MSRate_c"] <- "Microsaccade rate"
  out[out == "GazeDeviationBL_c"] <- "Gaze deviation baseline"
  out[out == "MSRateBL_c"] <- "Microsaccade rate baseline"

  # lme4 Load main effects (orthogonal polynomial or treatment coding)
  load_main <- grepl("^Load", out) & !grepl(":", out)
  if (any(load_main)) {
    lvls <- sub("^Load", "", out[load_main])
    out[load_main] <- vapply(
      lvls,
      function(lvl) {
        if (lvl %in% c(".L", "Load.L")) return("Load (linear trend)")
        if (lvl %in% c(".Q", "Load.Q")) return("Load (quadratic trend)")
        paste0("Load: ", format_load_level(lvl, load_ref))
      },
      character(1)
    )
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
        GazeDeviation = "Gaze deviation",
        MSRate = "Microsaccade rate",
        GazeDeviationBL = "Gaze deviation baseline",
        MSRateBL = "Microsaccade rate baseline",
        pred
      )
      load_txt <- if (lvl %in% c(".L", "Load.L")) {
        "linear trend"
      } else if (lvl %in% c(".Q", "Load.Q")) {
        "quadratic trend"
      } else {
        format_load_level(lvl, load_ref)
      }
      out[i] <- paste0("Interaction: ", pred, " x Load (", load_txt, ")")
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
    pred <- ifelse(pred == "GazeDeviation_c", "Gaze deviation", pred)
    pred <- ifelse(pred == "MSRate_c", "Microsaccade rate", pred)
    out[int_hits] <- paste0("Interaction: ", pred, " x Condition (", lvl, " vs ", ref, ")")
  }
  gsub(" \\(centered\\)", "", out, fixed = FALSE)
}

# Ordered-factor Load expands to polynomial contrasts (Load.L / Load.Q).
# Load effects are reported via emmeans pairwise contrasts, not these terms.
is_load_poly_term <- function(term) {
  grepl("Load\\.(L|Q)", as.character(term))
}

clean_fixed <- function(df, load_ref = NULL) {
  out <- df
  out <- out[!is_load_poly_term(out$Term), , drop = FALSE]
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

task_display_label <- function(task_name) {
  switch(
    tolower(as.character(task_name)),
    nback = "N-back",
    sternberg = "Sternberg",
    as.character(task_name)
  )
}

dv_unit <- function(dv_name) {
  switch(
    as.character(dv_name),
    Accuracy = "%",
    ReactionTime = "s",
    GazeDeviation = "px",
    MSRate = "Hz",
    GazeDeviationBL = "px",
    MSRateBL = "Hz",
    AlphaPower_baselineWindow = "\u03bcV\u00b2/Hz",
    ERSD = "dB",
    ""
  )
}

summary_table_title <- function(dv_name) {
  switch(
    as.character(dv_name),
    Accuracy = "Accuracy",
    ReactionTime = "Reaction time",
    GazeDeviation = "Gaze deviation",
    MSRate = "Microsaccade rate",
    GazeDeviationBL = "Gaze deviation (baseline-corrected)",
    MSRateBL = "Microsaccade rate (baseline-corrected)",
    AlphaPower_baselineWindow = "Baseline-window alpha power",
    ERSD = "ERS/ERD",
    display_dv(dv_name)
  )
}

desc_digits <- function(dv_name) {
  if (dv_name %in% c(
    "Accuracy", "ReactionTime", "GazeDeviation", "MSRate",
    "GazeDeviationBL", "MSRateBL"
  )) {
    2L
  } else {
    3L
  }
}

fmt_msd <- function(mean, sd, digits = 2L) {
  if (!is.finite(mean) || !is.finite(sd)) {
    return("")
  }
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean, sd)
}

PAIRWISE_CONTRAST_LABELS <- c(
  "Lowest vs middle load",
  "Lowest vs highest load",
  "Middle vs highest load"
)

wm_load_row_label <- function(level_idx) {
  paste0(
    NBACK_TASK_CONFIG$categories[level_idx],
    " / ",
    STERNBERG_TASK_CONFIG$categories[level_idx]
  )
}

task_descriptives <- function(dat, dv_name, categories) {
  vapply(
    categories,
    function(cat) {
      x <- dat[[dv_name]][dat$Load == cat]
      x <- x[!is.na(x)]
      if (length(x) == 0) {
        return("")
      }
      fmt_msd(mean(x), stats::sd(x), digits = desc_digits(dv_name))
    },
    character(1)
  )
}

pairwise_task_rows <- function(pair_df, eff_df, task_name, dv_name, comparisons) {
  if (is.null(pair_df) || nrow(pair_df) == 0) {
    return(NULL)
  }
  sub <- pair_df[pair_df$Variable == dv_name, , drop = FALSE]
  if (nrow(sub) == 0) {
    return(NULL)
  }
  sub <- merge_pairwise_effectsizes(sub, eff_df = eff_df, task_name = task_name, dv_name = dv_name)

  ci_low_col <- if ("CI95_low" %in% names(sub)) "CI95_low" else "CI_low"
  ci_high_col <- if ("CI95_high" %in% names(sub)) "CI95_high" else "CI_high"
  p_col <- if ("p_adj" %in% names(sub)) "p_adj" else "p"
  z_col <- if ("z" %in% names(sub)) "z" else "stat"

  rows <- vector("list", length(comparisons))
  for (i in seq_along(comparisons)) {
    low_lvl <- comparisons[[i]][1]
    high_lvl <- comparisons[[i]][2]
    row <- sub[sub$Group1 == high_lvl & sub$Group2 == low_lvl, , drop = FALSE]
    if (nrow(row) == 0) {
      rows[[i]] <- list(
        beta = NA_real_,
        se = NA_real_,
        statistic = NA_real_,
        padj = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        cohens_d = NA_real_
      )
    } else {
      dz <- if ("Cohens_dz" %in% names(row)) as.numeric(row$Cohens_dz[1]) else NA_real_
      rows[[i]] <- list(
        beta = -as.numeric(row$Estimate[1]),
        se = as.numeric(row$SE[1]),
        statistic = -as.numeric(row[[z_col]][1]),
        padj = as.numeric(row[[p_col]][1]),
        ci_low = -as.numeric(row[[ci_high_col]][1]),
        ci_high = -as.numeric(row[[ci_low_col]][1]),
        cohens_d = if (is.finite(dz)) dz else NA_real_
      )
    }
  }
  rows
}

pairwise_task_matrix <- function(rows, n_contrasts) {
  cols <- c("beta", "ci", "statistic", "padj", "cohens_d")
  if (is.null(rows)) {
    empty <- as.data.frame(
      matrix("", nrow = n_contrasts, ncol = length(cols)),
      stringsAsFactors = FALSE
    )
    names(empty) <- cols
    return(empty)
  }
  mat <- do.call(
    rbind,
    lapply(rows, function(row) {
      c(
        fmt_num(row$beta),
        fmt_ci(row$ci_low, row$ci_high),
        fmt_num(row$statistic),
        fmt_p(row$padj),
        fmt_num(row$cohens_d)
      )
    })
  )
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- cols
  df
}

build_descriptive_table <- function(nback_dat, stern_dat, dv_name) {
  unit <- dv_unit(dv_name)
  wm_labels <- vapply(seq_along(NBACK_TASK_CONFIG$categories), wm_load_row_label, character(1))
  nback_vals <- task_descriptives(nback_dat, dv_name, NBACK_TASK_CONFIG$categories)
  stern_vals <- task_descriptives(stern_dat, dv_name, STERNBERG_TASK_CONFIG$categories)

  list(
    body = data.frame(
      wm_load = wm_labels,
      nback = nback_vals,
      sternberg = stern_vals,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    unit = unit
  )
}

PAIRWISE_STAT_LABEL <- "t-value"

build_pairwise_table <- function(pair_nback, pair_stern, eff_nback, eff_stern, dv_name) {
  n_contrasts <- length(NBACK_TASK_CONFIG$comparisons)
  nback_rows <- pairwise_task_rows(
    pair_nback, eff_nback, "nback", dv_name, NBACK_TASK_CONFIG$comparisons
  )
  stern_rows <- pairwise_task_rows(
    pair_stern, eff_stern, "sternberg", dv_name, STERNBERG_TASK_CONFIG$comparisons
  )
  if (is.null(nback_rows) && is.null(stern_rows)) {
    return(NULL)
  }

  nback_mat <- pairwise_task_matrix(nback_rows, n_contrasts)
  stern_mat <- pairwise_task_matrix(stern_rows, n_contrasts)
  body <- cbind(
    data.frame(contrast = PAIRWISE_CONTRAST_LABELS[seq_len(n_contrasts)], stringsAsFactors = FALSE),
    nback_mat,
    stern_mat
  )
  names(body) <- c(
    "contrast",
    paste0("nback_", names(nback_mat)),
    paste0("sternberg_", names(stern_mat))
  )
  list(body = body)
}

# Width of the empty spacer column that separates the N-back and Sternberg
# blocks, creating the gap between the two task underlines.
SUMMARY_SPACER_WIDTH_IN <- 0.18

# Rule (line) thickness and row spacing for the combined summary tables.
# Half of the theme_booktabs default (1.5pt) for thinner lines.
SUMMARY_RULE_WIDTH <- 0.75
SUMMARY_LINE_SPACING <- 0.75
SUMMARY_ROW_HEIGHT_IN <- 0.14

add_task_header_rules <- function(ft, title_row = 1, j_nback, j_stern) {
  # Match the color of theme_booktabs rules; thickness is 0.85x the other rules.
  bdr <- fp_border(color = "#666666", width = SUMMARY_RULE_WIDTH * 0.85)
  # Draw the underline directly on the task-title row so no empty row sits
  # between the titles and their rules.
  ft <- valign(ft, i = title_row, valign = "bottom", part = "header")
  ft <- hline(ft, i = title_row, j = j_nback, border = bdr, part = "header")
  ft <- hline(ft, i = title_row, j = j_stern, border = bdr, part = "header")
  # Header font size is applied last: fontsize(part = "header") fails earlier
  # in the build because of the merged/spanned task-title row.
  fontsize(ft, size = SUMMARY_FONT_SIZE, part = "header")
}

make_descriptive_flextable <- function(desc_tbl) {
  unit <- desc_tbl$unit
  msd_lab <- paste0("M (SD) [", unit, "]")
  src <- desc_tbl$body
  # Spacer column between the two task columns keeps their underlines apart.
  body <- data.frame(
    wm_load = src[[1]],
    nback = src[[2]],
    .spacer = "",
    sternberg = src[[3]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ft <- flextable(body)
  ft <- delete_part(ft, part = "header")
  ft <- add_header_row(ft, values = c("WM load", msd_lab, "", msd_lab))
  ft <- add_header_row(ft, values = c("", "N-back", "", "Sternberg"), top = TRUE)
  ft <- style_summary_ft(
    ft,
    col_widths = c(2.4, 1.35, SUMMARY_SPACER_WIDTH_IN, 1.35)
  )
  add_task_header_rules(ft, title_row = 1, j_nback = 2, j_stern = 4)
}

make_pairwise_flextable <- function(pair_tbl) {
  stat_header <- c("\u03b2", "CI", PAIRWISE_STAT_LABEL, "padj", "Cohen's d")
  n_task_cols <- length(stat_header)
  src <- pair_tbl$body
  nback_block <- src[, 2:(1 + n_task_cols), drop = FALSE]
  stern_block <- src[, (2 + n_task_cols):(1 + 2 * n_task_cols), drop = FALSE]
  # Spacer column between the two task blocks keeps their underlines apart.
  body <- cbind(
    data.frame(contrast = src[[1]], stringsAsFactors = FALSE),
    nback_block,
    data.frame(.spacer = rep("", nrow(src)), stringsAsFactors = FALSE),
    stern_block
  )
  j_nback <- 2:(1 + n_task_cols)
  j_stern <- (3 + n_task_cols):(2 + 2 * n_task_cols)
  task_row <- c(
    "",
    "N-back",
    rep("", n_task_cols - 1L),
    "",
    "Sternberg",
    rep("", n_task_cols - 1L)
  )
  ft <- flextable(body)
  ft <- delete_part(ft, part = "header")
  ft <- add_header_row(ft, values = c("Contrast", stat_header, "", stat_header))
  ft <- add_header_row(ft, values = task_row, top = TRUE)
  ft <- merge_at(ft, i = 1, j = j_nback, part = "header")
  ft <- merge_at(ft, i = 1, j = j_stern, part = "header")
  # Relative widths (beta, CI, t-value, padj, Cohen's d); scaled to the page
  # width in style_summary_ft(). CI holds the confidence interval and the last
  # column is wide enough to keep the "Cohen's d" header on one line.
  per_task_widths <- c(0.40, 0.95, 0.45, 0.33, 0.55)
  ft <- style_summary_ft(
    ft,
    col_widths = c(1.15, per_task_widths, SUMMARY_SPACER_WIDTH_IN, per_task_widths),
    full_width = TRUE
  )
  add_task_header_rules(ft, title_row = 1, j_nback = j_nback, j_stern = j_stern)
}

doc_add_section_table <- function(doc, section_title, ft) {
  doc <- doc_add_title(doc, section_title, size = SUMMARY_FONT_SIZE)
  body_add_flextable(doc, ft)
}

build_measure_results_doc <- function(
  spec,
  nback_dat,
  stern_dat,
  pair_nback,
  pair_stern,
  eff_nback,
  eff_stern,
  out_path
) {
  desc_tbl <- build_descriptive_table(nback_dat, stern_dat, spec$dv)
  pair_tbl <- build_pairwise_table(pair_nback, pair_stern, eff_nback, eff_stern, spec$dv)
  if (is.null(pair_tbl)) {
    message("Skipping summary table ", spec$id, " (no pairwise contrasts available).")
    return(invisible(NULL))
  }

  csv_desc_path <- file.path(combined_dir, paste0(spec$id, "_descriptive.csv"))
  csv_pair_path <- file.path(combined_dir, paste0(spec$id, "_pairwise.csv"))
  write.csv(desc_tbl$body, csv_desc_path, row.names = FALSE)
  write.csv(pair_tbl$body, csv_pair_path, row.names = FALSE)

  doc <- read_docx()
  doc <- doc_add_title(doc, summary_table_title(spec$dv), size = SUMMARY_FONT_SIZE)
  doc <- body_add_par(doc, "")
  doc <- doc_add_section_table(doc, "Descriptive Statistics", make_descriptive_flextable(desc_tbl))
  doc <- body_add_par(doc, "")
  doc <- doc_add_section_table(doc, "Pairwise Contrasts", make_pairwise_flextable(pair_tbl))
  doc <- body_add_par(doc, "")
  doc <- body_add_fpar(
    doc,
    fpar(
      ftext(
        paste(
          "Note. Descriptive statistics are means with standard deviations in parentheses.",
          "Pairwise contrasts are from linear mixed-effects models with random intercepts for Subject.",
          "\u03b2 reflects the higher-load minus lower-load contrast within each task.",
          "padj = FDR-adjusted p-value.",
          "Full model estimates are reported in the supplementary tables."
        ),
        fp_text(font.size = SUMMARY_FONT_SIZE)
      )
    )
  )

  print(doc, target = out_path)
  message("Wrote summary table -> ", out_path)
  invisible(out_path)
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

style_summary_ft <- function(ft, col_widths, full_width = FALSE) {
  ft <- theme_booktabs(ft)
  # Halve the thickness of the booktabs rules (top, header, bottom).
  rule_border <- fp_border(color = "#666666", width = SUMMARY_RULE_WIDTH)
  ft <- hline_top(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "body")
  ft <- fontsize(ft, size = SUMMARY_FONT_SIZE, part = "body")
  ft <- bold(ft, part = "header")
  ft <- padding(ft, padding = 1, part = "all")
  ft <- align(ft, align = "left", j = 1, part = "all")
  if (ncol(ft$body$dataset) > 1) {
    ft <- align(ft, align = "center", j = 2:ncol(ft$body$dataset), part = "all")
  }
  ft <- line_spacing(ft, space = SUMMARY_LINE_SPACING, part = "all")
  ft <- hrule(ft, rule = "exact", part = "all")
  ft <- height_all(ft, height = SUMMARY_ROW_HEIGHT_IN, part = "all")
  if (full_width && !is.null(col_widths)) {
    # Scale columns so they exactly fill the page width, keeping the font at
    # 8pt. fit_to_width() is avoided here because it shrinks the font instead.
    col_widths <- col_widths * (SUMMARY_PAIRWISE_TABLE_WIDTH_IN / sum(col_widths))
  }
  if (!is.null(col_widths) && length(col_widths) == ncol(ft$body$dataset)) {
    ft <- width(ft, j = seq_along(col_widths), width = col_widths)
    ft <- set_table_properties(ft, layout = "fixed")
  } else {
    ft <- set_table_properties(ft, layout = "fixed")
    ft <- autofit(ft)
  }
  ft <- align(ft, align = "center", part = "header")
  align(ft, j = 1, align = "left", part = "header")
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

append_model_section <- function(
  doc,
  section_title,
  model_text = NULL,
  fixed_tbl,
  random_tbl = NULL,
  lrt_tbl = NULL,
  pair_tbl = NULL
) {
  doc <- doc_add_title(doc, section_title, size = 11)
  if (!is.null(model_text) && nzchar(model_text)) {
    doc <- body_add_fpar(
      doc,
      fpar(ftext(paste0("Model: ", model_text), fp_text(font.size = 9)))
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
  doc
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
    id = "nback_pow_baselineWindow",
    title = "N-back Baseline-Window Alpha Power (exploratory)",
    fixed = "AOC_mixedlm_fixed_pow_baselineWindow_nback.csv",
    pair_src = if (nrow(pair_nback) > 0) pair_nback[pair_nback$Variable == "AlphaPower_baselineWindow", ] else NULL
  ),
  list(
    id = "sternberg_pow_baselineWindow",
    title = "Sternberg Baseline-Window Alpha Power (exploratory)",
    fixed = "AOC_mixedlm_fixed_pow_baselineWindow_sternberg.csv",
    pair_src = if (nrow(pair_stern) > 0) pair_stern[pair_stern$Variable == "AlphaPower_baselineWindow", ] else NULL
  ),
  list(
    id = "nback_ersd_by_gazedeviation",
    title = "N-back ERS/ERD by GazeDeviation",
    fixed = "AOC_ersd_fixed_gazedeviation_nback.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_nback.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_nback.csv"
  ),
  list(
    id = "nback_ersd_by_gazedeviation_full",
    title = "N-back ERS/ERD by GazeDeviation [full interaction]",
    fixed = "AOC_ersd_fixed_gazedeviation_nback_full.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_nback.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_nback_full.csv"
  ),
  list(
    id = "nback_ersd_by_gazedeviation_reduced",
    title = "N-back ERS/ERD by GazeDeviation [additive]",
    fixed = "AOC_ersd_fixed_gazedeviation_nback_reduced.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_nback.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_nback_reduced.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedeviation",
    title = "Sternberg ERS/ERD by GazeDeviation",
    fixed = "AOC_ersd_fixed_gazedeviation_sternberg.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_sternberg.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_sternberg.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedeviation_full",
    title = "Sternberg ERS/ERD by GazeDeviation [full interaction]",
    fixed = "AOC_ersd_fixed_gazedeviation_sternberg_full.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_sternberg.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_sternberg_full.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedeviation_reduced",
    title = "Sternberg ERS/ERD by GazeDeviation [additive]",
    fixed = "AOC_ersd_fixed_gazedeviation_sternberg_reduced.csv",
    lrt = "AOC_ersd_lrtTbl_gazedeviation_sternberg.csv",
    pair = "AOC_ersd_contrasts_gazedeviation_sternberg_reduced.csv"
  ),
  list(
    id = "nback_ersd_by_msrate",
    title = "N-back ERS/ERD by MSRate",
    fixed = "AOC_ersd_fixed_msrate_nback.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_nback.csv",
    pair = "AOC_ersd_contrasts_msrate_nback.csv"
  ),
  list(
    id = "nback_ersd_by_msrate_full",
    title = "N-back ERS/ERD by MSRate [full interaction]",
    fixed = "AOC_ersd_fixed_msrate_nback_full.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_nback.csv",
    pair = "AOC_ersd_contrasts_msrate_nback_full.csv"
  ),
  list(
    id = "nback_ersd_by_msrate_reduced",
    title = "N-back ERS/ERD by MSRate [additive]",
    fixed = "AOC_ersd_fixed_msrate_nback_reduced.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_nback.csv",
    pair = "AOC_ersd_contrasts_msrate_nback_reduced.csv"
  ),
  list(
    id = "sternberg_ersd_by_msrate",
    title = "Sternberg ERS/ERD by MSRate",
    fixed = "AOC_ersd_fixed_msrate_sternberg.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_sternberg.csv",
    pair = "AOC_ersd_contrasts_msrate_sternberg.csv"
  ),
  list(
    id = "sternberg_ersd_by_msrate_full",
    title = "Sternberg ERS/ERD by MSRate [full interaction]",
    fixed = "AOC_ersd_fixed_msrate_sternberg_full.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_sternberg.csv",
    pair = "AOC_ersd_contrasts_msrate_sternberg_full.csv"
  ),
  list(
    id = "sternberg_ersd_by_msrate_reduced",
    title = "Sternberg ERS/ERD by MSRate [additive]",
    fixed = "AOC_ersd_fixed_msrate_sternberg_reduced.csv",
    lrt = "AOC_ersd_lrtTbl_msrate_sternberg.csv",
    pair = "AOC_ersd_contrasts_msrate_sternberg_reduced.csv"
  ),
  list(
    id = "nback_ersd_by_gazedev_bl",
    title = "N-back ERS/ERD by GazeDeviationBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_nback.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_nback.csv"
  ),
  list(
    id = "nback_ersd_by_gazedev_bl_full",
    title = "N-back ERS/ERD by GazeDeviationBL [full interaction] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_nback_full.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_nback_full.csv"
  ),
  list(
    id = "nback_ersd_by_gazedev_bl_reduced",
    title = "N-back ERS/ERD by GazeDeviationBL [additive] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_nback_reduced.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_nback_reduced.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedev_bl",
    title = "Sternberg ERS/ERD by GazeDeviationBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_sternberg.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_sternberg.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedev_bl_full",
    title = "Sternberg ERS/ERD by GazeDeviationBL [full interaction] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_sternberg_full.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_sternberg_full.csv"
  ),
  list(
    id = "sternberg_ersd_by_gazedev_bl_reduced",
    title = "Sternberg ERS/ERD by GazeDeviationBL [additive] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_gazedev_bl_sternberg_reduced.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_gazedev_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_gazedev_bl_sternberg_reduced.csv"
  ),
  list(
    id = "nback_ersd_by_ms_bl",
    title = "N-back ERS/ERD by MSRateBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_nback.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_nback.csv"
  ),
  list(
    id = "nback_ersd_by_ms_bl_full",
    title = "N-back ERS/ERD by MSRateBL [full interaction] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_nback_full.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_nback_full.csv"
  ),
  list(
    id = "nback_ersd_by_ms_bl_reduced",
    title = "N-back ERS/ERD by MSRateBL [additive] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_nback_reduced.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_nback.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_nback_reduced.csv"
  ),
  list(
    id = "sternberg_ersd_by_ms_bl",
    title = "Sternberg ERS/ERD by MSRateBL (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_sternberg.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_sternberg.csv"
  ),
  list(
    id = "sternberg_ersd_by_ms_bl_full",
    title = "Sternberg ERS/ERD by MSRateBL [full interaction] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_sternberg_full.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_sternberg_full.csv"
  ),
  list(
    id = "sternberg_ersd_by_ms_bl_reduced",
    title = "Sternberg ERS/ERD by MSRateBL [additive] (exploratory)",
    fixed = "AOC_ersd_exploratory_fixed_ms_bl_sternberg_reduced.csv",
    lrt = "AOC_ersd_exploratory_lrtTbl_ms_bl_sternberg.csv",
    pair = "AOC_ersd_exploratory_contrasts_ms_bl_sternberg_reduced.csv"
  ),
  list(
    id = "nback_alpha_by_gazedeviation",
    title = "N-back AlphaPower by GazeDeviation (legacy)",
    fixed = "AOC_alpha_fixed_alphapower_gazedeviation_nback.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_gazedeviation_nback.csv",
    pair = "AOC_alpha_contrasts_alphapower_gazedeviation_nback.csv"
  ),
  list(
    id = "sternberg_alpha_by_gazedeviation",
    title = "Sternberg AlphaPower by GazeDeviation (legacy)",
    fixed = "AOC_alpha_fixed_alphapower_gazedeviation_sternberg.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_gazedeviation_sternberg.csv",
    pair = "AOC_alpha_contrasts_alphapower_gazedeviation_sternberg.csv"
  ),
  list(
    id = "nback_alpha_by_msrate",
    title = "N-back AlphaPower by MSRate (legacy)",
    fixed = "AOC_alpha_fixed_alphapower_msrate_nback.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_msrate_nback.csv",
    pair = "AOC_alpha_contrasts_alphapower_msrate_nback.csv"
  ),
  list(
    id = "sternberg_alpha_by_msrate",
    title = "Sternberg AlphaPower by MSRate (legacy)",
    fixed = "AOC_alpha_fixed_alphapower_msrate_sternberg.csv",
    lrt = "AOC_alpha_lrtTbl_alphapower_msrate_sternberg.csv",
    pair = "AOC_alpha_contrasts_alphapower_msrate_sternberg.csv"
  )
)

built <- list()

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

  formula_text <- unique(fixed_raw$ModelLabel)
  dv_text <- unique(fixed_raw$DV)
  dv_text <- if (length(dv_text) > 0 && !is.na(dv_text[1])) display_dv(dv_text[1]) else "Outcome"
  model_text <- NULL
  if (length(formula_text) > 0 && !is.na(formula_text[1])) {
    model_text <- prettify_model_label(formula_text[1], dv_label = dv_text)
  }

  built[[spec$id]] <- list(
    title = spec$title,
    model_text = model_text,
    fixed_tbl = fixed_tbl,
    random_tbl = random_tbl,
    lrt_tbl = lrt_tbl,
    pair_tbl = pair_tbl
  )

  doc <- read_docx()
  doc <- append_model_section(
    doc = doc,
    section_title = spec$title,
    model_text = model_text,
    fixed_tbl = fixed_tbl,
    random_tbl = random_tbl,
    lrt_tbl = lrt_tbl,
    pair_tbl = pair_tbl
  )

  print(doc, target = file.path(out_dir, paste0(spec$id, "_full_model_table.docx")))
}

# ---- Co-variation (ERS/ERD coupling) summary tables ------------------------
# Builds combined N-back | Sternberg summary documents for the H5 models that
# test whether ERS/ERD coupling with an oculomotor predictor differs across WM
# load. One document per predictor (Gaze Deviation, Microsaccade Rate).

read_covariation_csv <- function(kind, predictor_key, task, exploratory = FALSE, variant = NULL) {
  prefix <- if (isTRUE(exploratory)) "AOC_ersd_exploratory" else "AOC_ersd"
  if (is.null(variant) || !nzchar(variant)) {
    fname <- sprintf("%s_%s_%s_%s.csv", prefix, kind, predictor_key, task)
  } else if (kind == "lrtTbl") {
    # LRT is shared across full / reduced / selected variants.
    fname <- sprintf("%s_%s_%s_%s.csv", prefix, kind, predictor_key, task)
  } else {
    fname <- sprintf("%s_%s_%s_%s_%s.csv", prefix, kind, predictor_key, task, variant)
  }
  safe_read_csv(file.path(stats_dir, fname))
}

build_covariation_fixed_matrix <- function(fixed_df, terms) {
  n <- length(terms)
  beta <- character(n)
  ci <- character(n)
  p <- character(n)
  for (i in seq_along(terms)) {
    row <- if (!is.null(fixed_df)) {
      fixed_df[fixed_df$Term == terms[[i]]$term, , drop = FALSE]
    } else {
      NULL
    }
    if (!is.null(row) && nrow(row) > 0) {
      beta[i] <- fmt_num(row$beta[1])
      ci[i] <- fmt_ci(row$CI_low[1], row$CI_high[1])
      p[i] <- fmt_p(row$p[1])
    }
  }
  data.frame(beta = beta, ci = ci, p = p, stringsAsFactors = FALSE)
}

build_covariation_lrt_matrix <- function(lrt_df) {
  if (is.null(lrt_df) || nrow(lrt_df) == 0) {
    return(data.frame(chi2 = "", df = "", p = "", stringsAsFactors = FALSE))
  }
  row <- lrt_df[lrt_df$Term == "Interaction", , drop = FALSE]
  if (nrow(row) == 0) row <- lrt_df[1, , drop = FALSE]
  data.frame(
    chi2 = fmt_num(row$LR_Chi2[1]),
    df = as.character(row$df[1]),
    p = fmt_p(row$p[1]),
    stringsAsFactors = FALSE
  )
}

build_covariation_pairwise_matrix <- function(contr_df, comparisons) {
  n <- length(comparisons)
  beta <- character(n)
  ci <- character(n)
  p <- character(n)
  if (!is.null(contr_df) && nrow(contr_df) > 0) {
    for (i in seq_along(comparisons)) {
      low <- comparisons[[i]][1]
      high <- comparisons[[i]][2]
      # Contrast rows are stored as Group1 = higher load, Group2 = lower load,
      # so Estimate already reflects the higher-load minus lower-load direction.
      row <- contr_df[contr_df$Group1 == high & contr_df$Group2 == low, , drop = FALSE]
      if (nrow(row) > 0) {
        beta[i] <- fmt_num(row$Estimate[1])
        ci[i] <- fmt_ci(row$CI95_low[1], row$CI95_high[1])
        p[i] <- fmt_p(row$p_adj[1])
      }
    }
  }
  data.frame(beta = beta, ci = ci, p = p, stringsAsFactors = FALSE)
}

make_grouped_two_task_ft <- function(first_header, stat_header, row_labels,
                                     nback_mat, stern_mat,
                                     first_col_width = 1.8,
                                     per_task_widths = c(0.6, 1.0, 0.5)) {
  n_task_cols <- length(stat_header)
  names(nback_mat) <- paste0("nb_", seq_len(ncol(nback_mat)))
  names(stern_mat) <- paste0("st_", seq_len(ncol(stern_mat)))
  body <- cbind(
    data.frame(rowlab = row_labels, stringsAsFactors = FALSE),
    nback_mat,
    data.frame(.spacer = rep("", length(row_labels)), stringsAsFactors = FALSE),
    stern_mat
  )
  j_nback <- 2:(1 + n_task_cols)
  j_stern <- (3 + n_task_cols):(2 + 2 * n_task_cols)
  task_row <- c(
    "",
    "N-back",
    rep("", n_task_cols - 1L),
    "",
    "Sternberg",
    rep("", n_task_cols - 1L)
  )
  ft <- flextable(body)
  ft <- delete_part(ft, part = "header")
  ft <- add_header_row(ft, values = c(first_header, stat_header, "", stat_header))
  ft <- add_header_row(ft, values = task_row, top = TRUE)
  ft <- merge_at(ft, i = 1, j = j_nback, part = "header")
  ft <- merge_at(ft, i = 1, j = j_stern, part = "header")
  ft <- style_summary_ft(
    ft,
    col_widths = c(first_col_width, per_task_widths, SUMMARY_SPACER_WIDTH_IN, per_task_widths),
    full_width = TRUE
  )
  add_task_header_rules(ft, title_row = 1, j_nback = j_nback, j_stern = j_stern)
}

build_covariation_results_doc <- function(spec, out_path, variant = NULL) {
  exploratory <- isTRUE(spec$exploratory)
  # Gaze / MS main effect only. Load is reported via emmeans pairwise contrasts below.
  terms <- list(
    list(term = spec$fixed_term, label = spec$label)
  )
  term_labels <- vapply(terms, function(t) t$label, character(1))

  fixed_nb <- read_covariation_csv("fixed", spec$key, "nback", exploratory = exploratory, variant = variant)
  fixed_st <- read_covariation_csv("fixed", spec$key, "sternberg", exploratory = exploratory, variant = variant)
  lrt_nb <- read_covariation_csv("lrtTbl", spec$key, "nback", exploratory = exploratory)
  lrt_st <- read_covariation_csv("lrtTbl", spec$key, "sternberg", exploratory = exploratory)
  contr_nb <- read_covariation_csv("contrasts", spec$key, "nback", exploratory = exploratory, variant = variant)
  contr_st <- read_covariation_csv("contrasts", spec$key, "sternberg", exploratory = exploratory, variant = variant)

  fixed_nb_mat <- build_covariation_fixed_matrix(fixed_nb, terms)
  fixed_st_mat <- build_covariation_fixed_matrix(fixed_st, terms)
  lrt_nb_mat <- build_covariation_lrt_matrix(lrt_nb)
  lrt_st_mat <- build_covariation_lrt_matrix(lrt_st)
  pw_nb_mat <- build_covariation_pairwise_matrix(contr_nb, NBACK_TASK_CONFIG$comparisons)
  pw_st_mat <- build_covariation_pairwise_matrix(contr_st, STERNBERG_TASK_CONFIG$comparisons)

  id_suffix <- if (is.null(variant) || !nzchar(variant)) {
    ""
  } else {
    paste0("_", variant)
  }

  write.csv(
    cbind(term = term_labels, nback = fixed_nb_mat, sternberg = fixed_st_mat),
    file.path(combined_dir, paste0(spec$id, id_suffix, "_fixed.csv")),
    row.names = FALSE
  )
  write.csv(
    cbind(
      effect = paste0("Load x ", spec$label),
      nback = lrt_nb_mat,
      sternberg = lrt_st_mat
    ),
    file.path(combined_dir, paste0(spec$id, id_suffix, "_interaction_lrt.csv")),
    row.names = FALSE
  )
  write.csv(
    cbind(
      contrast = PAIRWISE_CONTRAST_LABELS,
      nback = pw_nb_mat,
      sternberg = pw_st_mat
    ),
    file.path(combined_dir, paste0(spec$id, id_suffix, "_pairwise.csv")),
    row.names = FALSE
  )

  beta_sym <- "\u03b2"
  chi_sym <- "\u03c7\u00b2"
  interaction_label <- paste0("Load \u00d7 ", spec$label)

  model_predictor <- sub("_c$", "", spec$fixed_term)
  model_labels <- c()
  if (!is.null(fixed_nb) && nrow(fixed_nb) > 0 && "ModelLabel" %in% names(fixed_nb)) {
    model_labels <- c(model_labels, as.character(fixed_nb$ModelLabel[1]))
  }
  if (!is.null(fixed_st) && nrow(fixed_st) > 0 && "ModelLabel" %in% names(fixed_st)) {
    model_labels <- c(model_labels, as.character(fixed_st$ModelLabel[1]))
  }
  model_labels <- unique(model_labels[!is.na(model_labels) & nzchar(model_labels)])
  model_text <- if (length(model_labels) == 1L) {
    prettify_model_label(model_labels[1], dv_label = "ERS/ERD")
  } else if (identical(variant, "reduced")) {
    paste0("ERS/ERD ~ ", model_predictor, " + Load + (1|Subject)")
  } else {
    paste0("ERS/ERD ~ ", model_predictor, " * Load + (1|Subject)")
  }
  full_model_text <- paste0("ERS/ERD ~ ", model_predictor, " * Load + (1|Subject)")
  reduced_model_text <- paste0("ERS/ERD ~ ", model_predictor, " + Load + (1|Subject)")

  variant_tag <- if (identical(variant, "full")) {
    " [full interaction model]"
  } else if (identical(variant, "reduced")) {
    " [additive model]"
  } else {
    ""
  }
  doc_title <- paste0("ERS/ERD Co-Variation with ", spec$label, variant_tag)
  if (exploratory) {
    doc_title <- paste0(doc_title, " (exploratory)")
  }

  note_body <- if (identical(variant, "full")) {
    paste0(
      "Fixed effects and pairwise load contrasts are reported from the full interaction model (",
      full_model_text, "). "
    )
  } else if (identical(variant, "reduced")) {
    paste0(
      "Fixed effects and pairwise load contrasts are reported from the additive model (",
      reduced_model_text, "). "
    )
  } else {
    paste0(
      "When the interaction was not significant (p \u2265 .05), the interaction term was removed; ",
      "fixed effects and pairwise load contrasts are then reported from the selected model (",
      model_text, "). ",
      "Separate tables also report the full interaction model and the additive model. "
    )
  }

  doc <- read_docx()
  doc <- doc_add_title(
    doc,
    doc_title,
    size = SUMMARY_FONT_SIZE
  )
  doc <- body_add_fpar(
    doc,
    fpar(ftext(paste0("Model: ", model_text), fp_text(font.size = SUMMARY_FONT_SIZE)))
  )
  doc <- body_add_par(doc, "")
  doc <- doc_add_section_table(
    doc,
    "Interaction (Likelihood-Ratio Test)",
    make_grouped_two_task_ft(
      "Effect",
      c(chi_sym, "df", "p"),
      interaction_label,
      lrt_nb_mat,
      lrt_st_mat,
      per_task_widths = c(0.7, 0.7, 0.7)
    )
  )
  doc <- body_add_par(doc, "")
  doc <- doc_add_section_table(
    doc,
    "Fixed Effects",
    make_grouped_two_task_ft(
      "Term",
      c(beta_sym, "CI", "p"),
      term_labels,
      fixed_nb_mat,
      fixed_st_mat
    )
  )
  doc <- body_add_par(doc, "")
  doc <- doc_add_section_table(
    doc,
    "Pairwise Contrasts",
    make_grouped_two_task_ft(
      "Contrast",
      c(beta_sym, "CI", "padj"),
      PAIRWISE_CONTRAST_LABELS,
      pw_nb_mat,
      pw_st_mat
    )
  )
  doc <- body_add_par(doc, "")
  exploratory_note <- if (exploratory) {
    "This analysis is exploratory. "
  } else {
    ""
  }
  doc <- body_add_fpar(
    doc,
    fpar(
      ftext(
        paste0(
          "Note. ", exploratory_note,
          "ERS/ERD coupling with ", tolower(spec$label),
          " was tested with linear mixed-effects models (random intercepts for Subject). ",
          "The preregistered full model was ", full_model_text, ". ",
          "A likelihood-ratio test compared this model to the additive model (",
          model_predictor, " + Load). ",
          note_body,
          "\u03b2 in the pairwise section reflects the higher-load minus lower-load contrast. ",
          "padj = FDR-adjusted p-value. ",
          "Full model estimates are reported in the supplementary tables."
        ),
        fp_text(font.size = SUMMARY_FONT_SIZE)
      )
    )
  )

  print(doc, target = out_path)
  message("Wrote co-variation table -> ", out_path)
  invisible(out_path)
}

covariation_specs <- list(
  list(
    id = "ersd_gaze_deviation",
    key = "gazedeviation",
    fixed_term = "GazeDeviation_c",
    label = "Gaze Deviation"
  ),
  list(
    id = "ersd_microsaccade_rate",
    key = "msrate",
    fixed_term = "MSRate_c",
    label = "Microsaccade Rate"
  ),
  list(
    id = "ersd_gaze_deviation_bl",
    key = "gazedev_bl",
    fixed_term = "GazeDeviationBL_c",
    label = "Gaze Deviation (baseline-corrected)",
    exploratory = TRUE
  ),
  list(
    id = "ersd_microsaccade_rate_bl",
    key = "ms_bl",
    fixed_term = "MSRateBL_c",
    label = "Microsaccade Rate (baseline-corrected)",
    exploratory = TRUE
  )
)

combined_specs <- list(
  list(id = "reaction_time", dv = "ReactionTime"),
  list(id = "accuracy", dv = "Accuracy"),
  list(id = "gaze_deviation", dv = "GazeDeviation"),
  list(id = "microsaccade_rate", dv = "MSRate"),
  list(id = "gaze_deviation_bl", dv = "GazeDeviationBL"),
  list(id = "microsaccade_rate_bl", dv = "MSRateBL"),
  list(id = "ers_erd", dv = "ERSD"),
  list(id = "alpha_power_baselineWindow", dv = "AlphaPower_baselineWindow")
)

nback_dat <- load_task_data(NBACK_TASK_CONFIG)
stern_dat <- load_task_data(STERNBERG_TASK_CONFIG)

for (spec in combined_specs) {
  if (!spec$dv %in% names(nback_dat) || !spec$dv %in% names(stern_dat)) {
    message("Skipping summary table ", spec$id, " (variable missing from merged data).")
    next
  }
  build_measure_results_doc(
    spec = spec,
    nback_dat = nback_dat,
    stern_dat = stern_dat,
    pair_nback = pair_nback,
    pair_stern = pair_stern,
    eff_nback = eff_nback,
    eff_stern = eff_stern,
    out_path = file.path(combined_dir, paste0(spec$id, "_results_table.docx"))
  )
}

for (spec in covariation_specs) {
  for (variant in list(NULL, "full", "reduced")) {
    id_suffix <- if (is.null(variant)) "" else paste0("_", variant)
    build_covariation_results_doc(
      spec = spec,
      out_path = file.path(combined_dir, paste0(spec$id, id_suffix, "_results_table.docx")),
      variant = variant
    )
  }
}

message("Full model CSVs and DOCX tables generated in: ", out_dir)
message("Summary results tables generated in: ", combined_dir)
