# Shared helpers for AOC GLMM fitting (condition- and trial-level).

# Grand-mean z-score a continuous predictor: (x - mean) / SD over the analysis sample.
# Slopes are then in DV units per 1 SD increase in the predictor.
grand_mean_zscore <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sigma <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) {
    warning("grand_mean_zscore: non-positive SD; returning mean-centered values only.")
    return(list(z = x - mu, mean = mu, sd = NA_real_))
  }
  list(z = (x - mu) / sigma, mean = mu, sd = sigma)
}

gaze_z_colname <- function(gaze_var) {
  paste0(gaze_var, "_z")
}

add_gaze_z_column <- function(dat, gaze_var) {
  gs <- grand_mean_zscore(dat[[gaze_var]])
  dat[[gaze_z_colname(gaze_var)]] <- gs$z
  dat
}

STANDARDIZED_GAZE_BETA_NOTE <- paste0(
  "Gaze predictor coefficients are grand-mean standardized. ",
  "Pairwise load contrasts on ERS/ERD are in dB. "
)

# Variance components appended to exported fixed-effects CSVs.
RE_TERM_SUBJECT <- "Group Var"
RE_TERM_TRIAL <- "Trial Var"
RE_TERM_RESIDUAL <- "Residual SD"

# Profile SD CIs for random intercepts and residual. lme4 may reorder grouping
# factors vs the formula (e.g. Trial before Subject), so .sig01 is NOT always
# Subject. Map .sig0k rows to group names via getME(model, "theta") order.
profile_random_sd_cis <- function(model) {
  ci_theta <- tryCatch(
    stats::confint(model, method = "profile", parm = "theta_", quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(ci_theta)) {
    return(NULL)
  }
  theta_names <- names(lme4::getME(model, "theta"))
  out <- list()
  for (i in seq_along(theta_names)) {
    rn <- sprintf(".sig%02d", i)
    if (!rn %in% rownames(ci_theta)) {
      next
    }
    grp <- sub("\\..*$", "", theta_names[i])
    out[[grp]] <- as.numeric(ci_theta[rn, ])
  }
  if (".sigma" %in% rownames(ci_theta)) {
    out[["__sigma__"]] <- as.numeric(ci_theta[".sigma", ])
  }
  out
}

variance_component_from_varcorr <- function(model, grp) {
  vc <- as.data.frame(lme4::VarCorr(model))
  re_row <- vc[vc$grp == grp & vc$var1 == "(Intercept)", , drop = FALSE]
  if (nrow(re_row) < 1 && identical(grp, "Subject")) {
    re_row <- vc[vc$grp == "ID" & vc$var1 == "(Intercept)", , drop = FALSE]
  }
  if (nrow(re_row) < 1) {
    return(NULL)
  }
  list(
    beta = re_row$vcov[1],
    SE = NA_real_,
    stat = NA_real_,
    p = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_
  )
}

# Variance-component point estimate plus profile CI for one grouping factor.
# profile_cis: optional output of profile_random_sd_cis() to avoid re-profiling.
variance_component_stats <- function(model, grp, profile_cis = NULL, include_rlrt = FALSE) {
  vc <- as.data.frame(lme4::VarCorr(model))
  re_row <- vc[vc$grp == grp & vc$var1 == "(Intercept)", , drop = FALSE]
  if (nrow(re_row) < 1 && identical(grp, "Subject")) {
    re_row <- vc[vc$grp == "ID" & vc$var1 == "(Intercept)", , drop = FALSE]
  }
  if (nrow(re_row) < 1) {
    return(NULL)
  }

  var_est <- re_row$vcov[1]
  sd_est <- re_row$sdcor[1]
  ci_low <- ci_high <- se_var <- stat <- p <- NA_real_

  if (is.null(profile_cis)) {
    profile_cis <- profile_random_sd_cis(model)
  }
  sd_ci <- NULL
  if (!is.null(profile_cis)) {
    sd_ci <- profile_cis[[grp]]
    if (is.null(sd_ci) && identical(grp, "Subject")) {
      sd_ci <- profile_cis[["ID"]]
    }
  }
  if (!is.null(sd_ci) && length(sd_ci) >= 2 && all(is.finite(sd_ci))) {
    ci_low <- sd_ci[1]^2
    ci_high <- sd_ci[2]^2
    sd_se <- (sd_ci[2] - sd_ci[1]) / (2 * stats::qnorm(0.975))
    se_var <- 2 * sd_est * sd_se
  }

  if (isTRUE(include_rlrt) && requireNamespace("RLRsim", quietly = TRUE)) {
    lrt <- tryCatch(
      suppressMessages(RLRsim::exactRLRT(model)),
      error = function(e) NULL
    )
    if (!is.null(lrt)) {
      stat <- as.numeric(lrt$stat)
      p <- as.numeric(lrt$p.value)
    }
  }

  list(
    beta = var_est,
    SE = se_var,
    stat = stat,
    p = p,
    CI_low = ci_low,
    CI_high = ci_high
  )
}

residual_sd_component <- function(model) {
  sig <- tryCatch(stats::sigma(model), error = function(e) NA_real_)
  if (!is.finite(sig)) {
    return(NULL)
  }
  list(
    beta = sig,
    SE = NA_real_,
    stat = NA_real_,
    p = NA_real_,
    CI_low = NA_real_,
    CI_high = NA_real_
  )
}

is_random_export_term <- function(term) {
  term %in% c(RE_TERM_SUBJECT, RE_TERM_TRIAL, RE_TERM_RESIDUAL) ||
    grepl("Var$|Cov$", as.character(term))
}

# Display token for adjusted p-value columns/headers (rendered as p with subscript adj).
PADJ_HEADER_TOKEN <- "padj"

padj_flextable_header <- function() {
  flextable::as_paragraph(
    flextable::as_chunk("p"),
    flextable::as_sub("adj")
  )
}

apply_padj_column_headers <- function(ft) {
  cols <- names(ft$body$dataset)
  j <- which(cols %in% c("padj", "p_adj"))
  if (length(j) == 0) {
    return(ft)
  }
  for (jj in j) {
    ft <- flextable::compose(ft, part = "header", j = jj, value = padj_flextable_header())
  }
  ft
}

padj_note_chunks <- function(suffix = " = FDR-adjusted p-value.", font_size = 10) {
  list(
    officer::ftext("p", officer::fp_text(font.size = font_size)),
    officer::ftext("adj", officer::fp_text(font.size = font_size, vertical.align = "subscript")),
    officer::ftext(suffix, officer::fp_text(font.size = font_size))
  )
}

note_label_chunks <- function(font_size = 10) {
  list(
    officer::ftext("Note.", officer::fp_text(bold = TRUE, font.size = font_size)),
    officer::ftext(" ", officer::fp_text(font.size = font_size))
  )
}

table_section_spacer_fpar <- function(font_size = 6) {
  officer::fpar(
    officer::ftext("\u00a0", officer::fp_text(font.size = font_size)),
    fp_p = officer::fp_par(line_spacing = 1, padding.top = 0, padding.bottom = 0)
  )
}

table_note_fp_par <- function() {
  officer::fp_par(line_spacing = 1, padding.top = 0, padding.bottom = 0)
}

assemble_table_note_fpar <- function(chunks) {
  do.call(officer::fpar, c(chunks, list(fp_p = table_note_fp_par())))
}

interaction_kept_from_lrt <- function(lrt_path) {
  if (!file.exists(lrt_path)) {
    return(FALSE)
  }
  lrt <- tryCatch(
    read.csv(lrt_path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(lrt) || nrow(lrt) < 1L || !"InteractionKept" %in% names(lrt)) {
    return(FALSE)
  }
  isTRUE(lrt$InteractionKept[1])
}

covariation_lrt_path <- function(
    predictor_key,
    task,
    stats_dir,
    exploratory = FALSE,
    trials = FALSE) {
  prefix <- if (isTRUE(exploratory)) {
    "AOC_ersd_exploratory_lrtTbl"
  } else {
    "AOC_ersd_lrtTbl"
  }
  if (isTRUE(trials)) {
    file.path(stats_dir, sprintf("%s_%s_%s_trials.csv", prefix, predictor_key, task))
  } else {
    file.path(stats_dir, sprintf("%s_%s_%s.csv", prefix, predictor_key, task))
  }
}

covariation_dual_interaction_retained <- function(
    predictor_key,
    stats_dir,
    exploratory = FALSE,
    trials = FALSE) {
  tasks <- c("nback", "sternberg")
  all(vapply(
    tasks,
    function(task) {
      interaction_kept_from_lrt(
        covariation_lrt_path(predictor_key, task, stats_dir, exploratory, trials)
      )
    },
    logical(1)
  ))
}

filter_interaction_full_model_specs <- function(specs, stats_dir) {
  Filter(function(spec) {
    if (!grepl("_full$", spec$id)) {
      return(TRUE)
    }
    if (is.null(spec$lrt)) {
      return(FALSE)
    }
    interaction_kept_from_lrt(file.path(stats_dir, spec$lrt))
  }, specs)
}

remove_interaction_full_exports <- function(
    gaze_sname,
    task,
    stats_dir,
    exploratory = FALSE,
    trials = FALSE) {
  prefix_fixed <- if (exploratory) "AOC_ersd_exploratory_fixed" else "AOC_ersd_fixed"
  prefix_con <- if (exploratory) "AOC_ersd_exploratory_contrasts" else "AOC_ersd_contrasts"
  if (isTRUE(trials)) {
    files <- c(
      file.path(stats_dir, sprintf("%s_%s_%s_full_trials.csv", prefix_fixed, gaze_sname, task)),
      file.path(stats_dir, sprintf("%s_%s_%s_full_trials.csv", prefix_con, gaze_sname, task))
    )
  } else {
    files <- c(
      file.path(stats_dir, sprintf("%s_%s_%s_full.csv", prefix_fixed, gaze_sname, task)),
      file.path(stats_dir, sprintf("%s_%s_%s_full.csv", prefix_con, gaze_sname, task))
    )
  }
  for (path in files) {
    if (file.exists(path)) {
      file.remove(path)
      message("Removed stale full-interaction export: ", basename(path))
    }
  }
  invisible(NULL)
}

remove_dual_covariation_variant_artifacts <- function(
    spec_id,
    variant,
    full_dir,
    combined_dir) {
  if (is.null(variant) || !identical(variant, "full")) {
    return(invisible(NULL))
  }
  id_suffix <- paste0("_", variant)
  files <- c(
    file.path(full_dir, paste0(spec_id, id_suffix, "_full_model_table.docx")),
    file.path(full_dir, paste0(spec_id, id_suffix, "_fixed.csv")),
    file.path(full_dir, paste0(spec_id, id_suffix, "_random.csv")),
    file.path(full_dir, paste0(spec_id, id_suffix, "_interaction_lrt.csv")),
    file.path(full_dir, paste0(spec_id, id_suffix, "_pairwise.csv")),
    file.path(combined_dir, paste0(spec_id, id_suffix, "_results_table.docx"))
  )
  for (path in files) {
    if (file.exists(path)) {
      file.remove(path)
      message("Removed stale covariation artifact: ", basename(path))
    }
  }
  invisible(NULL)
}
