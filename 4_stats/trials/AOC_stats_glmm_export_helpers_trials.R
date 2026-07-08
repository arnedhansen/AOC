# Export trial-level GLMM results (RDS) to CSV files.
# Outputs go under data/stats/trials (separate from condition-level stats).

suppressPackageStartupMessages({
  library(broom.mixed)
  library(lme4)
  library(emmeans)
  if (!requireNamespace("RLRsim", quietly = TRUE)) {
    stop("Package 'RLRsim' required for random-effects inference. Install with install.packages('RLRsim').")
  }
  library(RLRsim)
})

export_helpers_dir <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
}
source(file.path(export_helpers_dir, "AOC_stats_glmm_preprocess_trials.R"))

lr_effect_sizes <- function(LR, df_diff, n_obs) {
  if (!is.finite(LR) || LR <= 0 || !is.finite(n_obs) || n_obs <= 0) {
    return(c(R2_LR = NA_real_, f2_LR = NA_real_))
  }
  R2_LR <- 1 - exp(-LR / n_obs)
  f2_LR <- if (R2_LR < 1) R2_LR / (1 - R2_LR) else NA_real_
  c(R2_LR = R2_LR, f2_LR = f2_LR)
}

LOAD_SNAME <- c(
  Accuracy = "acc",
  ReactionTime = "rt",
  GazeDeviation = "gazedev",
  MSRate = "ms",
  GazeDeviationBL = "gazedev_bl",
  MSRateBL = "ms_bl",
  ERSD = "ersd"
)

INTERACTION_SNAME <- c(
  GazeDeviation = "gazedeviation",
  MSRate = "msrate",
  GazeDeviationBL = "gazedev_bl",
  MSRateBL = "ms_bl"
)

pairwise_to_export_df <- function(pw, task, variable, model_label, n_obs) {
  z_vals <- if ("z.ratio" %in% names(pw)) pw$z.ratio else pw$t.ratio
  if ("lower.CL" %in% names(pw)) {
    ci_lo <- pw$lower.CL
    ci_hi <- pw$upper.CL
  } else {
    ci_lo <- pw$estimate - 1.96 * pw$SE
    ci_hi <- pw$estimate + 1.96 * pw$SE
  }
  p_vals <- if ("p.value" %in% names(pw)) pw$p.value else pw$p
  p_adj_vals <- if ("adj.p.value" %in% names(pw)) pw$adj.p.value else p_vals

  split_contrast <- strsplit(as.character(pw$contrast), " - ", fixed = TRUE)
  strip_lvl <- function(s) gsub("[()]", "", trimws(s))
  group2 <- vapply(split_contrast, function(x) strip_lvl(x[1]), character(1))
  group1 <- vapply(split_contrast, function(x) strip_lvl(x[2]), character(1))

  data.frame(
    Task = task,
    Variable = variable,
    ModelLabel = model_label,
    N_obs = n_obs,
    Group1 = group1,
    Group2 = group2,
    Estimate = pw$estimate,
    SE = pw$SE,
    z = z_vals,
    p = p_vals,
    CI95_low = ci_lo,
    CI95_high = ci_hi,
    p_adj = p_adj_vals,
    stringsAsFactors = FALSE
  )
}

# Trial-level effect sizes from the fitted mixed model (no subject averaging).
# Cohen's d = pairwise Load contrast / residual SD (emmeans::eff_size).
compute_cohens_d_trials <- function(res, task, comparisons, padj_lookup = NULL) {
  if (is.null(res) || is.null(res$model) || is.null(res$data)) return(NULL)

  model <- res$model
  dat <- res$data
  dv <- res$dv
  n_obs <- nrow(dat)
  n_subj <- length(unique(dat$Subject))

  sig <- tryCatch(stats::sigma(model), error = function(e) NA_real_)
  if (!is.finite(sig) || sig <= 0) return(NULL)

  edf <- tryCatch(stats::df.residual(model), error = function(e) NA_real_)
  if (!is.finite(edf) || edf <= 0) {
    edf <- max(n_obs - length(unique(dat$Load)), 1)
  }

  emm <- tryCatch(
    emmeans::emmeans(model, ~Load, data = dat),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)

  eff <- tryCatch(
    as.data.frame(emmeans::eff_size(emm, sigma = sig, edf = edf)),
    error = function(e) NULL
  )
  if (is.null(eff) || nrow(eff) == 0) return(NULL)

  # Pairwise means / CIs from the same emmeans object (trial-level contrasts).
  pw <- tryCatch(
    as.data.frame(pairs(emm, adjust = "fdr")),
    error = function(e) NULL
  )

  strip_lvl <- function(s) gsub("[()]", "", trimws(s))
  parse_contrast <- function(contrast) {
    parts <- strsplit(as.character(contrast), " - ", fixed = TRUE)[[1]]
    list(g2 = strip_lvl(parts[1]), g1 = strip_lvl(parts[2]))
  }

  rows <- list()
  for (i in seq_len(nrow(eff))) {
    parsed <- parse_contrast(eff$contrast[i])
    g1 <- parsed$g1
    g2 <- parsed$g2

    md <- ci_lo <- ci_hi <- NA_real_
    if (!is.null(pw) && "contrast" %in% names(pw)) {
      j <- which(as.character(pw$contrast) == as.character(eff$contrast[i]))
      if (length(j) == 1) {
        md <- pw$estimate[j]
        if ("lower.CL" %in% names(pw)) {
          ci_lo <- pw$lower.CL[j]
          ci_hi <- pw$upper.CL[j]
        } else if ("SE" %in% names(pw)) {
          ci_lo <- md - 1.96 * pw$SE[j]
          ci_hi <- md + 1.96 * pw$SE[j]
        }
      }
    }

    padj <- NA_real_
    if (!is.null(padj_lookup)) {
      key <- paste(g1, g2, sep = "|")
      if (key %in% names(padj_lookup)) {
        padj <- padj_lookup[[key]]
      } else {
        key_rev <- paste(g2, g1, sep = "|")
        if (key_rev %in% names(padj_lookup)) padj <- padj_lookup[[key_rev]]
      }
    }

    d_col <- if ("effect.size" %in% names(eff)) {
      eff$effect.size[i]
    } else if ("eff.size" %in% names(eff)) {
      eff$eff.size[i]
    } else {
      NA_real_
    }

    rows[[length(rows) + 1]] <- data.frame(
      Task = task,
      Variable = dv,
      Group1 = g1,
      Group2 = g2,
      N = n_subj,
      N_obs = n_obs,
      MeanDiff = md,
      Cohens_dz = d_col,
      CI95_low = ci_lo,
      CI95_high = ci_hi,
      p_adj = padj,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)

  # Keep only requested load comparisons when provided.
  if (!is.null(comparisons) && length(comparisons) > 0) {
    keep <- vapply(seq_len(nrow(out)), function(i) {
      any(vapply(comparisons, function(cmp) {
        (out$Group1[i] == cmp[1] && out$Group2[i] == cmp[2]) ||
          (out$Group1[i] == cmp[2] && out$Group2[i] == cmp[1])
      }, logical(1)))
    }, logical(1))
    out <- out[keep, , drop = FALSE]
  }
  out
}

random_intercept_stats <- function(model) {
  vc <- as.data.frame(VarCorr(model))
  re_row <- vc[
    vc$grp %in% c("Subject", "ID") & vc$var1 == "(Intercept)",
    ,
    drop = FALSE
  ]
  if (nrow(re_row) < 1) {
    return(NULL)
  }

  var_est <- re_row$vcov[1]
  sd_est <- re_row$sdcor[1]
  ci_low <- ci_high <- se_var <- stat <- p <- NA_real_

  ci_theta <- tryCatch(
    confint(model, method = "profile", parm = "theta_", quiet = TRUE),
    error = function(e) NULL
  )
  if (!is.null(ci_theta) && ".sig01" %in% rownames(ci_theta)) {
    sd_ci <- ci_theta[".sig01", ]
    ci_low <- sd_ci[1]^2
    ci_high <- sd_ci[2]^2
    sd_se <- (sd_ci[2] - sd_ci[1]) / (2 * stats::qnorm(0.975))
    se_var <- 2 * sd_est * sd_se
  }

  lrt <- tryCatch(
    suppressMessages(exactRLRT(model)),
    error = function(e) NULL
  )
  if (!is.null(lrt)) {
    stat <- as.numeric(lrt$stat)
    p <- as.numeric(lrt$p.value)
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

export_fixed_effects <- function(res, task, stats_dir, model = NULL, model_label = NULL) {
  if (is.null(model)) model <- res$model
  if (is.null(model_label)) model_label <- res$model_label
  fx <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)

  out <- data.frame(
    Task = task,
    DV = res$dv,
    ModelLabel = model_label,
    Term = fx$term,
    beta = fx$estimate,
    SE = fx$std.error,
    stat = fx$statistic,
    p = fx$p.value,
    CI_low = fx$conf.low,
    CI_high = fx$conf.high,
    N_obs = res$n_obs,
    stringsAsFactors = FALSE
  )

  re_stats <- random_intercept_stats(model)
  if (!is.null(re_stats)) {
    out <- rbind(
      out,
      data.frame(
        Task = task, DV = res$dv, ModelLabel = model_label,
        Term = "Group Var",
        beta = re_stats$beta,
        SE = re_stats$SE,
        stat = re_stats$stat,
        p = re_stats$p,
        CI_low = re_stats$CI_low,
        CI_high = re_stats$CI_high,
        N_obs = res$n_obs,
        stringsAsFactors = FALSE
      )
    )
  }

  out
}

export_load_model <- function(res, task, sname, stats_dir, comparisons) {
  if (is.null(res)) return(invisible(NULL))

  fixed_path <- file.path(
    stats_dir, sprintf("AOC_mixedlm_fixed_%s_%s_trials.csv", sname, task)
  )
  write.csv(export_fixed_effects(res, task, stats_dir), fixed_path, row.names = FALSE)
  message("Saved fixed effects → ", basename(fixed_path))

  pw_df <- pairwise_to_export_df(
    res$pairwise, task, res$dv, res$model_label, res$n_obs
  )
  padj_map <- setNames(pw_df$p_adj, paste(pw_df$Group1, pw_df$Group2, sep = "|"))
  eff <- compute_cohens_d_trials(res, task, comparisons, padj_map)

  list(pairwise = pw_df, effectsizes = eff)
}

export_interaction_model <- function(res, task, gaze_sname, stats_dir, exploratory = FALSE) {
  if (is.null(res)) return(invisible(NULL))

  prefix_fixed <- if (exploratory) "AOC_ersd_exploratory_fixed" else "AOC_ersd_fixed"
  prefix_lrt <- if (exploratory) "AOC_ersd_exploratory_lrtTbl" else "AOC_ersd_lrtTbl"
  prefix_con <- if (exploratory) "AOC_ersd_exploratory_contrasts" else "AOC_ersd_contrasts"

  # Selected model (manuscript default: full if interaction kept, else additive).
  fixed_path <- file.path(
    stats_dir, sprintf("%s_%s_%s_trials.csv", prefix_fixed, gaze_sname, task)
  )
  write.csv(export_fixed_effects(res, task, stats_dir), fixed_path, row.names = FALSE)
  message("Saved ERSD fixed effects → ", basename(fixed_path))

  label_full <- if (!is.null(res$model_label_full)) {
    res$model_label_full
  } else {
    res$model_label
  }
  label_red <- if (!is.null(res$model_label_reduced)) {
    res$model_label_reduced
  } else {
    res$model_label
  }

  model_full <- if (!is.null(res$model_full)) res$model_full else res$model
  model_red <- if (!is.null(res$model_reduced)) res$model_reduced else res$model

  fixed_full_path <- file.path(
    stats_dir, sprintf("%s_%s_%s_full_trials.csv", prefix_fixed, gaze_sname, task)
  )
  write.csv(
    export_fixed_effects(res, task, stats_dir, model = model_full, model_label = label_full),
    fixed_full_path,
    row.names = FALSE
  )
  message("Saved ERSD fixed effects (full) → ", basename(fixed_full_path))

  fixed_red_path <- file.path(
    stats_dir, sprintf("%s_%s_%s_reduced_trials.csv", prefix_fixed, gaze_sname, task)
  )
  write.csv(
    export_fixed_effects(res, task, stats_dir, model = model_red, model_label = label_red),
    fixed_red_path,
    row.names = FALSE
  )
  message("Saved ERSD fixed effects (reduced) → ", basename(fixed_red_path))

  lrt <- res$lrt_interaction
  eff <- lr_effect_sizes(lrt$LR_Chi2, lrt$df, res$n_obs)
  lrt_tbl <- data.frame(
    Term = lrt$Term,
    df = lrt$df,
    LR_Chi2 = lrt$LR_Chi2,
    p = lrt$p,
    R2_LR = eff["R2_LR"],
    f2_LR = eff["f2_LR"],
    N_obs = res$n_obs,
    Task = task,
    DV = res$dv,
    GazePredictor = res$gaze_var,
    InteractionKept = isTRUE(res$interaction_kept),
    stringsAsFactors = FALSE
  )
  lrt_path <- file.path(
    stats_dir, sprintf("%s_%s_%s_trials.csv", prefix_lrt, gaze_sname, task)
  )
  write.csv(lrt_tbl, lrt_path, row.names = FALSE)
  message("Saved interaction LRT (drop1) → ", basename(lrt_path))

  write_contrasts <- function(pw, model_label, suffix = NULL) {
    pw_df <- pairwise_to_export_df(pw, task, res$dv, model_label, res$n_obs)
    pw_df$GazePredictor <- res$gaze_var
    fname <- if (is.null(suffix)) {
      sprintf("%s_%s_%s_trials.csv", prefix_con, gaze_sname, task)
    } else {
      sprintf("%s_%s_%s_%s_trials.csv", prefix_con, gaze_sname, task, suffix)
    }
    con_path <- file.path(stats_dir, fname)
    write.csv(pw_df, con_path, row.names = FALSE)
    message("Saved ERSD contrasts", if (!is.null(suffix)) paste0(" (", suffix, ")") else "", " → ", basename(con_path))
  }

  write_contrasts(res$pairwise, res$model_label)
  pw_full <- if (!is.null(res$pairwise_full)) res$pairwise_full else res$pairwise
  pw_red <- if (!is.null(res$pairwise_reduced)) res$pairwise_reduced else res$pairwise
  write_contrasts(pw_full, label_full, "full")
  write_contrasts(pw_red, label_red, "reduced")
}

export_glmm_results <- function(results, stats_dir = AOC_STATS_DIR_TRIALS) {
  task <- results$task
  comparisons <- results$task_config$comparisons

  if (!dir.exists(stats_dir)) {
    dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
  }

  all_pw <- list()
  all_eff <- list()

  for (dv in names(results$load_models)) {
    sname <- LOAD_SNAME[[dv]]
    if (is.null(sname)) next
    out <- export_load_model(results$load_models[[dv]], task, sname, stats_dir, comparisons)
    if (!is.null(out)) {
      all_pw[[length(all_pw) + 1]] <- out$pairwise
      if (!is.null(out$effectsizes)) {
        all_eff[[length(all_eff) + 1]] <- out$effectsizes
      }
    }
  }

  for (nm in names(results$interaction_confirmatory)) {
    gaze_sname <- INTERACTION_SNAME[[nm]]
    export_interaction_model(
      results$interaction_confirmatory[[nm]], task, gaze_sname, stats_dir, exploratory = FALSE
    )
  }

  for (nm in names(results$interaction_exploratory)) {
    gaze_sname <- INTERACTION_SNAME[[nm]]
    export_interaction_model(
      results$interaction_exploratory[[nm]], task, gaze_sname, stats_dir, exploratory = TRUE
    )
  }

  if (length(all_pw) > 0) {
    pw_path <- file.path(stats_dir, sprintf("AOC_pairwise_mixedlm_%s_trials.csv", task))
    write.csv(do.call(rbind, all_pw), pw_path, row.names = FALSE)
    message("Saved pooled pairwise contrasts → ", basename(pw_path))
  }

  if (length(all_eff) > 0) {
    eff_path <- file.path(stats_dir, sprintf("AOC_pairwise_effectsizes_%s_trials.csv", task))
    write.csv(do.call(rbind, all_eff), eff_path, row.names = FALSE)
    message("Saved pairwise effect sizes → ", basename(eff_path))
  }

  invisible(NULL)
}

export_all_tasks_trials <- function(stats_dir = AOC_STATS_DIR_TRIALS, tasks = c("nback", "sternberg")) {
  if (!dir.exists(stats_dir)) {
    stop("Trial stats directory not found: ", stats_dir)
  }
  for (task in tasks) {
    rds_path <- file.path(stats_dir, sprintf("AOC_glmm_results_%s_trials.rds", task))
    if (!file.exists(rds_path)) {
      warning(
        "Missing ", rds_path,
        " — run AOC_stats_glmm_", task, "_trials.R first."
      )
      next
    }
    message("\n=== Exporting trial-level ", task, " ===")
    results <- readRDS(rds_path)
    export_glmm_results(results, stats_dir)
  }
  invisible(NULL)
}
