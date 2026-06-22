# Export GLMM results (RDS) to CSV files for raincloud brackets and full model tables.

suppressPackageStartupMessages({
  library(broom.mixed)
  library(lme4)
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
source(file.path(export_helpers_dir, "AOC_stats_glmm_preprocess.R"))

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

compute_cohens_dz <- function(dat, dv, task, comparisons, padj_lookup = NULL) {
  rows <- list()
  for (cmp in comparisons) {
    g1 <- cmp[1]
    g2 <- cmp[2]
    d1 <- dat[dat$Load == g1, c("Subject", dv), drop = FALSE]
    d2 <- dat[dat$Load == g2, c("Subject", dv), drop = FALSE]
    names(d1)[2] <- "v1"
    names(d2)[2] <- "v2"
    m <- merge(d1, d2, by = "Subject")
    diffs <- m$v2 - m$v1
    diffs <- diffs[is.finite(diffs)]
    n <- length(diffs)
    if (n >= 2 && stats::sd(diffs) > 0) {
      md <- mean(diffs)
      sd_d <- stats::sd(diffs)
      dz <- md / sd_d
      tcrit <- stats::qt(0.975, df = n - 1)
      se_md <- sd_d / sqrt(n)
      ci_lo <- md - tcrit * se_md
      ci_hi <- md + tcrit * se_md
    } else {
      md <- dz <- ci_lo <- ci_hi <- NA_real_
    }
    padj <- NA_real_
    if (!is.null(padj_lookup)) {
      key <- paste(g1, g2, sep = "|")
      if (key %in% names(padj_lookup)) padj <- padj_lookup[[key]]
    }
    rows[[length(rows) + 1]] <- data.frame(
      Task = task,
      Variable = dv,
      Group1 = g1,
      Group2 = g2,
      N = n,
      MeanDiff = md,
      Cohens_dz = dz,
      CI95_low = ci_lo,
      CI95_high = ci_hi,
      p_adj = padj,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
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

export_fixed_effects <- function(res, task, stats_dir) {
  model <- res$model
  fx <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)

  out <- data.frame(
    Task = task,
    DV = res$dv,
    ModelLabel = res$model_label,
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
        Task = task, DV = res$dv, ModelLabel = res$model_label,
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

  fixed_path <- file.path(stats_dir, sprintf("AOC_mixedlm_fixed_%s_%s.csv", sname, task))
  write.csv(export_fixed_effects(res, task, stats_dir), fixed_path, row.names = FALSE)
  message("Saved fixed effects → ", basename(fixed_path))

  pw_df <- pairwise_to_export_df(
    res$pairwise, task, res$dv, res$model_label, res$n_obs
  )
  padj_map <- setNames(pw_df$p_adj, paste(pw_df$Group1, pw_df$Group2, sep = "|"))
  eff <- compute_cohens_dz(res$data, res$dv, task, comparisons, padj_map)

  list(pairwise = pw_df, effectsizes = eff)
}

export_interaction_model <- function(res, task, gaze_sname, stats_dir, exploratory = FALSE) {
  if (is.null(res)) return(invisible(NULL))

  prefix_fixed <- if (exploratory) "AOC_ersd_exploratory_fixed" else "AOC_ersd_fixed"
  prefix_lrt <- if (exploratory) "AOC_ersd_exploratory_lrtTbl" else "AOC_ersd_lrtTbl"
  prefix_con <- if (exploratory) "AOC_ersd_exploratory_contrasts" else "AOC_ersd_contrasts"

  fixed_path <- file.path(stats_dir, sprintf("%s_%s_%s.csv", prefix_fixed, gaze_sname, task))
  write.csv(export_fixed_effects(res, task, stats_dir), fixed_path, row.names = FALSE)
  message("Saved ERSD fixed effects → ", basename(fixed_path))

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
    stringsAsFactors = FALSE
  )
  lrt_path <- file.path(stats_dir, sprintf("%s_%s_%s.csv", prefix_lrt, gaze_sname, task))
  write.csv(lrt_tbl, lrt_path, row.names = FALSE)
  message("Saved interaction LRT (drop1) → ", basename(lrt_path))

  pw_df <- pairwise_to_export_df(
    res$pairwise, task, res$dv, res$model_label, res$n_obs
  )
  pw_df$GazePredictor <- res$gaze_var
  con_path <- file.path(stats_dir, sprintf("%s_%s_%s.csv", prefix_con, gaze_sname, task))
  write.csv(pw_df, con_path, row.names = FALSE)
  message("Saved ERSD contrasts → ", basename(con_path))
}

export_glmm_results <- function(results, stats_dir = AOC_STATS_DIR) {
  task <- results$task
  comparisons <- results$task_config$comparisons

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
    pw_path <- file.path(stats_dir, sprintf("AOC_pairwise_mixedlm_%s.csv", task))
    write.csv(do.call(rbind, all_pw), pw_path, row.names = FALSE)
    message("Saved pooled pairwise contrasts → ", basename(pw_path))
  }

  if (length(all_eff) > 0) {
    eff_path <- file.path(stats_dir, sprintf("AOC_pairwise_effectsizes_%s.csv", task))
    write.csv(do.call(rbind, all_eff), eff_path, row.names = FALSE)
    message("Saved pairwise effect sizes → ", basename(eff_path))
  }

  invisible(NULL)
}

export_all_tasks <- function(stats_dir = AOC_STATS_DIR, tasks = c("nback", "sternberg")) {
  if (!dir.exists(stats_dir)) {
    stop("Stats directory not found: ", stats_dir)
  }
  for (task in tasks) {
    rds_path <- file.path(stats_dir, sprintf("AOC_glmm_results_%s.rds", task))
    if (!file.exists(rds_path)) {
      warning("Missing ", rds_path, " — run AOC_stats_glmm_", task, ".R first.")
      next
    }
    message("\n=== Exporting ", task, " ===")
    results <- readRDS(rds_path)
    export_glmm_results(results, stats_dir)
  }
  invisible(NULL)
}
