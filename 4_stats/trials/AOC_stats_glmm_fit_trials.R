# Trial-level GLMM fitting (lme4). Separate from condition-level AOC_stats_glmm_fit.R.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(car)
  library(emmeans)
})

fit_script_dir_trials <- {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("--file=", "", args[grep("^--file=", args)])
  if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg, winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
}
source(file.path(fit_script_dir_trials, "AOC_stats_glmm_preprocess_trials.R"))

EEG_DV_TRIALS <- "ERSD"

fit_lmer_gaussian_trials <- function(formula, data) {
  mod <- tryCatch(
    lmer(
      formula, data = data, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) {
      lmer(formula, data = data, REML = FALSE)
    }
  )

  fam <- "gaussian"
  res <- residuals(mod)
  res <- res[is.finite(res)]
  if (length(res) >= 20) {
    sw <- tryCatch(shapiro.test(res), error = function(e) NULL)
    sk <- if (stats::sd(res) > 0) {
      mean(((res - mean(res)) / stats::sd(res))^3)
    } else {
      0
    }
    if (!is.null(sw) && sw$p.value < 0.05 && is.finite(sk) && abs(sk) > 1) {
      mod_gamma <- tryCatch(
        glmer(
          formula, data = data, family = Gamma(link = "log"), REML = FALSE,
          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        ),
        error = function(e) NULL
      )
      if (!is.null(mod_gamma)) {
        mod <- mod_gamma
        fam <- "gamma"
      }
    }
  }

  mod@call$data <- data
  attr(mod, "family_used") <- fam
  mod
}

emmeans_load_pairs_trials <- function(model, data, at = NULL) {
  emm_args <- list(object = model, specs = ~Load, data = data)
  if (!is.null(at)) emm_args$at <- at
  emm <- do.call(emmeans, emm_args)
  as.data.frame(pairs(emm, adjust = "fdr"))
}

fit_dv_load_trials <- function(dat, dv, model_label) {
  sub <- dat[!is.na(dat[[dv]]), c("Subject", "Load", dv), drop = FALSE]
  if (nrow(sub) < 10 || length(unique(sub$Subject)) < 2) {
    warning(sprintf("Skipping %s ~ Load: insufficient data.", dv))
    return(NULL)
  }

  form <- as.formula(sprintf("%s ~ Load + (1 | Subject)", dv))
  model <- fit_lmer_gaussian_trials(form, sub)
  anova_tbl <- car::Anova(model, type = 2)
  confint_tbl <- tryCatch(
    confint(model, method = "Wald", parm = "beta_"),
    error = function(e) confint(model, method = "profile", parm = "beta_")
  )
  pairwise <- emmeans_load_pairs_trials(model, sub)

  list(
    kind = "load",
    dv = dv,
    model_label = model_label,
    model = model,
    family = attr(model, "family_used"),
    anova = anova_tbl,
    anova_type = 2L,
    confint = confint_tbl,
    pairwise = pairwise,
    n_obs = nrow(sub),
    data = sub
  )
}

fit_eeg_gaze_interaction_trials <- function(dat, gaze_var, model_label_full, exploratory = FALSE) {
  cols <- c("Subject", "Load", EEG_DV_TRIALS, gaze_var)
  sub <- dat[complete.cases(dat[, cols]), cols, drop = FALSE]
  if (nrow(sub) < 10 || length(unique(sub$Subject)) < 2) {
    warning(sprintf("Skipping %s ~ %s * Load: insufficient data.", EEG_DV_TRIALS, gaze_var))
    return(NULL)
  }

  gaze_c <- paste0(gaze_var, "_c")
  sub[[gaze_c]] <- sub[[gaze_var]] - mean(sub[[gaze_var]], na.rm = TRUE)
  sub$Load <- factor(sub$Load, levels = levels(dat$Load), ordered = TRUE)

  form_full <- as.formula(sprintf("%s ~ %s * Load + (1 | Subject)", EEG_DV_TRIALS, gaze_c))
  form_red <- as.formula(sprintf("%s ~ %s + Load + (1 | Subject)", EEG_DV_TRIALS, gaze_c))

  full <- fit_lmer_gaussian_trials(form_full, sub)
  red <- fit_lmer_gaussian_trials(form_red, sub)

  lrt_int <- anova(red, full)
  int_p <- as.numeric(lrt_int$`Pr(>Chisq)`[2])
  int_lr <- as.numeric(lrt_int$Chisq[2])
  int_df <- as.numeric(lrt_int$Df[2])
  interaction_kept <- isTRUE(!is.na(int_p) && int_p < 0.05)
  final <- if (interaction_kept) full else red

  drop1_tbl <- tryCatch(as.data.frame(drop1(full)), error = function(e) NULL)

  anova_type <- if (interaction_kept) 3L else 2L
  anova_tbl <- car::Anova(final, type = anova_type)
  confint_tbl <- tryCatch(
    confint(final, method = "Wald", parm = "beta_"),
    error = function(e) confint(final, method = "profile", parm = "beta_")
  )

  at_list <- setNames(list(0), gaze_c)
  pairwise <- emmeans_load_pairs_trials(final, sub, at = at_list)

  model_label <- if (interaction_kept) {
    model_label_full
  } else {
    sprintf("%s ~ %s + Load + (1|Subject)", EEG_DV_TRIALS, gaze_var)
  }

  lrt_interaction <- data.frame(
    Term = "Interaction",
    df = int_df,
    LR_Chi2 = int_lr,
    p = int_p,
    stringsAsFactors = FALSE
  )

  list(
    kind = "interaction",
    dv = EEG_DV_TRIALS,
    gaze_var = gaze_var,
    exploratory = exploratory,
    model_label = model_label,
    model_label_full = model_label_full,
    model_full = full,
    model = final,
    drop1 = drop1_tbl,
    lrt_interaction = lrt_interaction,
    interaction_kept = interaction_kept,
    family = attr(final, "family_used"),
    anova = anova_tbl,
    anova_type = anova_type,
    confint = confint_tbl,
    pairwise = pairwise,
    n_obs = nrow(sub),
    data = sub
  )
}

run_paper_glmm_trials <- function(task_config) {
  task <- task_config$name
  dat <- load_task_data_trials(task_config)

  message("Task (trials): ", task, " | N rows: ", nrow(dat))

  load_models <- list(
    Accuracy = fit_dv_load_trials(dat, "Accuracy", "Accuracy ~ Load + (1|Subject)"),
    ReactionTime = fit_dv_load_trials(dat, "ReactionTime", "ReactionTime ~ Load + (1|Subject)"),
    GazeDeviation = fit_dv_load_trials(dat, "GazeDeviation", "GazeDeviation ~ Load + (1|Subject)"),
    MSRate = fit_dv_load_trials(dat, "MSRate", "MSRate ~ Load + (1|Subject)"),
    GazeDeviationBL = fit_dv_load_trials(
      dat, "GazeDeviationBL", "GazeDeviationBL ~ Load + (1|Subject)"
    ),
    MSRateBL = fit_dv_load_trials(dat, "MSRateBL", "MSRateBL ~ Load + (1|Subject)"),
    ERSD = fit_dv_load_trials(dat, EEG_DV_TRIALS, sprintf("%s ~ Load + (1|Subject)", EEG_DV_TRIALS))
  )

  interaction_confirmatory <- list(
    GazeDeviation = fit_eeg_gaze_interaction_trials(
      dat, "GazeDeviation",
      sprintf("%s ~ GazeDeviation * Load + (1|Subject)", EEG_DV_TRIALS),
      exploratory = FALSE
    ),
    MSRate = fit_eeg_gaze_interaction_trials(
      dat, "MSRate",
      sprintf("%s ~ MSRate * Load + (1|Subject)", EEG_DV_TRIALS),
      exploratory = FALSE
    )
  )

  interaction_exploratory <- list(
    GazeDeviationBL = fit_eeg_gaze_interaction_trials(
      dat, "GazeDeviationBL",
      sprintf("%s ~ GazeDeviationBL * Load + (1|Subject)", EEG_DV_TRIALS),
      exploratory = TRUE
    ),
    MSRateBL = fit_eeg_gaze_interaction_trials(
      dat, "MSRateBL",
      sprintf("%s ~ MSRateBL * Load + (1|Subject)", EEG_DV_TRIALS),
      exploratory = TRUE
    )
  )

  results <- list(
    task = task,
    task_config = task_config,
    load_models = load_models,
    interaction_confirmatory = interaction_confirmatory,
    interaction_exploratory = interaction_exploratory,
    fitted_at = Sys.time(),
    level = "trials"
  )

  for (nm in names(load_models)) {
    res <- load_models[[nm]]
    if (!is.null(res)) {
      message(sprintf(
        "  %s ~ Load: family=%s, Anova type II",
        res$dv, res$family
      ))
    }
  }

  for (block in list(interaction_confirmatory, interaction_exploratory)) {
    for (nm in names(block)) {
      res <- block[[nm]]
      if (!is.null(res)) {
        message(sprintf(
          "  %s ~ %s * Load: interaction_kept=%s, Anova type %d, drop1 p=%.4f",
          res$dv, res$gaze_var, res$interaction_kept, res$anova_type,
          res$lrt_interaction$p[1]
        ))
      }
    }
  }

  invisible(results)
}

save_glmm_results_trials <- function(results, stats_dir = AOC_STATS_DIR_TRIALS) {
  if (!dir.exists(stats_dir)) {
    dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
  }
  path <- file.path(stats_dir, sprintf("AOC_glmm_results_%s_trials.rds", results$task))
  saveRDS(results, path)
  message("Saved trial-level GLMM results → ", path)
  invisible(path)
}
