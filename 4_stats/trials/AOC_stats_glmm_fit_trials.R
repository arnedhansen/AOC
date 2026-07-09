# Trial-level GLMM fitting (lme4). Separate from condition-level AOC_stats_glmm_fit.R.
# Random structure: (1|Subject) + (1|Trial) when Trial is present (crossed session order).

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
TRIAL_RE_LABEL <- "(1|Subject) + (1|Trial)"

has_trial_re_trials <- function(data) {
  "Trial" %in% names(data) && length(unique(stats::na.omit(data$Trial))) > 1
}

re_terms_trials <- function(data) {
  if (has_trial_re_trials(data)) {
    "(1 | Subject) + (1 | Trial)"
  } else {
    "(1 | Subject)"
  }
}

re_label_trials <- function(data) {
  if (has_trial_re_trials(data)) TRIAL_RE_LABEL else "(1|Subject)"
}

fit_lmer_gaussian_trials <- function(formula, data) {
  mod <- tryCatch(
    lmer(
      formula, data = data, REML = FALSE,
      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) {
      ftxt <- paste(deparse(formula), collapse = " ")
      if (grepl("\\(1 \\| Trial\\)", ftxt)) {
        f_fallback <- stats::as.formula(gsub("\\s*\\+\\s*\\(1 \\| Trial\\)", "", ftxt))
        warning("Falling back to (1|Subject) only for: ", ftxt, " | reason: ", conditionMessage(e))
        lmer(
          f_fallback, data = data, REML = FALSE,
          control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
        )
      } else {
        stop(e)
      }
    }
  )

  fam <- "gaussian"
  res <- residuals(mod)
  res <- res[is.finite(res)]
  # Shapiro is unreliable for large N; skip gamma switch when n > 5000.
  if (length(res) >= 20 && length(res) <= 5000) {
    sw <- tryCatch(shapiro.test(res), error = function(e) NULL)
    sk <- if (stats::sd(res) > 0) {
      mean(((res - mean(res)) / stats::sd(res))^3)
    } else {
      0
    }
    if (!is.null(sw) && sw$p.value < 0.05 && is.finite(sk) && abs(sk) > 1) {
      mod_gamma <- tryCatch(
        glmer(
          formula, data = data, family = Gamma(link = "log"),
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

fit_dv_load_trials <- function(dat, dv, model_label = NULL) {
  cols <- c("Subject", "Load", dv)
  if ("Trial" %in% names(dat)) cols <- c(cols, "Trial")
  sub <- dat[!is.na(dat[[dv]]), cols, drop = FALSE]
  if (nrow(sub) < 10 || length(unique(sub$Subject)) < 2) {
    warning(sprintf("Skipping %s ~ Load: insufficient data.", dv))
    return(NULL)
  }

  re <- re_terms_trials(sub)
  if (is.null(model_label)) {
    model_label <- sprintf("%s ~ Load + %s", dv, re_label_trials(sub))
  }
  form <- as.formula(sprintf("%s ~ Load + %s", dv, re))
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

fit_eeg_gaze_interaction_trials <- function(dat, gaze_var, model_label_full = NULL, exploratory = FALSE) {
  cols <- c("Subject", "Load", EEG_DV_TRIALS, gaze_var)
  if ("Trial" %in% names(dat)) cols <- c(cols, "Trial")
  sub <- dat[complete.cases(dat[, cols]), cols, drop = FALSE]
  if (nrow(sub) < 10 || length(unique(sub$Subject)) < 2) {
    warning(sprintf("Skipping %s ~ %s * Load: insufficient data.", EEG_DV_TRIALS, gaze_var))
    return(NULL)
  }

  gaze_c <- paste0(gaze_var, "_c")
  sub[[gaze_c]] <- sub[[gaze_var]] - mean(sub[[gaze_var]], na.rm = TRUE)
  sub$Load <- factor(sub$Load, levels = levels(dat$Load), ordered = FALSE)

  re <- re_terms_trials(sub)
  re_lab <- re_label_trials(sub)
  if (is.null(model_label_full)) {
    model_label_full <- sprintf("%s ~ %s * Load + %s", EEG_DV_TRIALS, gaze_var, re_lab)
  }

  form_full <- as.formula(sprintf("%s ~ %s * Load + %s", EEG_DV_TRIALS, gaze_c, re))
  form_red <- as.formula(sprintf("%s ~ %s + Load + %s", EEG_DV_TRIALS, gaze_c, re))

  full <- fit_lmer_gaussian_trials(form_full, sub)
  red <- fit_lmer_gaussian_trials(form_red, sub)

  lrt_int <- anova(red, full)
  int_p <- as.numeric(lrt_int$`Pr(>Chisq)`[2])
  int_lr <- as.numeric(lrt_int$Chisq[2])
  int_df <- as.numeric(lrt_int$Df[2])
  interaction_kept <- isTRUE(!is.na(int_p) && int_p < 0.05)
  final <- if (interaction_kept) full else red

  drop1_tbl <- tryCatch(as.data.frame(drop1(full)), error = function(e) NULL)

  at_list <- setNames(list(0), gaze_c)
  model_label_reduced <- sprintf("%s ~ %s + Load + %s", EEG_DV_TRIALS, gaze_var, re_lab)

  # Always summarize both models so tables can report full and additive fits.
  summarize_fit <- function(model, anova_type) {
    anova_tbl <- car::Anova(model, type = anova_type)
    confint_tbl <- tryCatch(
      confint(model, method = "Wald", parm = "beta_"),
      error = function(e) confint(model, method = "profile", parm = "beta_")
    )
    pairwise <- emmeans_load_pairs_trials(model, sub, at = at_list)
    list(anova = anova_tbl, anova_type = anova_type, confint = confint_tbl, pairwise = pairwise)
  }

  sum_full <- summarize_fit(full, 3L)
  sum_red <- summarize_fit(red, 2L)
  sum_final <- if (interaction_kept) sum_full else sum_red

  model_label <- if (interaction_kept) model_label_full else model_label_reduced

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
    model_label_reduced = model_label_reduced,
    model_full = full,
    model_reduced = red,
    model = final,
    drop1 = drop1_tbl,
    lrt_interaction = lrt_interaction,
    interaction_kept = interaction_kept,
    family = attr(final, "family_used"),
    family_full = attr(full, "family_used"),
    family_reduced = attr(red, "family_used"),
    anova = sum_final$anova,
    anova_type = sum_final$anova_type,
    anova_full = sum_full$anova,
    anova_type_full = sum_full$anova_type,
    anova_reduced = sum_red$anova,
    anova_type_reduced = sum_red$anova_type,
    confint = sum_final$confint,
    confint_full = sum_full$confint,
    confint_reduced = sum_red$confint,
    pairwise = sum_final$pairwise,
    pairwise_full = sum_full$pairwise,
    pairwise_reduced = sum_red$pairwise,
    n_obs = nrow(sub),
    data = sub
  )
}

run_paper_glmm_trials <- function(task_config) {
  task <- task_config$name
  dat <- load_task_data_trials(task_config)

  message(
    "Task (trials): ", task, " | N rows: ", nrow(dat),
    " | RE: ", re_label_trials(dat),
    " | IQR: ", IQR_K_TRIALS, "× within Subject×Load"
  )

  load_models <- list(
    Accuracy = fit_dv_load_trials(dat, "Accuracy"),
    ReactionTime = fit_dv_load_trials(dat, "ReactionTime"),
    GazeDeviation = fit_dv_load_trials(dat, "GazeDeviation"),
    MSRate = fit_dv_load_trials(dat, "MSRate"),
    GazeDeviationBL = fit_dv_load_trials(dat, "GazeDeviationBL"),
    MSRateBL = fit_dv_load_trials(dat, "MSRateBL"),
    ERSD = fit_dv_load_trials(dat, EEG_DV_TRIALS)
  )

  interaction_confirmatory <- list(
    GazeDeviation = fit_eeg_gaze_interaction_trials(
      dat, "GazeDeviation", exploratory = FALSE
    ),
    MSRate = fit_eeg_gaze_interaction_trials(
      dat, "MSRate", exploratory = FALSE
    )
  )

  interaction_exploratory <- list(
    GazeDeviationBL = fit_eeg_gaze_interaction_trials(
      dat, "GazeDeviationBL", exploratory = TRUE
    ),
    MSRateBL = fit_eeg_gaze_interaction_trials(
      dat, "MSRateBL", exploratory = TRUE
    )
  )

  results <- list(
    task = task,
    task_config = task_config,
    load_models = load_models,
    interaction_confirmatory = interaction_confirmatory,
    interaction_exploratory = interaction_exploratory,
    fitted_at = Sys.time(),
    level = "trials",
    random_effects = re_label_trials(dat),
    iqr_k = IQR_K_TRIALS
  )

  for (nm in names(load_models)) {
    res <- load_models[[nm]]
    if (!is.null(res)) {
      message(sprintf(
        "  %s ~ Load: family=%s, Anova type II | %s",
        res$dv, res$family, res$model_label
      ))
    }
  }

  for (block in list(interaction_confirmatory, interaction_exploratory)) {
    for (nm in names(block)) {
      res <- block[[nm]]
      if (!is.null(res)) {
        message(sprintf(
          "  %s ~ %s * Load: interaction_kept=%s, Anova type %d, LRT p=%.4f",
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
