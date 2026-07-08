# Data loading and preprocessing for AOC trial-level GLMMs.
# Input: AOC_merged_data_*_trials.csv (separate from condition-level merged_data).

AOC_BASE_DIR <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC"
AOC_STATS_DIR_TRIALS <- file.path(AOC_BASE_DIR, "data", "stats", "trials")

# Trial-level EEG is ERSD only; map task-specific windows to canonical DV names.
ERSD_BY_TASK_TRIALS <- list(
  nback = "ERSD_full",
  sternberg = "ERSD_late"
)

# Unbaselined gaze / MS for confirmatory models (trial early/late/full variants).
GAZE_RAW_SOURCE_BY_TASK_TRIALS <- list(
  nback = list(
    GazeDeviation = "GazeDeviationFull",
    MSRate = "MSRateFull"
  ),
  sternberg = list(
    GazeDeviation = "GazeDeviationLate",
    MSRate = "MSRateLate"
  )
)

# Baselined gaze / MS for exploratory models.
GAZE_BL_SOURCE_BY_TASK_TRIALS <- list(
  nback = list(
    GazeDeviationBL = "GazeDeviationFullBL",
    MSRateBL = "MSRateFullBL"
  ),
  sternberg = list(
    GazeDeviationBL = "GazeDeviationLateBL",
    MSRateBL = "MSRateLateBL"
  )
)

GLMM_VARIABLES_TRIALS <- c(
  "Accuracy", "ReactionTime",
  "GazeDeviation", "MSRate",
  "GazeDeviationBL", "MSRateBL",
  "ERSD"
)

NBACK_TASK_CONFIG_TRIALS <- list(
  name = "nback",
  input_csv = file.path(AOC_BASE_DIR, "data", "features", "AOC_merged_data_nback_trials.csv"),
  cond_to_label_numeric = list(
    c(`1` = "1-back", `2` = "2-back", `3` = "3-back")
  ),
  categories = c("1-back", "2-back", "3-back"),
  comparisons = list(
    c("1-back", "2-back"),
    c("1-back", "3-back"),
    c("2-back", "3-back")
  )
)

STERNBERG_TASK_CONFIG_TRIALS <- list(
  name = "sternberg",
  input_csv = file.path(AOC_BASE_DIR, "data", "features", "AOC_merged_data_sternberg_trials.csv"),
  cond_to_label_numeric = list(
    c(`2` = "WM load 2", `4` = "WM load 4", `6` = "WM load 6")
  ),
  categories = c("WM load 2", "WM load 4", "WM load 6"),
  comparisons = list(
    c("WM load 2", "WM load 4"),
    c("WM load 2", "WM load 6"),
    c("WM load 4", "WM load 6")
  )
)

iqr_mask_series_trials <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || sum(!is.na(x)) < 3) return(x)
  q1 <- stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE)
  q3 <- stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE)
  iqr <- q3 - q1
  if (!is.finite(iqr) || iqr == 0) return(x)
  lower <- q1 - 1.5 * iqr
  upper <- q3 + 1.5 * iqr
  x[!is.na(x) & (x < lower | x > upper)] <- NA_real_
  x
}

iqr_outlier_filter_trials <- function(df, variables, by = "Load") {
  out <- df
  for (v in variables) {
    if (!v %in% names(out)) next
    out[[v]] <- as.numeric(out[[v]])
    out[[v]] <- stats::ave(out[[v]], out[[by]], FUN = iqr_mask_series_trials)
  }
  out
}

enforce_task_window_columns_trials <- function(df, task_name) {
  ersd_col <- ERSD_BY_TASK_TRIALS[[task_name]]
  if (!is.null(ersd_col) && ersd_col %in% names(df)) {
    df$ERSD <- as.numeric(df[[ersd_col]])
  }

  gaze_raw <- GAZE_RAW_SOURCE_BY_TASK_TRIALS[[task_name]]
  if (!is.null(gaze_raw)) {
    for (canonical in names(gaze_raw)) {
      src <- gaze_raw[[canonical]]
      if (src %in% names(df)) {
        df[[canonical]] <- as.numeric(df[[src]])
      }
    }
  }

  gaze_bl <- GAZE_BL_SOURCE_BY_TASK_TRIALS[[task_name]]
  if (!is.null(gaze_bl)) {
    for (canonical in names(gaze_bl)) {
      src <- gaze_bl[[canonical]]
      if (src %in% names(df)) {
        df[[canonical]] <- as.numeric(df[[src]])
      }
    }
  }
  df
}

label_load_trials <- function(dat, task_config) {
  cond <- dat$Condition
  categories <- task_config$categories
  maps <- task_config$cond_to_label_numeric

  if (is.numeric(cond)) {
    uniq <- sort(unique(cond[!is.na(cond)]))
    applied <- NULL
    for (cand in maps) {
      if (all(uniq %in% as.numeric(names(cand)))) {
        applied <- cand
        break
      }
    }
    if (is.null(applied)) {
      applied <- stats::setNames(
        categories[seq_along(uniq)],
        as.character(uniq)
      )
    }
    lab <- applied[as.character(cond)]
    dat$Load <- factor(lab, levels = categories, ordered = TRUE)
  } else {
    lab <- trimws(gsub("\\s+", " ", as.character(cond)))
    dat$Load <- factor(lab, levels = categories, ordered = TRUE)
  }

  dat$Subject <- as.factor(dat$ID)
  dat$ID <- dat$Subject
  dat
}

load_task_data_trials <- function(task_config) {
  if (!file.exists(task_config$input_csv)) {
    stop("Missing trial-level merged CSV: ", task_config$input_csv)
  }
  dat <- read.csv(task_config$input_csv, stringsAsFactors = FALSE)
  dat <- enforce_task_window_columns_trials(dat, task_config$name)
  if ("Accuracy" %in% names(dat)) {
    dat$Accuracy[dat$Accuracy > 100] <- NA_real_
  }
  dat <- label_load_trials(dat, task_config)
  vars_present <- intersect(GLMM_VARIABLES_TRIALS, names(dat))
  iqr_outlier_filter_trials(dat, vars_present, by = "Load")
}
