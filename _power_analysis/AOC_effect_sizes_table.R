# AOC Effect Sizes Table
# Literature effect sizes used for AOC power analysis (Sternberg, N-back, gaze on alpha).
#
# Key output (figures/power_analysis/):
#   AOC_effect_sizes_table.docx
#
# Run: Rscript AOC_effect_sizes_table.R

suppressPackageStartupMessages({
  library(flextable)
  library(officer)
})

out_dir <- "/Volumes/g_psyplafor_methlab$/Students/Arne/AOC/figures/power_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

TABLE_FONT <- "Arial"
TABLE_FONT_SIZE <- 10
EM_DASH <- "\u2014"

fmt_effect_size <- function(x) {
  ifelse(is.na(x), EM_DASH, format(x, trim = TRUE, scientific = FALSE))
}

# %% Literature data -----------------------------------------------------------

sternberg_alpha <- data.frame(
  reference = c(
    "Bashivan et al. (2014)",
    "Hu et al. (2019)",
    "Jensen (2002)",
    "Michels et al. (2010)",
    "Proskovec et al. (2019)",
    "Rominger et al. (2019)",
    "Scheeringa et al. (2009)",
    "Tuladhar et al. (2007)",
    "Wianda & Ross (2019)"
  ),
  n = c(15, 20, 10, 16, 26, 52, 20, 5, 25),
  conditions = c(
    "2, 4, 6, 8",
    "1, 3, 5",
    "2, 4, 6",
    "2, 5",
    "4, 6",
    "4",
    "0, 3, 5, 7",
    "1, 2, 3, 4",
    "1, 2, 3, 4, 5"
  ),
  stimuli = c(
    "letters",
    "digits",
    "letters",
    "letters",
    "letters",
    "shapes (mod)",
    "letters",
    "faces (mod)",
    "letters (mod)"
  ),
  effect_size = c(NA, 0.19, NA, NA, NA, NA, NA, NA, NA),
  stringsAsFactors = FALSE
)

nback_alpha <- data.frame(
  reference = c(
    "Brouwer et al. (2012)",
    "Chen & Huang (2016)",
    "Gevins & Smith (2000)",
    "Haegens et al. (2014)",
    "Scharinger et al. (2017)",
    "Pergher et al. (2019)",
    "Pesonen et al. (2007)"
  ),
  n = c(35, 18, 80, 51, 20, 20, 36),
  conditions = c(
    "0, 1, 2",
    "1, 2",
    "1, 2",
    "0, 1, 2",
    "1, 2, 3, 4",
    "0, 1, 2, 3",
    "0, 1, 2, 3"
  ),
  stimuli = c(
    "letters",
    "symbols",
    "letters",
    "letters",
    "digits",
    "objects",
    "letters"
  ),
  effect_size = c(NA, 1.74, 1.38, NA, 1.8, NA, NA),
  stringsAsFactors = FALSE
)

gaze_alpha <- data.frame(
  reference = c(
    "Popov et al. (2021a)",
    "Popov et al. (2021b)"
  ),
  n = c(143, 122),
  task = c(
    "prosaccade task\n(Plomecka et al., 2020)",
    "delayed matching to sample task\n(Luck & Vogel, 1997)"
  ),
  effect_size = c(1.5, 0.4),
  stringsAsFactors = FALSE
)

# %% Table formatting ----------------------------------------------------------

style_effect_size_ft <- function(ft, col_widths = NULL) {
  ft <- theme_booktabs(ft)
  rule_border <- fp_border(color = "#666666", width = 0.75)
  ft <- hline_top(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "header")
  ft <- hline_bottom(ft, border = rule_border, part = "body")
  ft <- font(ft, fontname = TABLE_FONT, part = "all")
  ft <- fontsize(ft, size = TABLE_FONT_SIZE, part = "all")
  ft <- bold(ft, part = "header")
  ft <- padding(ft, padding = 3, part = "all")
  ft <- align(ft, align = "left", j = 1, part = "all")
  if (ncol(ft$body$dataset) > 1) {
    ft <- align(ft, align = "center", j = 2:ncol(ft$body$dataset), part = "all")
  }
  ft <- valign(ft, valign = "center", part = "all")
  ft <- set_table_properties(ft, layout = "fixed")
  if (!is.null(col_widths) && length(col_widths) == ncol(ft$body$dataset)) {
    ft <- width(ft, j = seq_along(col_widths), width = col_widths)
  } else {
    ft <- autofit(ft)
  }
  ft
}

make_wm_load_flextable <- function(df) {
  body <- data.frame(
    Reference = df$reference,
    N = df$n,
    Conditions = df$conditions,
    Stimuli = df$stimuli,
    `Effect Size` = fmt_effect_size(df$effect_size),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ft <- flextable(body)
  style_effect_size_ft(
    ft,
    col_widths = c(2.2, 0.45, 1.0, 1.0, 0.75)
  )
}

make_gaze_flextable <- function(df) {
  body <- data.frame(
    Reference = df$reference,
    N = df$n,
    Task = df$task,
    `Effect Size` = fmt_effect_size(df$effect_size),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  ft <- flextable(body)
  style_effect_size_ft(
    ft,
    col_widths = c(2.0, 0.45, 2.4, 0.75)
  )
}

add_section <- function(doc, title, ft) {
  doc <- body_add_fpar(doc, fpar(
    ftext(title, prop = fp_text(bold = TRUE, font.family = TABLE_FONT, font.size = TABLE_FONT_SIZE))
  ))
  doc <- body_add_flextable(doc, ft)
  body_add_par(doc, "")
}

# %% Export --------------------------------------------------------------------

doc <- read_docx()
doc <- add_section(
  doc,
  "Studies reporting on the effect of WM load in the Sternberg task on alpha power",
  make_wm_load_flextable(sternberg_alpha)
)
doc <- add_section(
  doc,
  "Studies reporting on the effect of WM load in the N-back task on alpha power",
  make_wm_load_flextable(nback_alpha)
)
doc <- add_section(
  doc,
  "Studies reporting on the effect of gaze on alpha power",
  make_gaze_flextable(gaze_alpha)
)

docx_path <- file.path(out_dir, "AOC_effect_sizes_table.docx")
print(doc, target = docx_path)
message("Saved: ", docx_path)
