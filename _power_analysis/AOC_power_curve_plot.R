# Shared ggplot styling for AOC SIMR power curves.
save_power_curve_plot <- function(power_curve_df, outfile) {
  out_dir <- dirname(outfile)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  base_font_size <- 15
  axis_font_size <- base_font_size * 1.5

  p <- ggplot(power_curve_df, aes(x = Subjects, y = Power)) +
    geom_point(size = 2, color = "blue") +
    geom_line(color = "blue", linetype = "dotted") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 3, color = "blue", alpha = 0.5) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "grey", linewidth = 0.5) +
    scale_x_continuous(
      breaks = c(0, 25, 50, 75, 100, 125)
    ) +
    scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.50, 0.75, 0.90, 1)
    ) +
    labs(
      x = "Subjects",
      y = "Power"
    ) +
    theme_minimal(base_size = base_font_size) +
    theme(
      axis.text = element_text(size = axis_font_size),
      axis.title = element_text(size = axis_font_size),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      axis.line = element_line(colour = "black", linewidth = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(0.2, "cm"),
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.x = element_line(colour = "black", linewidth = 0.5, lineend = "square"),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.5, lineend = "square"),
      axis.ticks.margin = unit(0.3, "cm"),
      axis.text.x = element_text(margin = margin(t = 10)),
      axis.text.y = element_text(margin = margin(r = 10))
    )

  ggsave(
    filename = outfile,
    plot = p,
    width = 1800 / 300,
    height = 1400 / 300,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}
