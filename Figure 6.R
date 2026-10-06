library(ggplot2)

# Figure 6: Illustration of pipeline data in the DISC trial ---------------
#
# This figure provides a simplified illustration of how pipeline data may
# emerge in the DISC trial due to the delay between patient recruitment and
# completion of the 12-month follow-up. The recruitment trajectory is
# illustrative and does not represent the actual recruitment history of the
# trial.


# Parameters --------------------------------------------------------------

n_total <- 710
recruitment_months <- 39
follow_up_months <- 12
pause_month <- 33


# Recruitment and follow-up ----------------------------------------------

# Time from start of recruitment until completion of final follow-up
time <- seq(
  0,
  recruitment_months + follow_up_months,
  by = 0.1
)

# Illustrative recruitment curve with slow initial recruitment
# followed by acceleration
recruitment_curve <- function(time) {
  pmin(
    n_total,
    n_total * (time / recruitment_months)^2
  )
}

# Cumulative number of recruited patients
recruited <- ifelse(
  time <= recruitment_months,
  recruitment_curve(time),
  n_total
)

# Cumulative number of patients with completed 12-month follow-up
completed_fu <- ifelse(
  time <= follow_up_months,
  0,
  recruitment_curve(
    pmin(time - follow_up_months, recruitment_months)
  )
)

# Patients recruited but without completed 12-month follow-up
pipeline <- recruited - completed_fu

plot_data <- data.frame(
  month = time,
  recruited = recruited,
  completed_fu = completed_fu,
  pipeline = pipeline
)


# Figure 6 ----------------------------------------------------------------

fig6 <- ggplot(plot_data, aes(x = month)) +
  
  # Patients in the follow-up pipeline
  geom_ribbon(
    aes(
      ymin = completed_fu,
      ymax = recruited,
      fill = "Pipeline data"
    ),
    alpha = 0.5
  ) +
  
  # Cumulative recruitment
  geom_line(
    aes(
      y = recruited,
      color = "Recruited patients",
      linetype = "Recruited patients"
    ),
    linewidth = 1.3
  ) +
  
  # Completed follow-up
  geom_line(
    aes(
      y = completed_fu,
      color = "Completed follow-up",
      linetype = "Completed follow-up"
    ),
    linewidth = 1.3
  ) +
  
  # Illustrative trial pause
  geom_vline(
    xintercept = pause_month,
    aes(color = "Trial pause"),
    linetype = "dotted",
    linewidth = 1
  ) +
  
  annotate(
    "text",
    x = pause_month + 1,
    y = 0.12 * n_total,
    label = "Trial paused",
    hjust = 0,
    size = 4
  ) +
  
  labs(
    x = "Months since recruitment start",
    y = "Number of patients"
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Recruited patients" = "black",
      "Completed follow-up" = "black",
      "Trial pause" = "red3"
    )
  ) +
  
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Recruited patients" = "solid",
      "Completed follow-up" = "dashed"
    )
  ) +
  
  scale_fill_manual(
    name = NULL,
    values = c(
      "Pipeline data" = "grey70"
    )
  ) +
  
  guides(
    fill = guide_legend(order = 1, nrow = 1),
    color = guide_legend(order = 2, nrow = 1),
    linetype = guide_legend(order = 2, nrow = 1)
  ) +
  
  theme_bw(base_size = 17) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = grid::unit(2.2, "cm"),
    legend.spacing.x = grid::unit(0.6, "cm")
  )
