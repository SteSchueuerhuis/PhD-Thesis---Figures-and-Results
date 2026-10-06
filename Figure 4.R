# Figure 4
# --------
# This script reproduces Figure 4 of the PhD thesis.

# install.packages("ggplot2")
library(ggplot2)


# Parameter settings ------------------------------------------------------

I1 <- 0.5
I2 <- 1

I1_tilde_star <- seq(0.1, 0.4, 0.1)
gamma <- seq(0.5, 0.9, 0.01)

f1_star <- 0
d1_star <- 1.96


# Create parameter grid --------------------------------------------------

df <- expand.grid(
  I1_tilde_star = I1_tilde_star,
  gamma = gamma
)

df$I1 <- I1
df$I2 <- I2
df$I1_star <- df$I1_tilde_star + df$I1
df$f1_star <- f1_star
df$d1_star <- d1_star


# Compute boundaries -----------------------------------------------------

df$f1 <- (
  df$f1_star * sqrt(df$I1_star / df$I1_tilde_star) -
    qnorm(df$gamma)
) * sqrt(df$I1_tilde_star / df$I1)

df$e1 <- (
  df$d1_star * sqrt(df$I1_star / df$I1_tilde_star) -
    qnorm(1 - df$gamma)
) * (
  sqrt(df$I1 / df$I1_tilde_star) +
    sqrt(df$I1_tilde_star / df$I1)
)^(-1)


# Create Figure 4 --------------------------------------------------------

ggplot(df, aes(x = gamma)) +
  geom_ribbon(
    aes(ymin = f1, ymax = e1),
    fill = "grey70",
    alpha = 0.25
  ) +
  geom_ribbon(
    aes(ymin = e1, ymax = 3),
    fill = "darkgreen",
    alpha = 0.15
  ) +
  geom_ribbon(
    aes(ymin = -3, ymax = f1),
    fill = "darkred",
    alpha = 0.15
  ) +
  geom_line(
    aes(y = f1, color = "f1"),
    linewidth = 1
  ) +
  geom_line(
    aes(y = e1, color = "e1"),
    linewidth = 1
  ) +
  facet_wrap(
    ~ I1_star,
    labeller = as_labeller(
      function(x) paste0("I[1]^'*'/I[2] == ", x),
      label_parsed
    )
  ) +
  scale_color_manual(
    values = c(
      "f1" = "darkred",
      "e1" = "darkgreen"
    ),
    labels = c(
      "f1" = expression(f[1]),
      "e1" = expression(e[1])
    )
  ) +
  scale_y_continuous(
    limits = c(-3, 3),
    breaks = seq(-3, 3, 1)
  ) +
  labs(
    x = expression(gamma),
    y = NULL,
    color = NULL
  ) +
  annotate(
    "text",
    x = 0.60,
    y = 2.60,
    label = "Pause recruitment",
    color = "grey20"
  ) +
  annotate(
    "text",
    x = 0.65,
    y = -1.80,
    label = "Pause recruitment",
    color = "grey20"
  ) +
  annotate(
    "text",
    x = 0.70,
    y = 1.00,
    label = "Continue recruitment",
    color = "grey20"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(
      size = 12,
      margin = margin(t = 1, b = 1)
    ),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )