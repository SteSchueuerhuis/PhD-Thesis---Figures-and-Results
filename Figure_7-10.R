# Figures 7-10
# ------------
# This script reproduces Figures 7-10 of the PhD thesis.
#
# The script:
#   1. defines helper functions for the binomial simulations,
#   2. constructs the group-sequential design,
#   3. determines the blinded unblinding regions,
#   4. evaluates the operating characteristics under three scenarios, and
#   5. creates Figures 7-10.

library(rpact)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)

# Helper functions --------------------------------------------------------

# Simulate binary outcomes for a two-arm trial with treatment (T) and
# control (C) response probabilities piT and piC. The interim analysis is
# performed after the prespecified fraction of observations is available
# in each treatment group.
#
# For both the interim and final data, the treatment effect is evaluated
# using the pooled two-sample z-statistic for two binomial proportions.
# The function returns the simulated data, the interim subset, the
# corresponding test statistics, and the sample sizes.
simulateBinomial <- function(
    nT,
    nC,
    piT,
    piC,
    interim_frac = 0.5,
    seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Determine group-specific sample sizes available at interim
  nT_interim <- floor(nT * interim_frac)
  nC_interim <- floor(nC * interim_frac)
  
  # Simulate binary outcomes in the treatment and control groups
  data <- data.frame(
    id = seq_len(nT + nC),
    group = c(rep("T", nT), rep("C", nC)),
    outcome = c(
      rbinom(nT, size = 1, prob = piT),
      rbinom(nC, size = 1, prob = piC)
    )
  )
  
  # Identify observations available at the interim analysis
  data$interim_available <- with(
    data,
    (group == "T" & ave(id, group, FUN = seq_along) <= nT_interim) |
      (group == "C" & ave(id, group, FUN = seq_along) <= nC_interim)
  )
  
  # Compute the pooled two-sample z-statistic for two binomial proportions.
  # The statistic is defined such that positive values favor treatment
  # when a lower response probability represents a beneficial outcome.
  binom_z <- function(dat) {
    
    yT <- sum(dat$outcome[dat$group == "T"])
    yC <- sum(dat$outcome[dat$group == "C"])
    
    nT <- sum(dat$group == "T")
    nC <- sum(dat$group == "C")
    
    pT <- yT / nT
    pC <- yC / nC
    
    p_pool <- (yT + yC) / (nT + nC)
    
    se <- sqrt(
      p_pool * (1 - p_pool) * (1 / nT + 1 / nC)
    )
    
    # If all outcomes are either 1 or 0, the SE is zero. In those cases, we manually set
    # z = 0, reflecting no observed treatment difference.
    if (se == 0) {
      z <- 0
    } else {
      z <- (pC - pT) / se
    }
    
    pval <- 2 * pnorm(-abs(z))
    
    list(
      z = z,
      p.value = pval,
      pT = pT,
      pC = pC,
      nT = nT,
      nC = nC
    )
  }
  
  # Apply the test to the interim and complete trial data
  interim_data <- subset(data, interim_available)
  
  interim_test <- binom_z(interim_data)
  final_test <- binom_z(data)
  
  list(
    data = data,
    interim_data = interim_data,
    interim_test = interim_test,
    final_test = final_test,
    sample_sizes = list(
      treatment = nT,
      control = nC,
      treatment_interim = nT_interim,
      control_interim = nC_interim
    )
  )
}


# Compute the blinded estimate of interim power as a function of the
# observed pooled response probability. The prespecified treatment effect
# (piT - piC) is combined with the blinded pooled estimate to estimate
# treatment-specific response probabilities under the assumed alternative.
#
# The resulting conditional power estimate is used below to determine
# whether the trial is unblinded at the interim analysis.
computeInterimPower <- function(
    piT,
    piC,
    alpha1,
    pi_pooled,
    kappa,
    nT1
) {
  
  # Reconstruct treatment-specific probabilities from the blinded pooled
  # response probability while retaining the prespecified treatment effect
  piTstar <- pmin(
    pi_pooled - (piT - piC) / (1 + kappa),
    1
  )
  
  piCstar <- pmax(
    pi_pooled + kappa * (piT - piC) / (1 + kappa),
    0
  )
  
  # Blinded estimate of interim power
  numerator <- sqrt(nT1) * (piTstar - piCstar) -
    qnorm(1 - alpha1) *
    sqrt((1 + kappa) * pi_pooled * (1 - pi_pooled))
  
  denominator <- sqrt(
    piTstar * (1 - piTstar) +
      kappa * piCstar * (1 - piCstar)
  )
  
  pnorm(numerator / denominator)
}

# Evaluate the operating characteristics of the two hybrid designs,
# the single-stage design, and the conventional group-sequential design
# across a grid of treatment and control response probabilities.
runSimulation <- function(
    piT_grid = seq(0, 1, by = 0.05),
    piC_grid = seq(0, 1, by = 0.05),
    nsim = 100000,
    planning,
    d,
    pLU,
    pUU,
    pLU2,
    pUU2,
    interim_frac = 0.5,
    seed = 2905
) {
  
  if (length(piT_grid) != length(piC_grid)) {print("Number of values for treatment and control group need to align")}
  
  nC <- ceiling(planning$numberOfSubjects2[2, 1])
  nT <- ceiling(planning$numberOfSubjects2[2, 1])
  
  crit <- d$criticalValues
  crit_single <- qnorm(1 - 0.025)
  fut <- d$futilityBounds
  
  out <- data.frame()
  
  for (k in seq_along(piT_grid)) {
    
    rej1 <- rej2 <- rej3 <- rej4 <-  numeric(nsim)   
    unb1 <- unb2 <- unb3 <- unb4 <-  numeric(nsim)  
    eN1 <- eN2 <- eN3 <- eN4 <- numeric(nsim)
    rejInt1 <- rejInt2 <- rejInt3 <- rejInt4 <- numeric(nsim)
    fut1 <- fut2 <- fut3 <- fut4 <- numeric(nsim)
    
    for (i in seq_len(nsim)) {
      
      res <- simulateBinomial(
        piC = piC_grid[k],
        piT = piT_grid[k],
        nC = nC,
        nT = nT,
        interim_frac = interim_frac,
        seed = seed + i
      )
      
      pi_pooled <- mean(res$interim_data$outcome)
      
      z1 <- res$interim_test$z
      z2 <- res$final_test$z
      
      
      # Design 1
      unblind1 <- !((pLU < pi_pooled) & (pi_pooled < pUU))
      unb1[i] <- unblind1
      rejInt1[i] <- ifelse(unblind1 & z1 > crit[1], 1, 0)
      rej1[i] <- ifelse(
        unblind1,
        z1 > crit[1] | ((z1 < crit[1] & z1 > fut[1]) & z2 > crit[2]),
        z2 > crit_single
      )
      fut1[i] <- ifelse(unblind1 & z1 < fut[1], 1, 0)
      eN1[i] <- ifelse(rejInt1[i] == 1 | fut1[i] == 1, (nC + nT)*interim_frac, nC + nT)
      
      # Design 2
      unblind2 <- !((pLU2 < pi_pooled) & (pi_pooled < pUU2))
      unb2[i] <- unblind2
      rejInt2[i] <- ifelse(unblind2 & z1 > crit[1], 1, 0)
      rej2[i] <- ifelse(
        unblind2,
        z1 > crit[1] | ((z1 < crit[1] & z1 > fut[1]) & z2 > crit[2]),
        z2 > crit_single
      )
      fut2[i] <- ifelse(unblind2 & z1 < fut[1], 1, 0)
      eN2[i] <- ifelse(rejInt2[i] == 1 | fut2[i] == 1, (nC + nT)*interim_frac, nC + nT)
      
      # Design 3: Single stage
      unb3[i] <- FALSE
      rej3[i] <- z2 > crit_single
      rejInt3[i] <- 0
      fut3[i] <- 0
      eN3[i] <- ifelse(rejInt3[i] == 1 | fut3[i] == 1, (nC + nT)*interim_frac, nC + nT)
      
      # Design 4: Group-sequential design
      unb4[i] <- TRUE
      rej4[i] <- z1 > crit[1] | ((z1 < crit[1] & z1 > fut[1]) & z2 > crit[2])
      rejInt4[i] <- z1 > crit[1]
      fut4[i] <- z1 < fut[1]
      eN4[i] <- ifelse(rejInt4[i] == 1 | fut4[i] == 1, (nC + nT)*interim_frac, nC + nT)
      
    }
    
    out <- rbind(
      out,
      data.frame(
        piC = piC_grid[k],
        piT = piT_grid[k],
        seed_init = seed,
        nsim = nsim,
        nC = nC,
        nT = nT,
        rejection_rate_design1 = mean(rej1),
        rejection_rate_design2 = mean(rej2),
        rejection_rate_design3 = mean(rej3),
        rejection_rate_design4 = mean(rej4),
        rejection_rate_design1_interim = mean(rejInt1),
        rejection_rate_design2_interim = mean(rejInt2),
        rejection_rate_design3_interim = mean(rejInt3),
        rejection_rate_design4_interim = mean(rejInt4),
        fut_design1 = mean(fut1),
        fut_design2 = mean(fut2),
        fut_design3 = mean(fut3),
        fut_design4 = mean(fut4),
        unblinding_rate_design1 = mean(unb1),
        unblinding_rate_design2 = mean(unb2),
        unblinding_rate_design3 = mean(unb3),
        unblinding_rate_design4 = mean(unb4),
        eN1 = mean(eN1),
        eN2 = mean(eN2),
        eN3 = mean(eN3),
        eN4 = mean(eN4)
      )
    )
  }
  
  out
}

# Group-sequential design -------------------------------------------------

# Reconstruct the BACLOREA trial as the two-stage group-sequential design
# considered in Section 3.3.1 of the thesis. The assumed response
# probabilities are 0.27 in the treatment group and 0.42 in the control
# group.
piT <- 0.27
piC <- 0.42

d <- getDesignGroupSequential(
  kMax = 2,
  typeOfDesign = "asKD",
  futilityBounds = 0,
  bindingFutility = TRUE,
  gammaA = 2
)

planning <- d |>
  getSampleSizeRates(
    pi1 = piT,
    pi2 = piC
  )


# Blinded interim power ---------------------------------------------------

# Evaluate the blinded estimate of interim power over the full range of
# possible pooled response probabilities.
alpha1 <- d$alphaSpent[1]
pi_pooled <- seq(0.01, 0.99, 0.001)
kappa <- 1
nT1 <- 81

pwr <- computeInterimPower(
  piT = piT,
  piC = piC,
  alpha1 = alpha1,
  pi_pooled = pi_pooled,
  kappa = kappa,
  nT1 = nT1
)

# Define the two hybrid designs -------------------------------------------

# Hybrid Design 1 uses a symmetric interval of +/- 0.10 around the
# prespecified interim power. Hybrid Design 2 instead uses fixed lower
# and upper power thresholds of 0.10 and 0.70.
df <- data.frame(
  pi_pooled = pi_pooled,
  pwr = pwr,
  pwrU = planning$rejectPerStage[1] + 0.10,
  pwrL = planning$rejectPerStage[1] - 0.10,
  pwrU2 = 0.70,
  pwrL2 = 0.10
)

# Determine unblinding regions -------------------------------------------

# Translate the power thresholds of the two hybrid designs to the pooled
# response-probability scale. For each threshold, the grid value producing
# the closest blinded power is selected separately below and above 0.5.
#
# pLU / pUU: lower and upper boundaries corresponding to the upper
#             power threshold of Hybrid Design 1
# pLL / pUL: lower and upper boundaries corresponding to the lower
#             power threshold of Hybrid Design 1
# The quantities ending in "2" are the corresponding boundaries for
# Hybrid Design 2.

pLU <- df %>%
  filter(pi_pooled < 0.5) %>%
  mutate(diff = abs(pwr - pwrU)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()

pUU <- df %>%
  filter(pi_pooled > 0.5) %>%
  mutate(diff = abs(pwr - pwrU)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()

pLL <- df %>%
  filter(pi_pooled < 0.5) %>%
  mutate(diff = abs(pwr - pwrL)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()

pUL <- df %>%
  filter(pi_pooled > 0.5) %>%
  mutate(diff = abs(pwr - pwrL)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()

pLU2 <- df %>%
  filter(pi_pooled < 0.5) %>%
  mutate(diff = abs(pwr - pwrU2)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()

pUU2 <- df %>%
  filter(pi_pooled > 0.5) %>%
  mutate(diff = abs(pwr - pwrU2)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()

pLL2 <- df %>%
  filter(pi_pooled < 0.5) %>%
  mutate(diff = abs(pwr - pwrL2)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()

pUL2 <- df %>%
  filter(pi_pooled > 0.5) %>%
  mutate(diff = abs(pwr - pwrL2)) %>%
  slice_min(diff, n = 1) %>%
  select(pi_pooled) %>%
  pull()


# Figure 7 ----------------------------------------------------------------

# Figure 7 illustrates the blinded interim-power rule for the two hybrid
# designs. The green shaded regions correspond to values of the blinded
# pooled response probability for which the trial is unblinded. 
# In the red shaded region, recruitment continues without unblinding.

# Figure 7a: Hybrid Design 1
fig7a <- ggplot(df, aes(x = pi_pooled, y = pwr)) +
  geom_line(
    color = "#2C7BB6",
    linewidth = 1.2
  ) +
  geom_point(
    data = data.frame(
      x = mean(c(piT, piC)),
      y = planning$rejectPerStage[1]
    ),
    aes(x = x, y = y),
    size = 5,
    color = "#D7191C"
  ) +
  geom_point(
    data = data.frame(
      x = pLU,
      y = planning$rejectPerStage[1] + 0.1
    ),
    aes(x = x, y = y),
    size = 5,
    color = "#FDAE61"
  ) +
  geom_point(
    data = data.frame(
      x = pUU,
      y = planning$rejectPerStage[1] + 0.1
    ),
    aes(x = x, y = y),
    size = 5,
    color = "#FDAE61"
  ) +
  geom_hline(
    yintercept = planning$rejectPerStage[1],
    linetype = "dashed",
    color = "#D7191C",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = planning$rejectPerStage[1] - 0.1,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  geom_hline(
    yintercept = planning$rejectPerStage[1] + 0.1,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = mean(c(piT, piC)),
    linetype = "dashed",
    color = "#D7191C",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = pLU,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = pUU,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03)),
    breaks = c(
      0,
      planning$rejectPerStage[1] - 0.1,
      planning$rejectPerStage[1] + 0.1,
      0.6, 0.7, 0.8, 1
    ),
    labels = scales::label_number(accuracy = 0.01),
    guide = guide_axis(check.overlap = FALSE)
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03)),
    breaks = c(0, pLU, 0.5, pUU, 1),
    labels = scales::label_number(accuracy = 0.01),
    guide = guide_axis(check.overlap = FALSE)
  ) +
  labs(
    x = expression(hat(pi)[pooled]^(1)),
    y = "Blinded estimate of the interim power"
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = pLU,
    ymin = -Inf, ymax = Inf,
    fill = "#90EE90", alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = pUU, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "#90EE90", alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = pLU, xmax = pUU,
    ymin = -Inf, ymax = Inf,
    fill = "#FF474C", alpha = 0.15
  ) +
  annotate(
    "text",
    x = 0.05, y = 0.45,
    label = "'>'~0.41",
    parse = TRUE,
    size = 8
  ) +
  annotate(
    "text",
    x = 0.2, y = 0.75,
    label = "Unblinding Area",
    size = 8,
    angle = 90
  ) +
  annotate(
    "text",
    x = 0.8, y = 0.75,
    label = "Unblinding Area",
    size = 8,
    angle = 90
  ) +
  theme_bw(base_size = 25) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )


# Figure 7b: Hybrid Design 2
fig7b <- ggplot(df, aes(x = pi_pooled, y = pwr)) +
  geom_line(
    color = "#2C7BB6",
    linewidth = 1.2
  ) +
  geom_point(
    data = data.frame(
      x = mean(c(piT, piC)),
      y = planning$rejectPerStage[1]
    ),
    aes(x = x, y = y),
    size = 5,
    color = "#D7191C"
  ) +
  geom_point(
    data = data.frame(x = pLU2, y = 0.7),
    aes(x = x, y = y),
    size = 5,
    color = "#FDAE61"
  ) +
  geom_point(
    data = data.frame(x = pUU2, y = 0.7),
    aes(x = x, y = y),
    size = 5,
    color = "#FDAE61"
  ) +
  geom_hline(
    yintercept = planning$rejectPerStage[1],
    linetype = "dashed",
    color = "#D7191C",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0.1,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  geom_hline(
    yintercept = 0.7,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = mean(c(piT, piC)),
    linetype = "dashed",
    color = "#D7191C",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = pLU2,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = pUU2,
    linetype = "dotted",
    color = "#FDAE61",
    linewidth = 1.2
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03)),
    breaks = c(0, 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 1),
    labels = scales::label_number(accuracy = 0.01),
    guide = guide_axis(check.overlap = FALSE)
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.03)),
    breaks = c(0, pLU2, 0.25, 0.5, 0.75, pUU2, 1),
    labels = scales::label_number(accuracy = 0.01),
    guide = guide_axis(check.overlap = FALSE)
  ) +
  labs(
    x = expression(hat(pi)[pooled]^(1)),
    y = "Blinded estimate of the interim power"
  ) +
  annotate(
    "rect",
    xmin = -Inf, xmax = pLU2,
    ymin = -Inf, ymax = Inf,
    fill = "#90EE90", alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = pUU2, xmax = Inf,
    ymin = -Inf, ymax = Inf,
    fill = "#90EE90", alpha = 0.2
  ) +
  annotate(
    "rect",
    xmin = pLU2, xmax = pUU2,
    ymin = -Inf, ymax = Inf,
    fill = "#FF474C", alpha = 0.15
  ) +
  annotate(
    "text",
    x = 0.05, y = 0.74,
    label = "'>'~0.7",
    parse = TRUE,
    size = 8
  ) +
  annotate(
    "text",
    x = 0.09, y = 0.355,
    label = "Unblinding Area",
    size = 8,
    angle = 90
  ) +
  annotate(
    "text",
    x = 0.905, y = 0.355,
    label = "Unblinding Area",
    size = 8,
    angle = 90
  ) +
  theme_bw(base_size = 25) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )




# Figure 8: Null hypothesis -----------------------------------------------

sim_fig8 <- runSimulation(
  piT_grid = seq(0.05, 0.95, by = 0.05),
  piC_grid = seq(0.05, 0.95, by = 0.05),
  nsim = 10000,
  planning = planning,
  d = d,
  pLU = pLU,
  pUU = pUU,
  pLU2 = pLU2,
  pUU2 = pUU2,
  interim_frac = 0.5,
  seed = 1234
)

names(sim_fig8) <- gsub(
  pattern = "rejection_rate_design([1-4])_interim",
  replacement = "rejection_rate_interim_design\\1",
  x = names(sim_fig8)
)

sim_fig8_long <- sim_fig8 %>%
  rename(
    eN_design1 = eN1,
    eN_design2 = eN2,
    eN_design3 = eN3,
    eN_design4 = eN4
  ) %>%
  pivot_longer(
    cols = matches(
      "^(rejection_rate|rejection_rate_interim|unblinding_rate|eN|fut)_design[1-4]"
    ),
    names_to = c("performance_characteristic", "method"),
    names_pattern = "(.+)_design([1-4])",
    values_to = "value"
  ) %>%
  mutate(
    method = paste0("design", method)
  )

plot_fig8 <- sim_fig8_long %>%
  mutate(
    performance_characteristic = recode(
      performance_characteristic,
      rejection_rate = "Type-I Error Rate",
      rejection_rate_interim = "Interim Rejection Rate",
      unblinding_rate = "Unblinding Probability",
      eN = "Expected Sample Size",
      fut = "Futility Probability"
    ),
    method = recode(
      method,
      design1 = "Hybrid Design 1",
      design2 = "Hybrid Design 2",
      design3 = "Single-Stage Design",
      design4 = "Group-Sequential Design"
    )
  )

n_control <- ceiling(planning$numberOfSubjects2[2, 1])
n_treatment <- ceiling(planning$numberOfSubjects2[2, 1])

ref_lines_fig8 <- data.frame(
  performance_characteristic = c(
    "Type-I Error Rate",
    "Type-I Error Rate",
    "Type-I Error Rate",
    "Expected Sample Size",
    "Expected Sample Size"
  ),
  yintercept = c(
    0.025,
    0.025 + qnorm(1 - 0.025) *
      sqrt(0.025 * (1 - 0.025) / 1e4),
    0.025 - qnorm(1 - 0.025) *
      sqrt(0.025 * (1 - 0.025) / 1e4),
    n_control,
    n_control + n_treatment
  ),
  linetype = c(
    "dashed",
    "dotted",
    "dotted",
    "solid",
    "solid"
  )
)

fig8 <- plot_fig8 %>%
  filter(performance_characteristic != "Interim Rejection Rate") %>%
  ggplot(
    aes(
      x = piC,
      y = value,
      color = method
    )
  ) +
  geom_hline(
    data = ref_lines_fig8,
    aes(
      yintercept = yintercept,
      linetype = linetype
    ),
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ performance_characteristic,
    scales = "free_y",
    ncol = 2
  ) +
  facetted_pos_scales(
    y = list(
      performance_characteristic == "Type-I Error Rate" ~
        scale_y_continuous(limits = c(0, 0.05)),
      
      performance_characteristic == "Interim Rejection Rate" ~
        scale_y_continuous(limits = c(0, 0.025)),
      
      performance_characteristic == "Unblinding Probability" ~
        scale_y_continuous(limits = c(0, 1)),
      
      performance_characteristic == "Expected Sample Size" ~
        scale_y_continuous(limits = c(140, 350)),
      
      performance_characteristic == "Futility Probability" ~
        scale_y_continuous(limits = c(0, 1))
    )
  ) +
  labs(
    x = expression(
      "Common response probability " * pi[C] == pi[T]
    ),
    y = NULL,
    color = NULL,
    linetype = NULL
  ) +
  scale_x_continuous(
    breaks = seq(0.05, 0.95, by = 0.1),
    limits = c(0.05, 0.95)
  ) +
  scale_color_manual(
    values = c(
      "#0072B2",
      "#D55E00",
      "#009E73",
      "#CC79A7"
    )
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey40"
    ),
    strip.text = element_text(
      face = "bold",
      size = 15
    ),
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(
      size = 15,
      margin = margin(t = 10)
    ),
    axis.text = element_text(size = 15),
    plot.margin = margin(10, 10, 10, 10)
  )


# Figure 9: Fixed treatment effect ----------------------------------------

# Evaluate the operating characteristics while keeping the true treatment
# effect fixed at piC - piT = 0.15 and varying the underlying response rate.

delta <- 0.15

piT <- seq(0, 0.85, 0.05)
piC <- piT + delta

sim_fig9 <- runSimulation(
  piT_grid = piT,
  piC_grid = piC,
  nsim = 10000,
  planning = planning,
  d = d,
  pLU = pLU,
  pUU = pUU,
  pLU2 = pLU2,
  pUU2 = pUU2,
  interim_frac = 0.5,
  seed = 1234
)

names(sim_fig9) <- gsub(
  pattern = "rejection_rate_design([1-4])_interim",
  replacement = "rejection_rate_interim_design\\1",
  x = names(sim_fig9)
)

sim_fig9_long <- sim_fig9 %>%
  rename(
    eN_design1 = eN1,
    eN_design2 = eN2,
    eN_design3 = eN3,
    eN_design4 = eN4
  ) %>%
  pivot_longer(
    cols = matches(
      "^(rejection_rate|rejection_rate_interim|unblinding_rate|eN|fut)_design[1-4]"
    ),
    names_to = c("performance_characteristic", "method"),
    names_pattern = "(.+)_design([1-4])",
    values_to = "value"
  ) %>%
  mutate(
    method = paste0("design", method)
  )

plot_fig9 <- sim_fig9_long %>%
  mutate(
    performance_characteristic = recode(
      performance_characteristic,
      rejection_rate = "Power",
      rejection_rate_interim = "Interim Rejection Rate",
      unblinding_rate = "Unblinding Probability",
      eN = "Expected Sample Size",
      fut = "Futility Probability"
    ),
    method = recode(
      method,
      design1 = "Hybrid Design 1",
      design2 = "Hybrid Design 2",
      design3 = "Single-Stage Design",
      design4 = "Group-Sequential Design"
    )
  )

ref_lines_fig9 <- data.frame(
  performance_characteristic = c(
    "Expected Sample Size",
    "Expected Sample Size"
  ),
  yintercept = c(
    n_control,
    n_control + n_treatment
  ),
  linetype = c(
    "dotted",
    "dotted"
  )
)

fig9 <- plot_fig9 %>%
  filter(performance_characteristic != "Futility Probability") %>%
  ggplot(
    aes(
      x = piT,
      y = value,
      color = method
    )
  ) +
  geom_hline(
    data = ref_lines_fig9,
    aes(yintercept = yintercept),
    linetype = "dashed",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ performance_characteristic,
    scales = "free_y",
    ncol = 2
  ) +
  facetted_pos_scales(
    y = list(
      performance_characteristic == "Power" ~
        scale_y_continuous(limits = c(0.7, 1)),
      
      performance_characteristic == "Interim Rejection Rate" ~
        scale_y_continuous(limits = c(0, 1)),
      
      performance_characteristic == "Unblinding Probability" ~
        scale_y_continuous(limits = c(0, 1)),
      
      performance_characteristic == "Expected Sample Size" ~
        scale_y_continuous(limits = c(140, 350)),
      
      performance_characteristic == "Futility Probability" ~
        scale_y_continuous(limits = c(0, 1))
    )
  ) +
  labs(
    x = expression("Response probability " * pi[T]),
    y = NULL,
    color = NULL,
    linetype = NULL
  ) +
  scale_x_continuous(
    breaks = seq(0, 0.85, by = 0.1),
    limits = c(0, 0.85)
  ) +
  scale_color_manual(
    values = c(
      "#0072B2",
      "#D55E00",
      "#009E73",
      "#CC79A7"
    )
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey40"
    ),
    strip.text = element_text(
      face = "bold",
      size = 15
    ),
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(
      size = 15,
      margin = margin(t = 10)
    ),
    axis.text = element_text(size = 15),
    plot.margin = margin(10, 10, 10, 10)
  )


# Figure 10: Varying treatment effect -------------------------------------

# Evaluate the operating characteristics for different true treatment
# effects while fixing the control-group response probability at 0.42.

delta <- seq(-0.30, 0.30, by = 0.05)

piC <- rep(0.42, length(delta))
piT <- piC + delta

sim_fig10 <- runSimulation(
  piT_grid = piT,
  piC_grid = piC,
  nsim = 10000,
  planning = planning,
  d = d,
  pLU = pLU,
  pUU = pUU,
  pLU2 = pLU2,
  pUU2 = pUU2,
  interim_frac = 0.5,
  seed = 1234
)

names(sim_fig10) <- gsub(
  pattern = "rejection_rate_design([1-4])_interim",
  replacement = "rejection_rate_interim_design\\1",
  x = names(sim_fig10)
)

sim_fig10_long <- sim_fig10 %>%
  rename(
    eN_design1 = eN1,
    eN_design2 = eN2,
    eN_design3 = eN3,
    eN_design4 = eN4
  ) %>%
  pivot_longer(
    cols = matches(
      "^(rejection_rate|rejection_rate_interim|unblinding_rate|eN|fut)_design[1-4]"
    ),
    names_to = c("performance_characteristic", "method"),
    names_pattern = "(.+)_design([1-4])",
    values_to = "value"
  ) %>%
  mutate(
    method = paste0("design", method)
  )

plot_fig10 <- sim_fig10_long %>%
  mutate(
    performance_characteristic = recode(
      performance_characteristic,
      rejection_rate = "Power",
      rejection_rate_interim = "Interim Rejection Rate",
      unblinding_rate = "Unblinding Probability",
      eN = "Expected Sample Size",
      fut = "Futility Probability"
    ),
    method = recode(
      method,
      design1 = "Hybrid Design 1",
      design2 = "Hybrid Design 2",
      design3 = "Single-Stage Design",
      design4 = "Group-Sequential Design"
    )
  )

ref_lines_fig10 <- data.frame(
  performance_characteristic = c(
    "Expected Sample Size",
    "Expected Sample Size"
  ),
  yintercept = c(
    n_control,
    n_control + n_treatment
  ),
  linetype = c(
    "dotted",
    "dotted"
  )
)

fig10 <- plot_fig10 %>%
  ggplot(
    aes(
      x = piT - piC,
      y = value,
      color = method
    )
  ) +
  geom_hline(
    data = ref_lines_fig10,
    aes(yintercept = yintercept),
    linetype = "dashed",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ performance_characteristic,
    scales = "free_y",
    ncol = 2
  ) +
  facetted_pos_scales(
    y = list(
      performance_characteristic == "Power" ~
        scale_y_continuous(limits = c(0, 1)),
      
      performance_characteristic == "Interim Rejection Rate" ~
        scale_y_continuous(limits = c(0, 1)),
      
      performance_characteristic == "Unblinding Probability" ~
        scale_y_continuous(limits = c(0, 1)),
      
      performance_characteristic == "Expected Sample Size" ~
        scale_y_continuous(limits = c(140, 350)),
      
      performance_characteristic == "Futility Probability" ~
        scale_y_continuous(limits = c(0, 1))
    )
  ) +
  labs(
    x = expression("Rate Difference " * pi[T] - pi[C]),
    y = NULL,
    color = NULL,
    linetype = NULL
  ) +
  scale_x_continuous(
    breaks = round(seq(-0.3, 0.3, by = 0.1), 2),
    limits = c(-0.3, 0.3)
  ) +
  scale_color_manual(
    values = c(
      "#0072B2",
      "#D55E00",
      "#009E73",
      "#CC79A7"
    )
  ) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1.4, "cm"),
    strip.background = element_rect(
      fill = "grey92",
      color = "grey40"
    ),
    strip.text = element_text(
      face = "bold",
      size = 15
    ),
    panel.grid.major = element_line(
      color = "grey85",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(
      size = 15,
      margin = margin(t = 10)
    ),
    axis.text = element_text(size = 15),
    plot.margin = margin(10, 10, 10, 10)
  )

