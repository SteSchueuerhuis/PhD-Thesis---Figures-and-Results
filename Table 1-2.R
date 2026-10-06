# This script reproduces Table 1 and 2 of the PhD thesis.

# ============================================================
# Pairwise comparison function
# ============================================================

classify <- function(x1, x2) {
  
  stopifnot(is.numeric(x1), is.numeric(x2))
  
  # Mann-Whitney comparison:
  # 1   if X2 is favorable to X1
  # 0.5 if X1 and X2 are tied
  # 0   if X1 is favorable to X2
  
  cX1X2 <- ifelse(
    x2 > x1, 1,
    ifelse(x2 < x1, 0, 0.5)
  )
  return(cX1X2)
}


# ============================================================
# Construct GPC comparison matrix
# ============================================================

getGPCMatrix <- function(X1, X2) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  dX1 <- ncol(X1)
  dX2 <- ncol(X2)
  nX1 <- nrow(X1)
  nX2 <- nrow(X2)
  
  # ----------------------------------------------------------
  # Multiple prioritized outcomes
  # ----------------------------------------------------------
  
  if (dX1 > 1 || dX2 > 1) {
    
    if (dX1 != dX2) {
      stop("Number of outcomes should be equal for both samples.")
    }
    
    H_list <- list()
    
    # Pairwise comparison matrix for each outcome
    for (i in seq_len(dX1)) {
      
      x1 <- X1[, i]
      x2 <- X2[, i]
      
      x1x2 <- expand.grid(x1, x2)
      
      H <- matrix(
        classify(x1x2[, 1], x1x2[, 2]),
        nrow = nX1,
        ncol = nX2
      )
      
      H_list[[paste0("H", i)]] <- H
    }
    
    # Start with highest-priority outcome
    G <- H_list[["H1"]]
    
    # For pairs tied on all previous outcomes,
    # use the next prioritized outcome
    for (i in 2:dX1) {
      
      is_tie <- G == 0.5
      G[is_tie] <- H_list[[paste0("H", i)]][is_tie]
    }
    
  } else {
    
    # --------------------------------------------------------
    # Single outcome
    # --------------------------------------------------------
    
    nX1 <- length(X1)
    nX2 <- length(X2)
    
    xy <- expand.grid(X1, X2)
    
    G <- matrix(
      classify(xy[, 1], xy[, 2]),
      nrow = nX1,
      ncol = nX2
    )
    
    # Important:
    # retain the comparison matrix as the single component
    H_list <- list(H1 = G)
  }
  
  res <- list(
    GPC = G,
    components = H_list
  )
  
  return(res)
}


# ============================================================
# GPC test statistic and confidence interval
# ============================================================

getGPCStats <- function(X1, X2, alpha = 0.05) {
  
  GPC_matrices <- getGPCMatrix(X1, X2)
  G <- GPC_matrices$GPC
  
  # ----------------------------------------------------------
  # Sample sizes
  # ----------------------------------------------------------
  
  nX1 <- nrow(G)
  nX2 <- ncol(G)
  
  N <- nX1 + nX2
  m <- min(nX1, nX2)
  
  dn <- 1 / (
    nX1 * nX2 *
      (nX1 - 1) * (nX2 - 1)
  )
  
  # ----------------------------------------------------------
  # Mann-Whitney effect and tie probability
  # ----------------------------------------------------------
  
  gpp <- sum(G)
  
  # theta =
  # P(X1 < X2) + 0.5 * P(X1 = X2)
  theta <- gpp / (nX1 * nX2)
  
  # Estimated probability of pairwise ties
  tau <- mean(G == 0.5)
  
  # ----------------------------------------------------------
  # Row and column sums
  # ----------------------------------------------------------
  
  gip <- rowSums(G)
  gpj <- colSums(G)
  
  # ----------------------------------------------------------
  # Q quantities
  # ----------------------------------------------------------
  
  Q2 <- sum(
    (gip - gpp / nX1)^2
  )
  
  Q1 <- sum(
    (gpj - gpp / nX2)^2
  )
  
  # ----------------------------------------------------------
  # Unbiased variance estimator
  # ----------------------------------------------------------
  
  sigmaN <- dn * (
    Q1 +
      Q2 -
      nX1 * nX2 *
      (
        theta * (1 - theta) -
          0.25 * tau
      )
  )
  
  # ----------------------------------------------------------
  # C^2 test
  # ----------------------------------------------------------
  
  q <- sigmaN / (theta * (1 - theta))
  
  C2 <- (4 / q) * (theta - 0.5)^2
  
  # Degenerate case
  C2[sigmaN == 0] <-
    4 * m * (theta - 0.5)^2
  
  cAlpha <- qchisq(
    1 - alpha,
    df = 1
  )
  
  # ----------------------------------------------------------
  # Compatible confidence interval
  # ----------------------------------------------------------
  
  ctr <- (
    2 * theta +
      q * cAlpha
  ) / (
    2 * (1 + q * cAlpha)
  )
  
  margin <- (
    sqrt(cAlpha) *
      sqrt(
        4 * q * theta * (1 - theta) +
          q^2 * cAlpha
      )
  ) / (
    2 * (1 + q * cAlpha)
  )
  
  lower_C2 <- ctr - margin
  upper_C2 <- ctr + margin
  
  p_C2 <- 1 - pchisq(
    C2,
    df = 1
  )
  
  # ----------------------------------------------------------
  # Results
  # ----------------------------------------------------------
  
  results_C2 <- data.frame(
    test = "C2",
    nX1 = nX1,
    nX2 = nX2,
    Favorable = sum(G == 1),
    Comparable = sum(G == 0.5),
    Unfavorable = sum(G == 0),
    tau = tau,
    variance = sigmaN,
    SE = sqrt(sigmaN),
    theta = theta,
    lower = lower_C2,
    upper = upper_C2,
    statistic = C2,
    pval = p_C2
  )
  
  return(results_C2)
}


# ============================================================
# Descriptive GPC decomposition + statistical inference
# ============================================================

getGPC <- function(X1, X2, alpha = 0.05) {
  
  G_obj <- getGPCMatrix(X1, X2)
  
  G <- G_obj$GPC
  Hlist <- G_obj$components
  
  k <- length(Hlist)
  n_tot <- length(G)
  
  
  # ----------------------------------------------------------
  # Helper functions
  # ----------------------------------------------------------
  
  # Keep comparisons for the current outcome only where
  # all previous prioritized outcomes remained tied
  cont_filter <- function(prev, cur) {
    
    cur2 <- cur
    
    cur2[
      prev != 0.5 |
        is.na(prev)
    ] <- NA
    
    return(cur2)
  }
  
  # Resolve remaining ties using the next prioritized outcome
  fill05 <- function(base, src) {
    idx <- base == 0.5
    base[idx] <- src[idx]
    return(base)
  }
  
  
  # ----------------------------------------------------------
  # Storage
  # ----------------------------------------------------------
  
  n <- numeric(k)
  
  fav <- character(k)
  unfav <- character(k)
  comp <- character(k)
  
  theta_cumul <- numeric(k)
  theta_contribution <- numeric(k)
  
  
  # ----------------------------------------------------------
  # Cumulative Mann-Whitney effects
  # ----------------------------------------------------------
  base <- Hlist[[1]]
  
  theta_cumul[1] <- mean(base)
  
  if (k >= 2) {
    
    for (j in 2:k) {
      
      base <- fill05(
        base,
        Hlist[[j]]
      )
      
      theta_cumul[j] <- mean(base)
    }
  }
  
  # ----------------------------------------------------------
  # Outcome-specific contributions
  #
  # theta_final =
  # 0.5 + sum(theta_contribution)
  # ----------------------------------------------------------
  
  theta_contribution[1] <-
    theta_cumul[1] - 0.5
  
  # Matrix identifying the comparisons that reach
  # the current prioritized outcome
  prev_filt <- Hlist[[1]]
  
  
  for (j in seq_len(k)) {
    
    if (j == 1) {
      
      cur_filt <- Hlist[[1]]
      
    } else {
      
      cur_filt <- cont_filter(
        prev_filt,
        Hlist[[j]]
      )
      
      prev_filt <- cur_filt
    }
    
    # Number of pairwise comparisons reaching this outcome
    n[j] <- sum(!is.na(cur_filt))
    
    # Favorable comparisons
    fav[j] <- paste0(
      sum(cur_filt == 1, na.rm = TRUE),
      " (",
      round(
        sum(cur_filt == 1, na.rm = TRUE) /
          n_tot * 100,
        2
      ),
      "%)"
    )
    
    # Comparable comparisons
    comp[j] <- paste0(
      sum(cur_filt == 0.5, na.rm = TRUE),
      " (",
      round(
        sum(cur_filt == 0.5, na.rm = TRUE) /
          n_tot * 100,
        2
      ),
      "%)"
    )
    
    # Unfavorable comparisons
    unfav[j] <- paste0(
      sum(cur_filt == 0, na.rm = TRUE),
      " (",
      round(
        sum(cur_filt == 0, na.rm = TRUE) /
          n_tot * 100,
        2
      ),
      "%)"
    )
    
    
    # Contribution of outcomes k >= 2 to the final
    # Mann-Whitney effect
    if (j >= 2) {
      
      theta_contribution[j] <-
        sum(
          cur_filt - 0.5,
          na.rm = TRUE
        ) / n_tot
    }
  }
  
  
  # ----------------------------------------------------------
  # Estimation table
  # ----------------------------------------------------------
  
  estimation <- data.frame(
    k = seq_len(k),
    n = n,
    Favorable = fav,
    Comparable = comp,
    Unfavorable = unfav,
    theta_cumul = theta_cumul,
    theta_contribution = theta_contribution
  )
  
  rownames(estimation) <- NULL
  
  
  # ----------------------------------------------------------
  # Statistical inference
  # ----------------------------------------------------------
  
  testing <- getGPCStats(
    X1,
    X2,
    alpha = alpha
  )
  
  
  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------
  
  return(
    list(
      estimation = estimation,
      testing = testing
    )
  )
}


# ============================================================
# Example data from Evans et al. (2015), Table 2
# ============================================================

intervention <- matrix(
  c(
    2, 1, 1, 2, 3, 2, 3, 3, 1, 3, 2, 1, 2,
    5, 3, 4, 4, 3, 3, 5, 4, 7, 8, 6, 5, 8
  ),
  ncol = 2
)

control <- matrix(
  c(
    3, 2, 1, 2, 3, 2, 1, 2, 3, 1, 1, 2, 3,
    12, 7, 9, 8, 6, 11, 10, 9, 9, 6, 8, 10, 10
  ),
  ncol = 2
)


# =====================================================================
# Strategy 1: Co-Primary Endpoints with Bonferroni-Correction (Table 1)
# =====================================================================

clinical_outcome <- getGPC(
  intervention[, 1],
  control[, 1],
  alpha = 0.025
)

clinical_outcome$estimation
clinical_outcome$testing

antibiotic_use <- getGPC(
  intervention[, 2],
  control[, 2],
  alpha = 0.025
)

antibiotic_use$estimation
antibiotic_use$testing

# ============================================================
# Strategy 2: GPC (Table 2)
# ============================================================

gpc_analysis <- getGPC(
  intervention,
  control,
  alpha = 0.05
)

gpc_analysis$estimation
gpc_analysis$testing
