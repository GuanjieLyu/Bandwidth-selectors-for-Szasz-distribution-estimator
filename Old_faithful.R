##########################################
## Old faithful data application
##########################################
library(np)
data("faithful")
head(faithful)
summary(faithful$eruptions)
summary(faithful$waiting)

F_szasz_emp <- function(x, X, m) {
  K <- ceiling(m * X)
  mean(ifelse(K >= 1L, 1 - ppois(K - 1L, lambda = m * x), 1))
}

select_m_szasz_single <- function(
    X,
    x_mode = c("quantile", "value"),
    p,
    x0     = NULL,
    m_grid = NULL,
    B,
    seed = NULL,
    verbose = TRUE
) {
  stopifnot(is.numeric(X), all(is.finite(X)), all(X >= 0))
  n <- length(X)
  if (n < 10L) stop("Need at least 10 observations.")
  x_mode <- match.arg(x_mode)
  
  # ---- choose x0 ----
  if (x_mode == "quantile") {
    stopifnot(is.numeric(p), length(p) == 1L, p >= 0, p <= 1)
    x0 <- as.numeric(quantile(X, probs = p,   type = 1, names = FALSE))
  } else { # "value"
    stopifnot(is.numeric(x0), length(x0) == 1L, is.finite(x0), x0 >= 0)
  }
  x0 <- max(0, x0)  # nonnegative support
  stopifnot(length(x0) == 1L)  # F_szasz_emp expects scalar x
  
  # ---- m grid ----
  if (is.null(m_grid)) {
    m_min <- max(3L, floor(0.4 * sqrt(n)))
    m_max <- max(m_min + 2L, ceiling(n^(1/3)))
    m_grid <- unique(round(exp(seq(log(m_min), log(m_max), length.out = 20))))
  }
  m_grid <- sort(unique(pmax(1L, as.integer(m_grid))))
  J <- length(m_grid)
  
  if (!is.null(seed)) set.seed(seed)
  if (verbose) message(sprintf("n=%d, x0=%.6g, |m|=%d, B=%d", n, x0, J, B))
  
  # empirical cdf at x0
  Fn_x0 <- ecdf(X)(x0)
  
  # accumulators over bootstrap and m
  sum_vec   <- numeric(J)
  sumsq_vec <- numeric(J)
  
  for (b in seq_len(B)) {
    Xb <- sample(X, size = n, replace = TRUE)
    for (j in seq_len(J)) {
      m <- m_grid[j]
      # --- use your Poisson-form estimator directly ---
      Fb_x0 <- F_szasz_emp(x0, Xb, m)  # scalar
      sum_vec[j]   <- sum_vec[j]   + Fb_x0
      sumsq_vec[j] <- sumsq_vec[j] + Fb_x0^2
    }
    if (verbose && (b %% max(1L, floor(B/10)) == 0L)) {
      message(sprintf(" bootstrap %d/%d", b, B))
    }
  }
  
  boot_mean <- sum_vec / B
  boot_var  <- (sumsq_vec / B) - boot_mean^2
  boot_var[boot_var < 0] <- 0
  
  # design variables for regressions
  z_bias <- 1 / m_grid
  z_var  <- (1 / n) * (1 / sqrt(m_grid))
  
  # Bias regression through origin at x0
  b_fit  <- lm(boot_mean ~ z_bias)               # <-- intercept included
  b_hat  <- unname(coef(b_fit)["z_bias"])       # slope
  
  # Variance regression with intercept: Var^* ~ a + b * z_var  => V_hat = -b
  fit_var  <- lm(boot_var ~ z_var)
  V_hat <- -unname(coef(fit_var)["z_var"])
  if (!is.finite(V_hat) || V_hat <= 0) V_hat <- .Machine$double.eps
  
  # optimal m at this single x0
  m_star_x <- n^(2/3) * ((4 * b_hat^2 / V_hat)^(2/3))
  m_star   <- max(1L, as.integer(round(m_star_x)))
  
  list(
    m_star    = m_star,
    m_star_x  = m_star_x,
    x0        = x0,
    b_hat     = b_hat,
    V_hat     = V_hat,
    m_grid    = m_grid,
    boot_mean = boot_mean,  # length-J
    boot_var  = boot_var,   # length-J
    Fn_x0     = Fn_x0,
    fit_var   = fit_var
  )
}

F_kernel_cdf <- function(t, X, h, kernel = c("gaussian")) {
  kernel <- match.arg(kernel)
  if (kernel == "gaussian") {
    return( mean(pnorm((t - X) / h)) )
  }
}


select_m_estimation <- function(
    X,
    x_mode = c("quantile", "value"),
    p = p,                 # if x_mode = "quantile"
    x = NULL,                # if x_mode = "value"
    m_grid,
    h = 0.5 * length(X)^(-1/5),
    q_type = 1,
    tol = 1e-12              # Poisson tail tolerance
) {
  stopifnot(is.numeric(X), all(is.finite(X)), all(X >= 0))
  n <- length(X); if (n < 10L) stop("Need at least 10 observations.")
  x_mode <- match.arg(x_mode)
  
  # choose x
  if (x_mode == "quantile") {
    stopifnot(is.numeric(p), length(p) == 1L, p >= 0, p <= 1)
    x <- as.numeric(quantile(X, probs = p, type = q_type, names = FALSE))
  } else {
    stopifnot(is.numeric(x), length(x) == 1L, is.finite(x), x >= 0)
  }
  x <- max(0, x)
  
  # grid
  m_grid <- sort(unique(pmax(1L, as.integer(m_grid))))
  J <- length(m_grid)
  
  Fn    <- ecdf(X)
  Fhn_x <- F_kernel_cdf(x, X, h)
  
  Bias_hat <- Var_hat <- MSE_hat <- numeric(J)
  
  for (j in seq_len(J)) {
    m <- m_grid[j]
    lambda <- m * x
    
    # --------- F_{m,n}(x) via Poisson-tail identity (exact, no truncation) ---------
    K <- ceiling(m * X)
    tail_p <- ifelse(K >= 1L, 1 - ppois(K - 1L, lambda = lambda), 1.0)  # g_i
    Fmn <- mean(tail_p)
    
    # --------- Var_hat from sample (stable) ----------------------------------------
    # Var( mean(g_i) ) = ( E[g^2] - (E[g])^2 ) / n, plug in empirical E
    Var_hat[j] <- max(0, (mean(tail_p^2) - Fmn^2) / n)
    
    # --------- Bias_hat using truncated+renormalized Poisson on k/m ----------------
    # choose Kmax with a CLT buffer to avoid under-truncation at large lambda
    Kmax <- if (lambda > 0) {
      max(qpois(1 - tol, lambda = lambda), ceiling(lambda + 10 * sqrt(lambda)))
    } else 0L
    k  <- 0:Kmax
    pk <- dpois(k, lambda = lambda)
    s  <- sum(pk); if (s == 0) stop("Poisson underflow: increase tol or reduce m/x")
    pk <- pk / s                              # renormalize truncated weights
    
    tk   <- k / m
    Fhk  <- vapply(tk, function(tt) F_kernel_cdf(tt, X, h), numeric(1))
    
    Bias_hat[j] <- sum(Fhk * pk) - Fhn_x
    
    # --------- MSE ---------------------------------------------------------------
    MSE_hat[j] <- Var_hat[j] + Bias_hat[j]^2
  }
  
  idx <- which.min(MSE_hat)
  list(
    m_opt_check = m_grid[idx],
    index_min   = idx,
    x           = x,
    p           = if (x_mode == "quantile") p else NA_real_,
    x_mode      = x_mode,
    h           = h,
    table       = data.frame(m = m_grid,
                             Bias_hat = Bias_hat,
                             Var_hat  = Var_hat,
                             MSE_hat  = MSE_hat)
  )
}



## --- Data ---
X <- faithful$waiting
stopifnot(all(is.finite(X)), all(X >= 0))
n <- length(X)
x_grid <- quantile(X, probs = c(0.1, 0.3, 0.5, 0.7, 0.9), type = 1)
m_grid <- seq(5, 300, by = 10)
# x_grid: numeric vector of evaluation points (nonnegative)
m_hat_vec <- sapply(x_grid, function(x0) {
  select_m_szasz_single(
    X       = X,
    x_mode  = "value",
    p       = NA_real_,     # unused in "value" mode
    x0      = x0,
    m_grid  = m_grid,       # or NULL to auto-build
    B       = 500,
    seed    = NULL,         # same RNG stream; omit for independent draws
    verbose = FALSE
  )$m_star
})


m_hat_vec_est <- sapply(x_grid, function(x0) {
  select_m_estimation(
    X       = X,
    x_mode  = "value",
    x       = x0,          # <-- use x, not x0
    m_grid  = m_grid,
    h       = 0.5 * n^(-1/5)  # or your preferred bandwidth
  )$m_opt_check            # <-- field name returned by your function
})


F_szasz_emp_pairs <- function(x_grid, X, m_hat_vec) {
  stopifnot(is.numeric(x_grid), is.numeric(X), is.numeric(m_hat_vec),
            length(x_grid) == length(m_hat_vec),
            all(is.finite(x_grid)), all(x_grid >= 0),
            all(is.finite(m_hat_vec)), all(m_hat_vec > 0))
  as.numeric(mapply(function(x0, m) F_szasz_emp(x0, X, m),
                    x_grid, m_hat_vec))
}

vals <- F_szasz_emp_pairs(x_grid, X, m_hat_vec)
vals_est <- F_szasz_emp_pairs(x_grid, X, m_hat_vec_est)


F_szasz_emp_grid <- function(x, X, m) {
  stopifnot(is.numeric(x), is.numeric(X), is.numeric(m), length(m) == 1, m > 0)
  K  <- ceiling(m * X)
  Kmax <- if (length(K)) max(K) else 0L
  w  <- tabulate(K + 1L, nbins = Kmax + 1L) / length(X)  # weights w_0,...,w_Kmax
  
  # Build tail probabilities: for each x_j and k, t_{j,k} = 1 - P(Pois(m x_j) <= k-1)
  # Note: for k = 0, ppois(-1, λ) = 0, so tail = 1 (correct).
  xx  <- as.numeric(x)
  kk  <- 0:Kmax
  # Outer constructs an |x|-by-(Kmax+1) matrix of thresholds (k-1) and lambdas (m*x)
  thr <- outer(xx, kk, function(xx, kk) kk - 1)
  lam <- outer(xx, kk, function(xx, kk) m * xx)
  tails <- 1 - ppois(thr, lambda = lam)
  
  as.vector(tails %*% w)            # weighted sum over k for each x
}

vals_10 <- F_szasz_emp_grid(x_grid, X, m=10)
vals_50 <- F_szasz_emp_grid(x_grid, X, m=50)
vals_90 <- F_szasz_emp_grid(x_grid, X, m=90)
vals_130 <- F_szasz_emp_grid(x_grid, X, m=130)
vals_170 <- F_szasz_emp_grid(x_grid, X, m=170)
data.frame(x = x_grid, m_indirect = m_hat_vec, m_direct = m_hat_vec_est, mopt_indirect = vals, mopt_direct = vals_est, m10 = vals_10, 
           m50 = vals_50, m90 = vals_90, m130 = vals_130, m170 = vals_170)







##############################################################################
## Plot figures
##############################################################################
library(ggplot2)
library(tidyr)
library(dplyr)

## df is assumed to already exist like:
df <- data.frame(
  x = x_grid,
  m_indirect = m_hat_vec,
  m_direct   = m_hat_vec_est,
  mopt_indirect = vals,
  mopt_direct   = vals_est,
  m10 = vals_10, m50 = vals_50, m90 = vals_90, m130 = vals_130, m170 = vals_170
)

# Pretty labels for the five quantiles in your x_grid
df$q <- factor(c("10%","30%","50%","70%","90%"),
               levels = c("10%","30%","50%","70%","90%"))

# Long format for the CDF series to plot (dots only)
df_long <- df %>%
  pivot_longer(
    cols = c(mopt_indirect, mopt_direct, m10, m50, m90, m130, m170),
    names_to = "series", values_to = "Fhat"
  )

# Legend labels + shapes/colors (emphasize the two selectors)
series_lab <- c(
  mopt_indirect =  expression("Indirect " * m[opt]),
  mopt_direct   =  expression("Direct " * m[opt]),
  m10  = "m = 10", m50 = "m = 50", m90 = "m = 90", m130 = "m = 130", m170 = "m = 170"
)
shape_map <- c(mopt_indirect=16, mopt_direct=17, m10=1, m50=2, m90=0, m130=5, m170=6)
color_map <- c(mopt_indirect="black", mopt_direct="#333333",
               m10="#1b9e77", m50="#d95f02", m90="#7570b3", m130="#e7298a", m170="#66a61e")

ggplot(df_long, aes(x = q, y = Fhat, color = series, shape = series)) +
  geom_point(position = position_dodge(width = 0.45), size = 3, alpha = 0.95) +
  scale_shape_manual(values = shape_map, labels = series_lab, breaks = names(series_lab)) +
  scale_color_manual(values = color_map,  labels = series_lab, breaks = names(series_lab)) +
  labs(x = "Quantile",
       y = "Estimated CDF",
       color = "Series",
       shape = "Series",
       title = "Szász–Mirakyan CDF Estimates for Old Faithful Waiting Time") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold")  # adjust size here
  )






