###############################################################################
## Old Faithful data application
## Smoothing-parameter selection for the Szasz--Mirakyan distribution estimator
##
## This script implements:
##   1. the bootstrap-based indirect selector;
##   2. the estimation-based direct selector;
##   3. the nondegeneracy diagnostic required by Corollary 1; and
##   4. the Old Faithful quantile and 20-point-grid analyses.
##
## The indirect selector is not reported when the estimated leading constants
## fail the diagnostic. In particular, a nonpositive V_hat is never replaced by
## machine epsilon.
###############################################################################

## ---------------------------------------------------------------------------
## Basic validation and estimator functions
## ---------------------------------------------------------------------------

validate_sample <- function(X, min_n = 10L) {
  if (!is.numeric(X) || length(X) < min_n ||
      any(!is.finite(X)) || any(X < 0)) {
    stop(
      sprintf(
        "X must be a finite numeric vector on [0, infinity) with at least %d observations.",
        min_n
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_m_grid <- function(m_grid) {
  if (!is.numeric(m_grid) || length(m_grid) < 2L ||
      any(!is.finite(m_grid)) || any(m_grid <= 0)) {
    stop(
      "m_grid must contain at least two finite positive values.",
      call. = FALSE
    )
  }

  integer_error <- abs(m_grid - round(m_grid))
  if (any(integer_error > sqrt(.Machine$double.eps))) {
    stop("Every candidate in m_grid must be an integer.", call. = FALSE)
  }

  m_grid <- sort(unique(as.numeric(round(m_grid))))

  if (length(m_grid) < 2L) {
    stop(
      "m_grid must contain at least two distinct positive values.",
      call. = FALSE
    )
  }

  m_grid
}

F_szasz_emp <- function(x, X, m) {
  validate_sample(X, min_n = 1L)

  if (!is.numeric(x) || length(x) != 1L ||
      !is.finite(x) || x < 0) {
    stop("x must be one finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L ||
      !is.finite(m) || m <= 0) {
    stop("m must be one finite positive number.", call. = FALSE)
  }

  threshold <- ceiling(m * X) - 1
  mean(ppois(threshold, lambda = m * x, lower.tail = FALSE))
}

F_szasz_emp_vec <- function(x, X, m) {
  if (!is.numeric(x) || !is.numeric(m) || length(x) != length(m)) {
    stop("x and m must be numeric vectors of the same length.", call. = FALSE)
  }

  vapply(
    seq_along(x),
    function(i) {
      if (!is.finite(x[i]) || x[i] < 0 ||
          is.na(m[i]) || !is.finite(m[i]) || m[i] <= 0) {
        return(NA_real_)
      }

      F_szasz_emp(x = x[i], X = X, m = m[i])
    },
    numeric(1)
  )
}

F_kernel_cdf_vec <- function(t, X, h) {
  validate_sample(X, min_n = 1L)

  if (!is.numeric(t) || any(!is.finite(t))) {
    stop("t must be a finite numeric vector.", call. = FALSE)
  }

  if (!is.numeric(h) || length(h) != 1L ||
      !is.finite(h) || h <= 0) {
    stop("h must be one finite positive number.", call. = FALSE)
  }

  rowMeans(pnorm(outer(t, X, "-") / h))
}

## ---------------------------------------------------------------------------
## Bootstrap-based indirect selector
## ---------------------------------------------------------------------------

select_m_indirect_boot_grid <- function(
    X,
    x_grid,
    m_grid,
    B = 500L,
    seed = NULL,
    tol_b = 0,
    tol_V = 0,
    warn = TRUE,
    verbose = FALSE
) {
  validate_sample(X)
  m_grid <- validate_m_grid(m_grid)

  if (!is.numeric(x_grid) || length(x_grid) < 1L ||
      any(!is.finite(x_grid)) || any(x_grid <= 0)) {
    stop(
      "x_grid must contain finite strictly positive evaluation points.",
      call. = FALSE
    )
  }

  if (!is.numeric(B) || length(B) != 1L ||
      !is.finite(B) || B < 2 || B != floor(B)) {
    stop("B must be an integer of at least 2.", call. = FALSE)
  }
  B <- as.integer(B)

  for (z in list(tol_b = tol_b, tol_V = tol_V)) {
    if (!is.numeric(z) || length(z) != 1L ||
        !is.finite(z) || z < 0) {
      stop(
        "tol_b and tol_V must be finite nonnegative numbers.",
        call. = FALSE
      )
    }
  }

  n <- length(X)
  L <- length(x_grid)
  J <- length(m_grid)

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or one finite number.", call. = FALSE)
    }
    set.seed(seed)
  }

  ## Store bootstrap samples as multinomial counts. This is algebraically
  ## identical to resampling X with replacement, but avoids recomputing
  ## Poisson probabilities for every bootstrap replication.
  bootstrap_counts <- matrix(0L, nrow = n, ncol = B)
  for (r in seq_len(B)) {
    bootstrap_counts[, r] <- tabulate(
      sample.int(n, size = n, replace = TRUE),
      nbins = n
    )
  }
  bootstrap_weights <- bootstrap_counts / n

  boot_mean <- matrix(NA_real_, nrow = L, ncol = J)
  boot_var <- matrix(NA_real_, nrow = L, ncol = J)

  for (j in seq_along(m_grid)) {
    m <- m_grid[j]
    threshold <- ceiling(m * X) - 1

    tails <- outer(
      m * x_grid,
      threshold,
      function(lambda, cutoff) {
        ppois(cutoff, lambda = lambda, lower.tail = FALSE)
      }
    )

    bootstrap_values <- tails %*% bootstrap_weights
    current_mean <- rowMeans(bootstrap_values)
    centered <- sweep(bootstrap_values, 1L, current_mean, FUN = "-")

    boot_mean[, j] <- current_mean
    boot_var[, j] <- rowMeans(centered^2)

    if (verbose) {
      message(sprintf("Indirect selector: pilot value %d of %d", j, J))
    }
  }

  ## The bootstrap variance is nonnegative by construction. This guard only
  ## removes possible negative round-off at machine precision.
  numerical_floor <- 100 * .Machine$double.eps
  if (any(boot_var < -numerical_floor)) {
    warning(
      "A bootstrap variance was materially negative; inspect the calculation.",
      call. = FALSE
    )
  }
  boot_var <- pmax(boot_var, 0)

  z_bias <- 1 / m_grid
  z_var <- (1 / n) * (1 / sqrt(m_grid))

  out <- data.frame(
    x = as.numeric(x_grid),
    m_indirect = NA_real_,
    m_indirect_cont = NA_real_,
    b_hat = NA_real_,
    V_hat_raw = NA_real_,
    b_condition = FALSE,
    V_condition = FALSE,
    diagnostic_ok = FALSE,
    failure_reason = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(x_grid)) {
    b_fit <- lm.fit(
      x = cbind("(Intercept)" = 1, z_bias = z_bias),
      y = boot_mean[i, ]
    )
    b_hat <- unname(b_fit$coefficients[["z_bias"]])

    V_fit <- lm.fit(
      x = cbind("(Intercept)" = 1, z_var = z_var),
      y = boot_var[i, ]
    )
    V_hat_raw <- -unname(V_fit$coefficients[["z_var"]])

    b_ok <- is.finite(b_hat) && abs(b_hat) > tol_b
    V_ok <- is.finite(V_hat_raw) && V_hat_raw > tol_V
    diagnostic_ok <- b_ok && V_ok

    reason <- character(0)

    if (!is.finite(b_hat)) {
      reason <- c(reason, "non-finite b_hat")
    } else if (abs(b_hat) <= tol_b) {
      reason <- c(reason, "|b_hat| <= tol_b")
    }

    if (!is.finite(V_hat_raw)) {
      reason <- c(reason, "non-finite V_hat")
    } else if (V_hat_raw <= tol_V) {
      reason <- c(reason, "V_hat <= tol_V")
    }

    m_cont <- NA_real_
    m_selected <- NA_real_

    if (diagnostic_ok) {
      m_cont <- n^(2 / 3) *
        (4 * b_hat^2 / V_hat_raw)^(2 / 3)

      if (!is.finite(m_cont) || m_cont <= 0) {
        diagnostic_ok <- FALSE
        reason <- c(reason, "non-finite or nonpositive plug-in selector")
        m_cont <- NA_real_
      } else {
        ## Algorithm 1 rounds upward.
        m_selected <- max(1, ceiling(m_cont))
      }
    }

    out$m_indirect[i] <- m_selected
    out$m_indirect_cont[i] <- m_cont
    out$b_hat[i] <- b_hat
    out$V_hat_raw[i] <- V_hat_raw
    out$b_condition[i] <- b_ok
    out$V_condition[i] <- V_ok
    out$diagnostic_ok[i] <- diagnostic_ok
    out$failure_reason[i] <- if (length(reason) == 0L) {
      NA_character_
    } else {
      paste(reason, collapse = "; ")
    }
  }

  if (warn && any(!out$diagnostic_ok)) {
    failed_x <- paste(
      format(out$x[!out$diagnostic_ok], digits = 6, trim = TRUE),
      collapse = ", "
    )
    warning(
      sprintf(
        paste0(
          "The indirect-selector diagnostic failed at %d evaluation ",
          "point(s): %s. The selector is returned as NA at those points."
        ),
        sum(!out$diagnostic_ok),
        failed_x
      ),
      call. = FALSE
    )
  }

  attr(out, "boot_mean") <- boot_mean
  attr(out, "boot_var") <- boot_var
  attr(out, "m_grid") <- m_grid
  attr(out, "B") <- B
  attr(out, "tol_b") <- tol_b
  attr(out, "tol_V") <- tol_V

  out
}

## Backward-compatible single-point wrapper.
select_m_szasz_single <- function(
    X,
    x_mode = c("quantile", "value"),
    p = NULL,
    x0 = NULL,
    m_grid = seq(5, 300, by = 10),
    B = 500L,
    seed = NULL,
    tol_b = 0,
    tol_V = 0,
    verbose = FALSE
) {
  validate_sample(X)
  x_mode <- match.arg(x_mode)

  if (x_mode == "quantile") {
    if (!is.numeric(p) || length(p) != 1L ||
        !is.finite(p) || p < 0 || p > 1) {
      stop("p must be one probability in [0, 1].", call. = FALSE)
    }
    x0 <- as.numeric(
      quantile(X, probs = p, type = 1, names = FALSE)
    )
  } else if (!is.numeric(x0) || length(x0) != 1L ||
             !is.finite(x0) || x0 <= 0) {
    stop("x0 must be one finite strictly positive value.", call. = FALSE)
  }

  result <- select_m_indirect_boot_grid(
    X = X,
    x_grid = x0,
    m_grid = m_grid,
    B = B,
    seed = seed,
    tol_b = tol_b,
    tol_V = tol_V,
    warn = TRUE,
    verbose = verbose
  )

  list(
    m_star = result$m_indirect[1],
    m_star_x = result$m_indirect_cont[1],
    x0 = result$x[1],
    b_hat = result$b_hat[1],
    V_hat = result$V_hat_raw[1],
    diagnostic_ok = result$diagnostic_ok[1],
    failure_reason = result$failure_reason[1],
    m_grid = attr(result, "m_grid"),
    boot_mean = as.numeric(attr(result, "boot_mean")),
    boot_var = as.numeric(attr(result, "boot_var"))
  )
}

## ---------------------------------------------------------------------------
## Estimation-based direct selector
## ---------------------------------------------------------------------------

select_m_direct_grid <- function(
    X,
    x_grid,
    m_grid,
    h = 0.5 * length(X)^(-1 / 5),
    poisson_tol = 1e-12
) {
  validate_sample(X)
  m_grid <- validate_m_grid(m_grid)

  if (!is.numeric(x_grid) || length(x_grid) < 1L ||
      any(!is.finite(x_grid)) || any(x_grid < 0)) {
    stop(
      "x_grid must contain finite nonnegative evaluation points.",
      call. = FALSE
    )
  }

  if (!is.numeric(h) || length(h) != 1L ||
      !is.finite(h) || h <= 0) {
    stop("h must be one finite positive number.", call. = FALSE)
  }

  if (!is.numeric(poisson_tol) || length(poisson_tol) != 1L ||
      !is.finite(poisson_tol) ||
      poisson_tol <= 0 || poisson_tol >= 1) {
    stop("poisson_tol must lie strictly between 0 and 1.", call. = FALSE)
  }

  n <- length(X)
  L <- length(x_grid)
  J <- length(m_grid)

  Fhn_x <- F_kernel_cdf_vec(x_grid, X, h)

  MSE_hat <- matrix(NA_real_, nrow = L, ncol = J)
  Bias_hat <- matrix(NA_real_, nrow = L, ncol = J)
  Var_hat <- matrix(NA_real_, nrow = L, ncol = J)

  for (j in seq_along(m_grid)) {
    m <- m_grid[j]
    lambda_x <- m * x_grid
    threshold <- ceiling(m * X) - 1

    tails <- outer(
      lambda_x,
      threshold,
      function(lambda, cutoff) {
        ppois(cutoff, lambda = lambda, lower.tail = FALSE)
      }
    )

    Fmn <- rowMeans(tails)
    Var_hat[, j] <- pmax(
      0,
      (rowMeans(tails^2) - Fmn^2) / n
    )

    lambda_max <- max(lambda_x)
    Kmax <- if (lambda_max > 0) {
      max(
        qpois(1 - poisson_tol, lambda = lambda_max),
        ceiling(lambda_max + 10 * sqrt(lambda_max))
      )
    } else {
      0
    }

    k <- 0:as.integer(Kmax)
    tk <- k / m
    Fhk <- F_kernel_cdf_vec(tk, X, h)

    pk <- outer(k, lambda_x, dpois)
    pk_sum <- colSums(pk)

    if (any(!is.finite(pk_sum)) || any(pk_sum <= 0)) {
      stop(
        "Poisson weights underflowed in the direct selector.",
        call. = FALSE
      )
    }

    Bias_hat[, j] <- colSums(pk * Fhk) / pk_sum - Fhn_x
    MSE_hat[, j] <- Var_hat[, j] + Bias_hat[, j]^2
  }

  idx <- max.col(-MSE_hat, ties.method = "first")

  out <- data.frame(
    x = as.numeric(x_grid),
    m_direct = m_grid[idx],
    h = h,
    MSE_hat = MSE_hat[cbind(seq_len(L), idx)],
    Bias_hat = Bias_hat[cbind(seq_len(L), idx)],
    Var_hat = Var_hat[cbind(seq_len(L), idx)],
    at_grid_boundary = idx %in% c(1L, J),
    stringsAsFactors = FALSE
  )

  attr(out, "MSE_grid") <- MSE_hat
  attr(out, "Bias_grid") <- Bias_hat
  attr(out, "Var_grid") <- Var_hat
  attr(out, "m_grid") <- m_grid

  out
}

## Backward-compatible single-point wrapper.
select_m_estimation <- function(
    X,
    x_mode = c("quantile", "value"),
    p = NULL,
    x = NULL,
    m_grid = seq(5, 300, by = 10),
    h = 0.5 * length(X)^(-1 / 5),
    q_type = 1,
    poisson_tol = 1e-12
) {
  validate_sample(X)
  x_mode <- match.arg(x_mode)

  if (x_mode == "quantile") {
    if (!is.numeric(p) || length(p) != 1L ||
        !is.finite(p) || p < 0 || p > 1) {
      stop("p must be one probability in [0, 1].", call. = FALSE)
    }
    x <- as.numeric(
      quantile(X, probs = p, type = q_type, names = FALSE)
    )
  } else if (!is.numeric(x) || length(x) != 1L ||
             !is.finite(x) || x < 0) {
    stop("x must be one finite nonnegative value.", call. = FALSE)
  }

  result <- select_m_direct_grid(
    X = X,
    x_grid = x,
    m_grid = m_grid,
    h = h,
    poisson_tol = poisson_tol
  )

  MSE_grid <- attr(result, "MSE_grid")
  Bias_grid <- attr(result, "Bias_grid")
  Var_grid <- attr(result, "Var_grid")
  m_grid_used <- attr(result, "m_grid")

  list(
    m_opt_check = result$m_direct[1],
    index_min = match(result$m_direct[1], m_grid_used),
    x = result$x[1],
    p = if (x_mode == "quantile") p else NA_real_,
    x_mode = x_mode,
    h = h,
    table = data.frame(
      m = m_grid_used,
      Bias_hat = as.numeric(Bias_grid[1, ]),
      Var_hat = as.numeric(Var_grid[1, ]),
      MSE_hat = as.numeric(MSE_grid[1, ])
    )
  )
}

## ---------------------------------------------------------------------------
## Output helpers
## ---------------------------------------------------------------------------

format_m_tex <- function(z) {
  if (length(z) != 1L || is.na(z) || !is.finite(z)) {
    return("---")
  }

  if (abs(z) >= 1e5) {
    exponent <- floor(log10(abs(z)))
    mantissa <- z / 10^exponent
    return(sprintf("$%.2f\\times 10^{%d}$", mantissa, exponent))
  }

  formatC(z, format = "f", digits = 0)
}

write_old_faithful_table <- function(waiting, eruptions, filename) {
  if (nrow(waiting) != nrow(eruptions)) {
    stop("Waiting-time and eruption-duration tables must have equal rows.")
  }

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0(
      "\\caption{Pointwise selected smoothing parameters for the Old ",
      "Faithful data on 20 equally spaced grid points. The waiting-time ",
      "results are shown in the left panel and the eruption-duration ",
      "results in the right panel.}"
    ),
    "\\label{tab:old-faithful-20grid-selected-m}",
    "\\begin{threeparttable}",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular}{rrrrrrrr}",
    "\\toprule",
    "\\multicolumn{4}{c}{Waiting time} &",
    "\\multicolumn{4}{c}{Eruption duration} \\\\",
    "\\cmidrule(lr){1-4}\\cmidrule(lr){5-8}",
    paste0(
      "Grid point & $x$ & $\\widetilde m_{\\rm ind}(x)$ & ",
      "$\\widetilde m_{\\rm dir}(x)$ & ",
      "Grid point & $x$ & $\\widetilde m_{\\rm ind}(x)$ & ",
      "$\\widetilde m_{\\rm dir}(x)$ \\\\"
    ),
    "\\midrule"
  )

  for (i in seq_len(nrow(waiting))) {
    lines <- c(
      lines,
      sprintf(
        paste0(
          "%d & %.3f & %s & %s & ",
          "%d & %.3f & %s & %s \\\\"
        ),
        i,
        waiting$x[i],
        format_m_tex(waiting$m_indirect[i]),
        format_m_tex(waiting$m_direct[i]),
        i,
        eruptions$x[i],
        format_m_tex(eruptions$m_indirect[i]),
        format_m_tex(eruptions$m_direct[i])
      )
    )
  }

  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}",
    paste0(
      "\\item \\textit{Note}: A dash indicates that the estimated ",
      "nondegeneracy diagnostic for the bootstrap-based indirect ",
      "selector failed. The indirect selector is not reported at that ",
      "evaluation point."
    ),
    "\\end{tablenotes}",
    "\\end{threeparttable}",
    "\\end{table}"
  )

  writeLines(lines, con = filename)
}

plot_quantile_application <- function(
    quantile_results,
    X,
    outcome_label,
    filename,
    fixed_m = c(10, 50, 90, 130, 170)
) {
  probs <- quantile_results$prob
  x <- quantile_results$x

  fixed_values <- vapply(
    fixed_m,
    function(m) {
      vapply(x, F_szasz_emp, numeric(1), X = X, m = m)
    },
    numeric(length(x))
  )

  pdf(filename, width = 6.6, height = 4.8)
  old_par <- par(mar = c(4.2, 4.2, 1.0, 0.8), mgp = c(2.2, 0.7, 0))
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  x_index <- seq_along(probs)

  plot(
    x_index,
    quantile_results$Fhat_indirect,
    type = "n",
    xlim = c(0.75, length(probs) + 0.25),
    ylim = c(0, 1),
    xaxt = "n",
    xlab = "Quantile",
    ylab = "Estimated CDF",
    main = ""
  )
  axis(1, at = x_index, labels = paste0(100 * probs, "%"))

  fixed_pch <- c(1, 2, 0, 5, 6)
  fixed_col <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")

  for (j in seq_along(fixed_m)) {
    points(
      x_index,
      fixed_values[, j],
      pch = fixed_pch[j],
      col = fixed_col[j],
      cex = 1.0
    )
  }

  points(
    x_index,
    quantile_results$Fhat_indirect,
    pch = 16,
    col = "black",
    cex = 1.0
  )
  points(
    x_index,
    quantile_results$Fhat_direct,
    pch = 17,
    col = "gray25",
    cex = 1.0
  )

  legend(
    "topleft",
    legend = c(
      "Indirect selector",
      "Direct selector",
      paste0("m = ", fixed_m)
    ),
    pch = c(16, 17, fixed_pch),
    col = c("black", "gray25", fixed_col),
    bty = "n",
    cex = 0.78,
    ncol = 2
  )

  mtext(outcome_label, side = 3, line = 0.1, font = 2, cex = 0.95)
}

plot_grid_application <- function(results, X, outcome_label, filename) {
  pdf(filename, width = 6.6, height = 4.8)
  old_par <- par(mar = c(4.2, 4.2, 0.8, 0.8), mgp = c(2.2, 0.7, 0))
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  plot(
    ecdf(X),
    verticals = TRUE,
    do.points = FALSE,
    lwd = 1.2,
    col = "gray40",
    main = "",
    xlab = outcome_label,
    ylab = "CDF",
    xlim = range(results$x),
    ylim = c(0, 1)
  )

  ## NA values break the indirect curve at diagnostic failures.
  lines(
    results$x,
    results$Fhat_indirect,
    col = "#D55E00",
    lwd = 2.5
  )
  points(
    results$x[results$diagnostic_ok],
    results$Fhat_indirect[results$diagnostic_ok],
    col = "#D55E00",
    pch = 16,
    cex = 0.55
  )

  lines(
    results$x,
    results$Fhat_direct,
    col = "#0072B2",
    lwd = 2.5,
    lty = 2
  )
  points(
    results$x,
    results$Fhat_direct,
    col = "#0072B2",
    pch = 17,
    cex = 0.55
  )

  legend(
    "bottomright",
    legend = c(
      "Empirical CDF",
      "Indirect selector (valid points)",
      "Direct selector"
    ),
    col = c("gray40", "#D55E00", "#0072B2"),
    lty = c(1, 1, 2),
    lwd = c(1.2, 2.5, 2.5),
    pch = c(NA, 16, 17),
    bty = "n",
    cex = 0.82
  )
}

## ---------------------------------------------------------------------------
## Old Faithful analysis
## ---------------------------------------------------------------------------

analyze_old_faithful_outcome <- function(
    X,
    outcome_name,
    outcome_label,
    output_dir,
    seed,
    B = 500L,
    m_grid = seq(5, 300, by = 10),
    tol_b = 0,
    tol_V = 0
) {
  validate_sample(X)

  probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  x_quantile <- as.numeric(
    quantile(X, probs = probs, type = 1, names = FALSE)
  )
  x_grid <- seq(min(X), max(X), length.out = 20L)

  ## Run each selector once on the union of all evaluation points.
  x_all <- sort(unique(c(x_quantile, x_grid)))

  indirect_all <- select_m_indirect_boot_grid(
    X = X,
    x_grid = x_all,
    m_grid = m_grid,
    B = B,
    seed = seed,
    tol_b = tol_b,
    tol_V = tol_V,
    warn = TRUE
  )

  direct_all <- select_m_direct_grid(
    X = X,
    x_grid = x_all,
    m_grid = m_grid
  )

  extract_results <- function(target_x) {
    idx <- match(target_x, x_all)

    if (anyNA(idx)) {
      stop("Internal matching error for evaluation points.", call. = FALSE)
    }

    out <- cbind(
      indirect_all[idx, , drop = FALSE],
      direct_all[idx, setdiff(names(direct_all), "x"), drop = FALSE]
    )

    out$Fhat_indirect <- F_szasz_emp_vec(
      out$x,
      X,
      out$m_indirect
    )
    out$Fhat_direct <- F_szasz_emp_vec(
      out$x,
      X,
      out$m_direct
    )
    out$EDF <- ecdf(X)(out$x)
    out$outcome <- outcome_label

    out
  }

  quantile_results <- extract_results(x_quantile)
  quantile_results$prob <- probs
  quantile_results$quantile <- paste0(100 * probs, "%")

  quantile_results <- quantile_results[, c(
    "outcome", "quantile", "prob", "x",
    "m_indirect", "m_direct",
    "Fhat_indirect", "Fhat_direct", "EDF",
    "m_indirect_cont", "b_hat", "V_hat_raw",
    "b_condition", "V_condition", "diagnostic_ok",
    "failure_reason", "h", "MSE_hat",
    "Bias_hat", "Var_hat", "at_grid_boundary"
  )]

  grid_results <- extract_results(x_grid)
  grid_results$grid_point <- seq_len(nrow(grid_results))

  grid_results <- grid_results[, c(
    "outcome", "grid_point", "x",
    "m_indirect", "m_direct",
    "Fhat_indirect", "Fhat_direct", "EDF",
    "m_indirect_cont", "b_hat", "V_hat_raw",
    "b_condition", "V_condition", "diagnostic_ok",
    "failure_reason", "h", "MSE_hat",
    "Bias_hat", "Var_hat", "at_grid_boundary"
  )]

  write.csv(
    quantile_results,
    file.path(
      output_dir,
      paste0("OF_quantile_", outcome_name, "_values.csv")
    ),
    row.names = FALSE
  )

  write.csv(
    grid_results,
    file.path(
      output_dir,
      paste0("OF_20grid_", outcome_name, "_values.csv")
    ),
    row.names = FALSE
  )

  plot_quantile_application(
    quantile_results = quantile_results,
    X = X,
    outcome_label = outcome_label,
    filename = file.path(
      output_dir,
      if (outcome_name == "waiting") "OF_waiting.pdf" else "OF_eruption.pdf"
    )
  )

  plot_grid_application(
    results = grid_results,
    X = X,
    outcome_label = outcome_label,
    filename = file.path(
      output_dir,
      paste0("OF_20grid_", outcome_name, "_cdf.pdf")
    )
  )

  list(
    quantile = quantile_results,
    grid = grid_results
  )
}

run_old_faithful_analysis <- function(
    output_dir = file.path(getwd(), "output"),
    B = 500L,
    m_grid = seq(5, 300, by = 10),
    tol_b = 0,
    tol_V = 0,
    waiting_seed = 20260708L,
    eruption_seed = 20260709L
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  data("faithful", package = "datasets")

  waiting <- analyze_old_faithful_outcome(
    X = faithful$waiting,
    outcome_name = "waiting",
    outcome_label = "Waiting time",
    output_dir = output_dir,
    seed = waiting_seed,
    B = B,
    m_grid = m_grid,
    tol_b = tol_b,
    tol_V = tol_V
  )

  eruptions <- analyze_old_faithful_outcome(
    X = faithful$eruptions,
    outcome_name = "eruptions",
    outcome_label = "Eruption duration",
    output_dir = output_dir,
    seed = eruption_seed,
    B = B,
    m_grid = m_grid,
    tol_b = tol_b,
    tol_V = tol_V
  )

  all_grid <- rbind(waiting$grid, eruptions$grid)
  all_quantile <- rbind(waiting$quantile, eruptions$quantile)

  write.csv(
    all_grid,
    file.path(output_dir, "OF_20grid_selected_m_values.csv"),
    row.names = FALSE
  )

  write.csv(
    all_quantile,
    file.path(output_dir, "OF_quantile_selected_m_values.csv"),
    row.names = FALSE
  )

  failures <- all_grid[!all_grid$diagnostic_ok, c(
    "outcome", "grid_point", "x",
    "b_hat", "V_hat_raw", "failure_reason"
  )]

  write.csv(
    failures,
    file.path(output_dir, "OF_indirect_diagnostic_failures.csv"),
    row.names = FALSE
  )

  write_old_faithful_table(
    waiting = waiting$grid,
    eruptions = eruptions$grid,
    filename = file.path(
      output_dir,
      "OF_20grid_selected_m_values.tex"
    )
  )

  capture.output(
    sessionInfo(),
    file = file.path(output_dir, "sessionInfo.txt")
  )

  message("Analysis completed. Output directory: ", normalizePath(output_dir))
  message(
    "Indirect-selector diagnostic failures on the 20-point grids: ",
    nrow(failures)
  )

  invisible(list(
    waiting = waiting,
    eruptions = eruptions,
    failures = failures,
    settings = list(
      B = B,
      m_grid = m_grid,
      tol_b = tol_b,
      tol_V = tol_V,
      waiting_seed = waiting_seed,
      eruption_seed = eruption_seed
    )
  ))
}

## Execute only when this file is run as a script.
if (sys.nframe() == 0L) {
  run_old_faithful_analysis()
}
