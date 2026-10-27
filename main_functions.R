## =============================================================================
## main_functions.R
## Functions for the article "First Passage of Fractional Brownian Bridges:
## Duality, Error Bounds and Barrier Crossings" (Babaci & Guidoum).
##
## Notation: bridge from x0 to y on [0, T], Hurst index H, boundary
## S(t) = a + b t with x0 < a. Packages are loaded by packages.R.
## =============================================================================


## -----------------------------------------------------------------------------
## 1. Mean and covariance of the bridge
## -----------------------------------------------------------------------------

## Covariance of the fractional Brownian motion, R_H(s, t).
fbm_cov <- function(s, t, H) {
  0.5 * (s^(2 * H) + t^(2 * H) - abs(t - s)^(2 * H))
}

## Mean m(t), covariance K(s, t), variance sigma^2(t) and phi(t) of the bridge.
fbb_moments <- function(H, T, x0, y) {
  phi <- function(t) fbm_cov(t, T, H) / T^(2 * H)
  list(
    phi = phi,
    m   = function(t) x0 + phi(t) * (y - x0),
    K   = function(s, t) fbm_cov(s, t, H) - fbm_cov(s, T, H) * fbm_cov(t, T, H) / T^(2 * H),
    s2  = function(t) t^(2 * H) - fbm_cov(t, T, H)^2 / T^(2 * H)
  )
}


## -----------------------------------------------------------------------------
## 2. Exact simulation on a grid by circulant embedding (Davies-Harte)
## -----------------------------------------------------------------------------

## Eigenvalues of the circulant embedding of the fractional Gaussian noise with
## step T/n. The function stops if the embedding is not non-negative definite.
fgn_eigen <- function(H, n, T = 1) {
  stopifnot(n >= 2)
  j <- 0:n
  gam <- (T / n)^(2 * H) / 2 *
    (abs(j + 1)^(2 * H) - 2 * abs(j)^(2 * H) + abs(j - 1)^(2 * H))
  row <- c(gam[1:n], gam[n + 1], rev(gam[2:n]))
  lam <- Re(stats::fft(row))
  if (min(lam) < -1e-10 * max(lam))
    stop("circulant embedding is not non-negative definite for this (H, n)")
  pmax(lam, 0)
}

## m trajectories of the fractional Brownian motion on the grid k T / n,
## k = 0, ..., n, returned as an (n + 1) x m matrix. One complex FFT gives two
## independent trajectories (real and imaginary parts).
sim_fbm <- function(m, lam, n) {
  N2 <- length(lam)
  k  <- ceiling(m / 2)
  Z  <- matrix(complex(real = stats::rnorm(N2 * k), imaginary = stats::rnorm(N2 * k)),
               N2, k)
  Y  <- stats::mvfft(sqrt(lam / N2) * Z)
  inc <- cbind(Re(Y[1:n, , drop = FALSE]), Im(Y[1:n, , drop = FALSE]))[, 1:m, drop = FALSE]
  rbind(0, apply(inc, 2, cumsum))
}

## Bridge from x0 to y built from fBm trajectories B on the grid tg.
fbb_from_fbm <- function(B, tg, phi, x0, y) {
  x0 + B - outer(phi(tg), x0 + B[nrow(B), ] - y)
}

## Index k of the first grid point with X(t_k) >= a + b t_k, for each column of X;
## NA when the boundary is not crossed on the grid.
first_crossing <- function(X, tg, a, b) {
  above <- X >= (a + b * tg)
  idx <- max.col(t(above), ties.method = "first") - 1L
  idx[colSums(above) == 0] <- NA_integer_
  idx
}


## -----------------------------------------------------------------------------
## 3. Parallel Monte Carlo with nested grids and block streams
## -----------------------------------------------------------------------------

## One block of m trajectories. The trajectories are simulated on the finest grid
## of n * max(levels) steps; the first-passage times on the nested grids of
## n * lev steps are read from the same paths. All the scenarios
## (rows of scen: x0, y, a, b) use the same fBm paths (common random numbers).
## Returns, for each level, a matrix of counts: rows = scenarios, columns 1..nl =
## first crossing at grid index k, last column = no crossing.
mc_block <- function(m, lam, n, T, H, scen, levels) {
  nf <- n * max(levels)
  B  <- sim_fbm(m, lam, nf)
  tf <- T * (0:nf) / nf
  out <- vector("list", length(levels))
  for (l in seq_along(levels)) {
    nl   <- n * levels[l]
    rows <- seq(1, nf + 1, by = max(levels) / levels[l])
    tg   <- tf[rows]
    cnt  <- matrix(0L, nrow(scen), nl + 1)
    for (i in seq_len(nrow(scen))) {
      mo  <- fbb_moments(H, T, scen$x0[i], scen$y[i])
      X   <- fbb_from_fbm(B[rows, , drop = FALSE], tg, mo$phi, scen$x0[i], scen$y[i])
      idx <- first_crossing(X, tg, scen$a[i], scen$b[i])
      cnt[i, ] <- c(tabulate(idx[!is.na(idx)], nbins = nl), sum(is.na(idx)))
    }
    out[[l]] <- cnt
  }
  names(out) <- paste0("n", n * levels)
  out
}

## Streams of random numbers attached to the blocks: block k uses the k-th
## L'Ecuyer-CMRG stream derived from 'seed'.
block_seeds <- function(K, seed) {
  old_kind <- RNGkind()[1]
  on.exit(RNGkind(old_kind))
  set.seed(seed, kind = "L'Ecuyer-CMRG")
  s <- .Random.seed
  seeds <- vector("list", K)
  for (k in seq_len(K)) {
    seeds[[k]] <- s
    s <- parallel::nextRNGStream(s)
  }
  seeds
}

## Monte Carlo estimation of the first-passage law of the bridge.
##   H, T   : Hurst index and horizon
##   scen   : data frame with columns x0, y, a, b (one row per scenario)
##   n      : coarsest grid; the nested grids have n * levels steps
##   M, K   : number of trajectories and number of blocks (M multiple of 2K)
##   seed   : seed of the block streams
## The blocks are distributed with future.apply; the execution plan is chosen by
## the user, e.g. future::plan(future::multisession, workers = 4). The result does
## not depend on the plan or on the number of workers.
mc_fpt <- function(H, T = 1, scen, n, M, K = 100, seed = 1L, levels = c(1, 2, 4)) {
  stopifnot(all(scen$x0 < scen$a), M %% (2 * K) == 0)
  lam   <- fgn_eigen(H, n * max(levels), T)
  seeds <- block_seeds(K, seed)
  res <- future.apply::future_lapply(
    seq_len(K),
    function(k) mc_block(M / K, lam, n, T, H, scen, levels),
    future.seed = seeds,
    future.globals = list(mc_block = mc_block, sim_fbm = sim_fbm,
                          fbb_moments = fbb_moments, fbm_cov = fbm_cov,
                          fbb_from_fbm = fbb_from_fbm, first_crossing = first_crossing,
                          lam = lam, n = n, T = T, H = H, scen = scen,
                          levels = levels, M = M, K = K)
  )
  counts <- Reduce(function(u, v) Map(`+`, u, v), res)
  structure(list(counts = counts, H = H, T = T, scen = scen, n = n, M = M, K = K,
                 levels = levels, seed = seed), class = "fbb_mc")
}


## -----------------------------------------------------------------------------
## 4. Estimators: distribution function, confidence bands, histogram
## -----------------------------------------------------------------------------

## Empirical distribution function of tau_n on the grid of n * level steps,
## with pointwise (Wald) and uniform (DKW-Massart) bands. 'beta' is the bound of
## the discretization bias; it enlarges the upper limit only.
fpt_cdf <- function(fit, scenario = 1, level = 1, alpha = 0.05, beta = 0) {
  cnt <- fit$counts[[match(level, fit$levels)]][scenario, ]
  nl  <- length(cnt) - 1
  M   <- fit$M
  Fh  <- cumsum(cnt[1:nl]) / M
  se  <- sqrt(Fh * (1 - Fh) / M)
  z   <- stats::qnorm(1 - alpha / 2)
  eps <- sqrt(log(2 / alpha) / (2 * M))
  data.frame(t = fit$T * (1:nl) / nl, F = Fh,
             lower_pw = Fh - z * se, upper_pw = Fh + z * se + beta,
             lower_unif = Fh - eps, upper_unif = Fh + eps + beta)
}

## Estimated crossing probability P_T with its standard error.
fpt_prob <- function(fit, scenario = 1, level = 1) {
  cnt <- fit$counts[[match(level, fit$levels)]][scenario, ]
  p <- 1 - cnt[length(cnt)] / fit$M
  c(P_T = p, se = sqrt(p * (1 - p) / fit$M))
}

## Histogram estimator of the density with bins of 'mbin' grid steps
##; the value is attached to the bin centre.
fpt_hist <- function(fit, scenario = 1, level = 1, mbin = 4) {
  cdf <- fpt_cdf(fit, scenario, level)
  nl  <- nrow(cdf)
  h   <- fit$T / nl
  br  <- seq(0, nl - nl %% mbin, by = mbin)
  Fb  <- c(0, cdf$F)[br + 1]
  data.frame(t = (br[-1] - mbin / 2) * h, g = diff(Fb) / (mbin * h))
}


## -----------------------------------------------------------------------------
## 5. Fortet-type (renewal) approximations, for comparison
## -----------------------------------------------------------------------------

## Kernels of the two renewal approximations.
wfa_kernels <- function(H, T, x0, y, a, b) {
  mo <- fbb_moments(H, T, x0, y)
  S  <- function(t) a + b * t
  cond <- function(t, s) {
    k  <- mo$K(s, t)
    s2 <- mo$s2(s)
    list(mu = mo$m(t) + k / s2 * (S(s) - mo$m(s)), v = mo$s2(t) - k^2 / s2)
  }
  list(
    Psi = function(t, s) {
      cl <- cond(t, s); stats::pnorm((S(t) - cl$mu) / sqrt(cl$v), lower.tail = FALSE)
    },
    kappa = function(t, s) {
      cl <- cond(t, s); stats::dnorm(S(t), cl$mu, sqrt(cl$v)) * (t - s)^H
    },
    Q = function(t) stats::pnorm((S(t) - mo$m(t)) / sqrt(mo$s2(t)), lower.tail = FALSE),
    p = function(t) stats::dnorm(S(t), mo$m(t), sqrt(mo$s2(t)))
  )
}

## WFA-c: midpoint rule on the distribution-level equation.
## Returns the cell centres and the approximate density (cell means).
wfa_c <- function(H, T = 1, x0, y, a, b, eps, Tp, N) {
  kr <- wfa_kernels(H, T, x0, y, a, b)
  h  <- (Tp - eps) / N
  ti <- eps + h * (0:N)
  tm <- ti[-1] - h / 2
  g  <- numeric(N)
  for (i in 1:N) {
    j <- seq_len(i - 1)
    acc <- if (i > 1) h * sum(kr$Psi(ti[i + 1], tm[j]) * g[j]) else 0
    den <- h * kr$Psi(ti[i + 1], tm[i])
    g[i] <- if (den > 0) (kr$Q(ti[i + 1]) - acc) / den else 0
  }
  data.frame(t = tm, g = g)
}

## WFA-d: product integration of the density-level equation.
## Returns the grid points and the approximate density.
wfa_d <- function(H, T = 1, x0, y, a, b, eps, Tp, N) {
  kr <- wfa_kernels(H, T, x0, y, a, b)
  h  <- (Tp - eps) / N
  ti <- eps + h * (0:N)
  g  <- numeric(N)
  wii <- h^(1 - H) / (1 - H)
  for (i in 1:N) {
    j <- seq_len(i - 1)
    acc <- 0
    if (i > 1) {
      w <- wii * ((i - j + 1)^(1 - H) - (i - j)^(1 - H))
      acc <- sum(kr$kappa(ti[i + 1], ti[j + 1]) * w * g[j])
    }
    g[i] <- (kr$p(ti[i + 1]) - acc) / ((2 * pi)^(-0.5) * wii)
  }
  data.frame(t = ti[-1], g = g)
}


## Corrected approximation: logit P*_T = logit P_T^c - theta (H - 1/2)^rho, where
## P_T^c is given by wfa_c, and the density of wfa_c rescaled to P*_T.
## theta and rho are estimated in simulations.R.
wfa_corrected <- function(H, T = 1, x0, y, a, b, eps, Tp, N, theta, rho) {
  wc <- wfa_c(H, T, x0, y, a, b, eps, Tp, N)
  pc <- sum(wc$g) * diff(wc$t[1:2])
  dl <- if (H > 0.5) theta * (H - 0.5)^rho else 0
  ps <- if (pc >= 1) 1 else stats::plogis(stats::qlogis(pc) - dl)
  data.frame(t = wc$t, g = wc$g * ps / pc)
}


## -----------------------------------------------------------------------------
## 8. Theme of the figures (ggplot2)
## -----------------------------------------------------------------------------
theme_fig <- function(legend_size = 10, strip_size = 11) {
  theme_minimal() +
    theme(
      panel.grid       = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.75),
      panel.spacing    = unit(0, "lines"),
      legend.position      = "top",
      legend.justification = "centre",
      legend.box           = "horizontal",
      legend.spacing.x     = unit(1, "cm"),
      legend.key.width     = unit(1.2, "cm"),
      legend.text          = element_text(size = legend_size, face = "bold.italic"),
      legend.margin        = margin(t = 0, r = 0, b = -10, l = 0),
      legend.box.spacing   = unit(0.5, "cm"),
      axis.title.x     = element_text(face = "bold", size = 13),
      axis.title.y     = element_text(face = "bold", size = 13),
      axis.text        = element_text(),
      plot.margin      = unit(c(0, 0.05, 0, 0), "cm"),
      strip.text       = element_text(size = strip_size, face = "bold.italic",
                                      margin = margin(2, 1, 2, 1)),
      strip.background = element_rect(fill = "#f0f0f0", color = "black", linewidth = 0.75)
    )
}


## -----------------------------------------------------------------------------
## 7. Estimation of the Hurst index (second-order discrete variations)
## -----------------------------------------------------------------------------

## Estimator of H from a path x observed at equally spaced times, based on the
## second-order differences at lags 1 and 2 (Istas-Lang, Coeurjolly):
##   V_k = mean((x[i + 2k] - 2 x[i + k] + x[i])^2),  H = log2(V_2 / V_1) / 2.
hurst_dv <- function(x) {
  x <- as.numeric(x[is.finite(x)])
  v <- sapply(1:2, function(k) {
    n <- length(x); i <- seq_len(n - 2 * k)
    mean((x[i + 2 * k] - 2 * x[i + k] + x[i])^2)
  })
  0.5 * log2(v[2] / v[1])
}

## Robust version: the mean is replaced by the median of the squared second-order
## differences, which for a Gaussian path is also proportional to k^(2H) and is
## not affected by a few extreme days; H is half the slope of the least-squares
## regression of log(median) on log(k) for the lags k = 1, ..., kmax.
hurst_dv_robust <- function(x, kmax = 10) {
  x <- as.numeric(x[is.finite(x)])
  ks <- seq_len(kmax)
  v <- sapply(ks, function(k) {
    n <- length(x); i <- seq_len(n - 2 * k)
    stats::median((x[i + 2 * k] - 2 * x[i + k] + x[i])^2)
  })
  ok <- is.finite(v) & v > 0
  if (sum(ok) < 3) return(NA_real_)
  unname(stats::coef(stats::lm(log(v[ok]) ~ log(ks[ok])))[2] / 2)
}



## -----------------------------------------------------------------------------
## 6. Closed forms for H = 1/2
## -----------------------------------------------------------------------------

bb_density <- function(t, T = 1, x0, y, a, b) {
  yp <- y - b * T
  (a - x0) / sqrt(2 * pi) / t^1.5 * sqrt(T / (T - t)) *
    exp(-(a - x0)^2 / (2 * t) - (a - yp)^2 / (2 * (T - t)) + (yp - x0)^2 / (2 * T))
}

bb_prob <- function(T = 1, x0, y, a, b) {
  yp <- y - b * T
  if (yp >= a) 1 else exp(-2 * (a - x0) * (a - yp) / T)
}

bb_cdf <- function(t, T = 1, x0, y, a, b) {
  vapply(t, function(u) stats::integrate(bb_density, 0, u, T = T, x0 = x0, y = y,
                                         a = a, b = b)$value, numeric(1))
}
