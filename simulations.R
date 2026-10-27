## =============================================================================
## simulations.R
## Numerical results of the article (Section "Numerical results").
## Run from this folder:  Rscript simulations.R
## =============================================================================

source("packages.R")
source("main_functions.R")

workers <- max(1L, min(parallel::detectCores(), parallel::detectCores() - 1L))
plan(multisession, workers = workers)

## -----------------------------------------------------------------------------
## 1. Checks of the simulation and validation for H = 1/2
## -----------------------------------------------------------------------------

## Law of the simulated fBm: Var B(T) = T^(2H) and lag-1 covariance of the noise.
set.seed(1)
H <- 0.7; n <- 64; T <- 1
lam <- fgn_eigen(H, n, T)
B <- sim_fbm(20000, lam, n)
inc <- diff(B)
cat(sprintf("Var B(T) = %.4f (exact %.4f); lag-1 covariance %.3e (exact %.3e)\n",
            var(B[n + 1, ]), T^(2 * H), mean(inc[-1, ] * inc[-n, ]),
            (T / n)^(2 * H) / 2 * (2^(2 * H) - 2)))

## H = 1/2: Monte Carlo on nested grids against the closed forms.
scen <- data.frame(x0 = c(0, 0), y = c(0.5, 0.5), a = c(1, 1), b = c(0, -0.4))
fit <- mc_fpt(H = 0.5, T = 1, scen = scen, n = 256, M = 2e5, K = 100, seed = 2027)
for (i in seq_len(nrow(scen))) {
  s <- scen[i, ]
  cat(sprintf("scenario %d: exact P_T = %.5f; MC on n = 256, 512, 1024: %s\n", i,
              bb_prob(1, s$x0, s$y, s$a, s$b),
              paste(sprintf("%.5f", sapply(fit$levels, function(l) fpt_prob(fit, i, l)[1])),
                    collapse = ", ")))
}

## Reproducibility: identical counts with 1 and with several R sessions.
small <- function() mc_fpt(H = 0.7, T = 1, scen = scen[1, ], n = 64, M = 4000, K = 20, seed = 7)
plan(sequential); f1 <- small()
plan(multisession, workers = workers); f2 <- small()
cat("identical counts with 1 and", workers, "sessions:", identical(f1$counts, f2$counts), "\n")

## -----------------------------------------------------------------------------
## 2. Main experiment: H = 0.5, ..., 0.9 and three boundaries
##    A: y = 0, S(t) = 0.5;  C: y = 0, S(t) = 0.75 - 0.5 t;  D: y = 0, S(t) = 0.25 + 0.5 t
##    (x0 = 0, T = 1). C and D are dual: P_T(C) = P_T(D) for every H.
##    M = 2e5 paths, nested grids n = 512, 1024, 2048, 4096.
## -----------------------------------------------------------------------------
dir.create("results", showWarnings = FALSE)
Hs <- c(0.5, 0.6, 0.7, 0.8, 0.9)
sc <- data.frame(name = c("A", "C", "D"), x0 = 0, y = 0,
                 a = c(0.5, 0.75, 0.25), b = c(0, -0.5, 0.5))
levels <- c(1, 2, 4, 8)
if (file.exists("results/simulations_mc.rds")) {
  mc <- readRDS("results/simulations_mc.rds")        # reuse a previous run
} else {
  mc <- list()
  for (H in Hs) {
    t0 <- Sys.time()
    mc[[as.character(H)]] <- mc_fpt(H, T = 1, scen = sc[, c("x0", "y", "a", "b")],
                                    n = 512, M = 2e5, K = 400, seed = round(1e4 * H),
                                    levels = levels)
    cat(sprintf("H = %.1f done in %.1f min\n", H, as.numeric(Sys.time() - t0, units = "mins")))
  }
  saveRDS(mc, "results/simulations_mc.rds")
}

## Crossing probabilities: grid estimates, observed order of the bias and
## geometric extrapolation from the nested grids; Fortet-type approximations.
rows <- list()
for (H in Hs) {
  fit <- mc[[as.character(H)]]
  for (i in seq_len(nrow(sc))) {
    s  <- sc[i, ]
    p  <- sapply(levels, function(l) fpt_prob(fit, i, l)[1])
    d  <- diff(p); r <- d[3] / d[2]
    ord <- if (is.finite(r) && r > 0 && r < 1) -log2(r) else NA
    rx  <- if (is.na(ord)) 2^(-H) else r
    fc <- wfa_c(H, 1, s$x0, s$y, s$a, s$b, eps = 0.01, Tp = 0.999, N = 800)
    fd <- wfa_d(H, 1, s$x0, s$y, s$a, s$b, eps = 0.01, Tp = 0.999, N = 800)
    rows[[length(rows) + 1]] <- data.frame(
      H = H, boundary = s$name, P_512 = p[1], P_4096 = p[4],
      se = fpt_prob(fit, i, 8)[2], order = ord, P_ext = p[4] + d[3] * rx / (1 - rx),
      exact = if (H == 0.5) bb_prob(1, s$x0, s$y, s$a, s$b) else NA,
      P_fortet_c = sum(fc$g) * diff(fc$t[1:2]), P_fortet_d = sum(fd$g) * diff(fd$t[1:2]))
  }
}
res <- do.call(rbind, rows)
write.csv(res, "results/simulations_summary.csv", row.names = FALSE)
print(res, digits = 4)

## Figure: first-passage densities (histograms on the grid of 4096 steps, bins
## of 32 steps) for boundaries A and D and H = 0.5, ..., 0.9; exact density for H = 1/2.
## Top row: densities g(t), whose area is P_T; bottom row: conditional densities
## g(t) / P_T, whose area is 1, which show when the crossings occur.
titles <- c(A = "Constant boundary S(t) = 0.5", D = "Increasing boundary S(t) = 0.25 + 0.5t")
rows   <- c("density g(t)", "conditional density g(t) / P_T")
panel  <- function(nm, r) factor(paste0(titles[nm], ": ", rows[r]),
                                 levels = c(paste0(titles["A"], ": ", rows),
                                            paste0(titles["D"], ": ", rows))[c(1, 3, 2, 4)])
dens <- do.call(rbind, lapply(c("A", "D"), function(nm) {
  i <- match(nm, sc$name)
  do.call(rbind, lapply(Hs, function(H) {
    h  <- fpt_hist(mc[[as.character(H)]], i, level = 8, mbin = 32)
    h2 <- fpt_hist(mc[[as.character(H)]], i, level = 8, mbin = 64)
    mass <- sum(h2$g) * diff(h2$t[1:2])
    rbind(data.frame(panel = panel(nm, 1), H = sprintf("H = %.1f", H), t = h$t, g = h$g),
          data.frame(panel = panel(nm, 2), H = sprintf("H = %.1f", H), t = h2$t, g = h2$g / mass))
  }))
}))
tt <- seq(0.005, 0.995, by = 0.005)
exact <- do.call(rbind, lapply(c("A", "D"), function(nm) {
  i <- match(nm, sc$name)
  g <- bb_density(tt, 1, sc$x0[i], sc$y[i], sc$a[i], sc$b[i])
  pT <- bb_prob(1, sc$x0[i], sc$y[i], sc$a[i], sc$b[i])
  rbind(data.frame(panel = panel(nm, 1), t = tt, g = g),
        data.frame(panel = panel(nm, 2), t = tt, g = g / pT))
}))
## crossing probabilities P_T written in the top panels, in the colour of each H
cols <- c("black", "navy", "darkgreen", "darkorange", "firebrick")
ptab <- do.call(rbind, lapply(c("A", "D"), function(nm) {
  p <- res$P_ext[res$boundary == nm]
  top <- max(dens$g[dens$panel == panel(nm, 1)])
  data.frame(panel = panel(nm, 1), H = sprintf("H = %.1f", Hs), t = 0.97,
             g = top * (0.95 - 0.11 * (seq_along(Hs) - 1)),
             label = sprintf("P[T] == '%.3f'", p))
}))
fig1 <- ggplot(dens, aes(t, g, colour = H)) +
  geom_step(linewidth = 0.6) +
  geom_line(data = exact, aes(t, g, linetype = "Exact, H = 0.5"), colour = "black",
            linewidth = 0.5, inherit.aes = FALSE) +
  geom_text(data = ptab, aes(t, g, label = label, colour = H), parse = TRUE,
            hjust = 1, size = 3.2, show.legend = FALSE) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = cols) +
  scale_linetype_manual(values = "dashed") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(labels = function(x) sprintf("%.1f", x)) +
  guides(colour   = guide_legend(order = 1, override.aes = list(linewidth = 1.5)),
         linetype = guide_legend(order = 2, override.aes = list(linewidth = 0.8))) +
  labs(x = "Time t", y = NULL, colour = "Monte Carlo:", linetype = NULL) +
  theme_fig(strip_size = 9) +
  theme(legend.title = element_text(face = "bold", size = 10))
dir.create("../figures", showWarnings = FALSE)
ggsave("../figures/fig_densities.pdf", fig1, width = 9, height = 6.5)

## -----------------------------------------------------------------------------
## 3. Reference values for H > 1/2
##    P_T with an interval: lower end = P^(4096) - eps_M (guaranteed at level
##    1 - alpha by the bias sign and the DKW inequality), upper end = P_ext + eps_M
##    (extrapolated bias, heuristic); mean and median of tau given tau <= T.
## -----------------------------------------------------------------------------
alpha <- 0.05
epsM <- sqrt(log(2 / alpha) / (2 * 2e5))
ref <- do.call(rbind, lapply(Hs[-1], function(H) {
  fit <- mc[[as.character(H)]]
  do.call(rbind, lapply(seq_len(nrow(sc)), function(i) {
    r <- res[res$H == H & res$boundary == sc$name[i], ]
    cnt <- fit$counts[[match(8, fit$levels)]][i, ]
    nl <- length(cnt) - 1; tt <- (1:nl) / nl; w <- cnt[1:nl] / sum(cnt[1:nl])
    data.frame(H = H, boundary = sc$name[i], P_ext = r$P_ext,
               lower = r$P_4096 - epsM, upper = r$P_ext + epsM,
               mean_tau = sum(tt * w), median_tau = tt[which(cumsum(w) >= 0.5)[1]])
  }))
}))
write.csv(ref, "results/reference_values.csv", row.names = FALSE)
print(ref, digits = 4)

## -----------------------------------------------------------------------------
## 4. Corrected Fortet approximation (distribution form)
##    logit P* = logit P^c - theta (H - 1/2)^rho, density rescaled to P*.
##    Calibration: weighted nonlinear least squares on the logit errors of the
##    boundaries A, C, D for H = 0.6, ..., 0.9 (weights = 1 / Var(logit P)).
##    Validation: new Hurst indices and new boundaries, new Monte Carlo runs.
## -----------------------------------------------------------------------------
logit <- function(p) log(p / (1 - p))
cal <- res[res$H > 0.5, ]
cal$shift <- logit(cal$P_fortet_c) - logit(cal$P_ext)
cal$w <- 2e5 * cal$P_ext * (1 - cal$P_ext)
fit_d <- nls(shift ~ theta * (H - 0.5)^rho, data = cal, weights = w,
             start = list(theta = 1, rho = 0.6))
co <- summary(fit_d)$coefficients
cat(sprintf("calibration: theta = %.4f (se %.4f), rho = %.4f (se %.4f), residual se %.2f\n",
            co[1, 1], co[1, 2], co[2, 1], co[2, 2], summary(fit_d)$sigma))
write.csv(data.frame(theta = co[1, 1], se_theta = co[1, 2], rho = co[2, 1], se_rho = co[2, 2],
                     resid_se = summary(fit_d)$sigma), "results/correction_fit.csv", row.names = FALSE)

sv <- data.frame(name = c("V1", "V2", "V3", "V4"), x0 = c(0, 0, 0, 0.2),
                 y = c(-0.25, 0.1, 0, 0.4), a = c(0.5, 0.4, 1, 0.6), b = c(0, 0.3, -0.8, 0))
Hv <- c(0.55, 0.65, 0.75, 0.85)
if (file.exists("results/validation_mc.rds")) {
  mcv <- readRDS("results/validation_mc.rds")
} else {
  mcv <- list()
  for (H in Hv) mcv[[as.character(H)]] <- mc_fpt(H, T = 1, scen = sv[, c("x0", "y", "a", "b")],
                                                 n = 512, M = 2e5, K = 400,
                                                 seed = round(1e4 * H) + 7, levels = levels)
  saveRDS(mcv, "results/validation_mc.rds")
}
val <- do.call(rbind, lapply(Hv, function(H) {
  fit <- mcv[[as.character(H)]]
  do.call(rbind, lapply(seq_len(nrow(sv)), function(i) {
    s <- sv[i, ]
    p <- sapply(levels, function(l) fpt_prob(fit, i, l)[1]); d <- diff(p); r <- d[3] / d[2]
    if (!is.finite(r) || r <= 0 || r >= 1) r <- 2^(-H)
    pm <- p[4] + d[3] * r / (1 - r)
    fc <- wfa_c(H, 1, s$x0, s$y, s$a, s$b, eps = 0.01, Tp = 0.999, N = 800)
    fd <- wfa_d(H, 1, s$x0, s$y, s$a, s$b, eps = 0.01, Tp = 0.999, N = 800)
    fs <- wfa_corrected(H, 1, s$x0, s$y, s$a, s$b, eps = 0.01, Tp = 0.999, N = 800,
                        theta = co[1, 1], rho = co[2, 1])
    pc <- sum(fc$g) * diff(fc$t[1:2]); pd <- sum(fd$g) * diff(fd$t[1:2])
    ps <- sum(fs$g) * diff(fs$t[1:2])
    data.frame(H = H, boundary = s$name, P_mc = pm, se = fpt_prob(fit, i, 8)[2],
               P_c = pc, P_d = pd, P_corr = ps,
               err_c = (pc - pm) / pm, err_d = (pd - pm) / pm, err_corr = (ps - pm) / pm)
  }))
}))
write.csv(val, "results/correction_validation.csv", row.names = FALSE)
print(val, digits = 3)
cat(sprintf("mean |relative error|: Fortet-c %.3f, Fortet-d %.3f, corrected %.3f\n",
            mean(abs(val$err_c)), mean(abs(val$err_d)), mean(abs(val$err_corr))))
cat(sprintf("max  |relative error|: Fortet-c %.3f, Fortet-d %.3f, corrected %.3f\n",
            max(abs(val$err_c)), max(abs(val$err_d)), max(abs(val$err_corr))))

plan(sequential)
