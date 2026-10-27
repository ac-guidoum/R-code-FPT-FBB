## =============================================================================
## applications.R
## Application of the article: monthly volatility budgets (Section "Application").
## Run from this folder:  Rscript applications.R
##
## The data are not redistributed: they are downloaded by this script.
## =============================================================================

source("packages.R")
source("main_functions.R")
dir.create("results", showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. Screening of candidate series
##    H is estimated on the daily log-closes (log-levels for volatility indices,
##    levels for interest rates, which can be close to 0), over 2016-2025, and
##    separately on each calendar year to assess its stability.
## -----------------------------------------------------------------------------
cand <- data.frame(
  symbol = c("^IRX", "^FVX", "^TNX", "^TYX",                                   # US rates
             "^VIX", "^VXN",                                                   # volatility
             "^BVSP", "^MERV", "^MXX", "^JKSE", "^NSEI", "XU100.IS", "^KLSE",
             "PSEI.PS",                                                        # emerging
             "GC=F", "CL=F", "NG=F", "HG=F", "ZC=F", "ZW=F",                   # commodities
             "BTC-USD", "ETH-USD",                                             # crypto
             "^GSPC"),                                                         # control
  type   = c(rep("rate", 4), rep("volatility", 2), rep("emerging", 8),
             rep("commodity", 6), rep("crypto", 2), "control"),
  scale  = c(rep("level", 4), rep("log", 19)))

get_ohlc <- function(sym) {
  tryCatch(suppressWarnings(quantmod::getSymbols(sym, src = "yahoo", from = "2016-01-01",
                                                 to = "2025-12-31", auto.assign = FALSE)),
           error = function(e) NULL)
}
data <- lapply(cand$symbol, get_ohlc)
names(data) <- cand$symbol
saveRDS(data, "results/screening_data.rds")       # local copy, not redistributed

screen <- do.call(rbind, lapply(seq_len(nrow(cand)), function(i) {
  x <- data[[i]]
  if (is.null(x) || nrow(x) < 500)
    return(data.frame(symbol = cand$symbol[i], type = cand$type[i], days = NA,
                      H_all = NA, H_mean_year = NA, H_min_year = NA, H_max_year = NA,
                      acf1 = NA))
  cl <- as.numeric(quantmod::Cl(x)); ok <- is.finite(cl) & cl > 0
  z  <- if (cand$scale[i] == "log") log(cl[ok]) else cl[ok]
  yr <- as.integer(format(zoo::index(x)[ok], "%Y"))
  hy <- tapply(z, yr, hurst_dv)
  dz <- diff(z)
  data.frame(symbol = cand$symbol[i], type = cand$type[i], days = length(z),
             H_all = hurst_dv(z), H_mean_year = mean(hy), H_min_year = min(hy),
             H_max_year = max(hy),
             acf1 = stats::acf(dz, lag.max = 1, plot = FALSE)$acf[2])
}))
write.csv(screen, "results/screening_hurst.csv", row.names = FALSE)
print(screen, digits = 3)

## -----------------------------------------------------------------------------
## 1b. Screening of persistent daily series (FRED) and agricultural futures
##     Windows of 5 or 20 days will give the open, high and close used later.
##     Robust estimator of H (hurst_dv_robust), whole period and by year.
## -----------------------------------------------------------------------------
cand2 <- data.frame(
  symbol = c("T10Y2Y", "T10Y3M", "T10YIE", "DFII10", "BAA10Y", "DFF",
             "ZC=F", "ZS=F", "ZW=F"),
  src    = c(rep("FRED", 6), rep("yahoo", 3)),
  scale  = c(rep("level", 6), rep("log", 3)))
get_series <- function(sym, src) {
  x <- tryCatch(suppressWarnings(quantmod::getSymbols(sym, src = src, auto.assign = FALSE)),
                error = function(e) NULL)
  if (is.null(x)) return(NULL)
  x <- x[zoo::index(x) >= as.Date("2006-01-01") & zoo::index(x) <= as.Date("2025-12-31")]
  if (src == "yahoo") quantmod::Cl(x) else x[, 1]
}
data2 <- lapply(seq_len(nrow(cand2)), function(i) get_series(cand2$symbol[i], cand2$src[i]))
names(data2) <- cand2$symbol
saveRDS(data2, "results/screening_data2.rds")     # local copy, not redistributed

screen2 <- do.call(rbind, lapply(seq_len(nrow(cand2)), function(i) {
  x <- data2[[i]]
  if (is.null(x)) return(data.frame(symbol = cand2$symbol[i], days = NA, H_robust = NA,
                                    H_year_mean = NA, H_year_sd = NA, H_dv = NA, acf1 = NA,
                                    zero_changes = NA))
  v <- as.numeric(x); d <- zoo::index(x); ok <- is.finite(v) & (cand2$scale[i] == "level" | v > 0)
  z <- if (cand2$scale[i] == "log") log(v[ok]) else v[ok]
  hy <- tapply(z, as.integer(format(d[ok], "%Y")), hurst_dv_robust)
  dz <- diff(z)
  data.frame(symbol = cand2$symbol[i], days = length(z), H_robust = hurst_dv_robust(z),
             H_year_mean = mean(hy, na.rm = TRUE), H_year_sd = sd(hy, na.rm = TRUE), H_dv = hurst_dv(z),
             acf1 = stats::acf(dz, lag.max = 1, plot = FALSE)$acf[2],
             zero_changes = mean(dz == 0))
}))
write.csv(screen2, "results/screening_hurst2.csv", row.names = FALSE)
print(screen2, digits = 3)

## -----------------------------------------------------------------------------
## 1c. Screening of three further leads (robust estimator)
##     A. managed and crawling-peg currencies (Yahoo, log-rates);
##     B. cumulative de-meaned absolute returns of stock indices
##        z_t = sum_{s <= t} (|r_s| - mean|r|), long memory of volatility;
##     C. H at the weekly-to-quarterly scale: series sampled every 5 days,
##        lags of 1 to 12 weeks (5 to 60 days), for all screened series.
## -----------------------------------------------------------------------------
fx  <- c("CNY=X", "INR=X", "TRY=X", "ARS=X", "EGP=X", "VND=X", "NGN=X", "IDR=X",
         "RUB=X", "SAR=X", "HKD=X", "BRL=X", "EURUSD=X")             # BRL, EUR: floating controls
idx <- c("^GSPC", "^NDX", "^STOXX50E", "^N225", "^FTSE", "^HSI")
f3  <- "results/screening_data3.rds"
if (file.exists(f3)) data3 <- readRDS(f3) else {
  data3 <- lapply(c(fx, idx), get_series, src = "yahoo")
  names(data3) <- c(fx, idx)
  saveRDS(data3, f3)                                                 # local copy, not redistributed
}

h_weekly <- function(z) hurst_dv_robust(z[seq(1, length(z), by = 5)], kmax = 12)
screen_row <- function(z, d, symbol, lead) {
  dz <- diff(z)
  hy <- tapply(z, as.integer(format(d, "%Y")), hurst_dv_robust)
  data.frame(symbol = symbol, lead = lead, days = length(z),
             H_daily = hurst_dv_robust(z), H_year_mean = mean(hy, na.rm = TRUE),
             H_year_sd = sd(hy, na.rm = TRUE), H_weekly = h_weekly(z),
             acf1 = stats::acf(dz, lag.max = 1, plot = FALSE)$acf[2],
             zero_changes = mean(dz == 0),
             jump_share = max(abs(dz)) / sum(abs(dz)))
}

screen3 <- do.call(rbind, c(
  lapply(fx, function(s) {
    x <- data3[[s]]; if (is.null(x)) return(NULL)
    v <- as.numeric(x); ok <- is.finite(v) & v > 0
    screen_row(log(v[ok]), zoo::index(x)[ok], s, "A: managed FX")
  }),
  lapply(idx, function(s) {
    x <- data3[[s]]; if (is.null(x)) return(NULL)
    v <- as.numeric(x); ok <- is.finite(v) & v > 0
    a <- abs(diff(log(v[ok])))
    screen_row(cumsum(a - mean(a)), zoo::index(x)[ok][-1], s, "B: cumulated |returns|")
  }),
  lapply(names(data2), function(s) {
    x <- data2[[s]]; if (is.null(x)) return(NULL)
    v <- as.numeric(x); ok <- is.finite(v)
    z <- if (s %in% c("ZC=F", "ZS=F", "ZW=F")) log(v[ok & v > 0]) else v[ok]
    d <- zoo::index(x)[if (s %in% c("ZC=F", "ZS=F", "ZW=F")) ok & v > 0 else ok]
    screen_row(z, d, s, "C: earlier series")
  })))
write.csv(screen3, "results/screening_hurst3.csv", row.names = FALSE)
print(screen3, digits = 3)

## -----------------------------------------------------------------------------
## 2. Application:
## -----------------------------------------------------------------------------
plan(multisession, workers = max(1L, min(parallel::detectCores(), parallel::detectCores() - 1L)))

markets <- c(SP500 = "^GSPC", STOXX50 = "^STOXX50E")
ks    <- c(0.1, 0.2, 0.3, 0.5)   # budgets b = max(z_0, z_20) + k s
nstep <- 20                      # window of 21 trading days
nwin  <- 60                      # rolling window (trading days) for r_t
est_period  <- as.Date(c("2011-01-01", "2018-12-31"))
test_start  <- as.Date("2019-01-01")

## Cumulated excess volatility of a price series x (xts).
excess_vol <- function(x) {
  v <- as.numeric(x); ok <- is.finite(v) & v > 0
  a <- abs(diff(log(v[ok]))); d <- zoo::index(x)[ok][-1]
  m <- c(NA, stats::filter(a, rep(1 / 250, 250), sides = 1))[seq_along(a)]
  e <- a - m; keep <- is.finite(e)
  list(z = cumsum(e[keep]), d = d[keep])
}

## Windows of nstep + 1 days starting every 5 days.
windows <- function(z, d, from) {
  st <- seq(1, length(z) - nstep, by = 5)
  st <- st[d[st] >= from]
  data.frame(pos = st, date = d[st], month = format(d[st], "%Y-%m"),
             z0 = z[st], z1 = z[st + nstep],
             zmax = vapply(st, function(p) max(z[p:(p + nstep)]), 0))
}

## Probability that the maximum on the grid of nstep + 1 points exceeds the budget,
## for the fBB with Hurst index h and daily scale sig (one value per row of g).
prob_window <- function(g, h, sig) {
  sc <- sig * nstep^h
  scen <- data.frame(x0 = 0, y = (g$z1 - g$z0) / sc, a = (g$b - g$z0) / sc, b = 0)
  fit <- mc_fpt(h, T = 1, scen = scen, n = nstep, M = 2e4, K = 20, seed = 11, levels = 1)
  1 - fit$counts[[1]][, nstep + 1] / fit$M
}

app <- list()
for (nm in names(markets)) {
  ev <- excess_vol(data3[[markets[nm]]]); z <- ev$z; d <- ev$d
  est <- d >= est_period[1] & d <= est_period[2]
  H   <- hurst_dv_robust(z[est], kmax = nstep)
  hy  <- tapply(z[est], as.integer(format(d[est], "%Y")), hurst_dv_robust, kmax = nstep)
  seH <- sd(hy) / sqrt(length(hy))
  s   <- stats::mad(diff(z[est], lag = nstep))
  sig <- stats::mad(diff(z[est]))
  w <- windows(z, d, test_start)
  w$r <- vapply(w$pos, function(p) stats::mad(diff(z[(p - nwin):(p - 1)])), 0) / sig
  g <- do.call(rbind, lapply(ks, function(k) cbind(w, k = k)))
  g$top <- pmax(g$z0, g$z1)
  for (mode in c("constant", "rolling")) {
    rr <- if (mode == "constant") rep(1, nrow(g)) else g$r
    gg <- g; gg$b <- gg$top + gg$k * s * rr
    gg$hit    <- as.integer(gg$zmax >= gg$b)
    gg$p_H    <- prob_window(gg, H, sig * rr)
    gg$p_half <- prob_window(gg, 0.5, sig * rr)
    gg$p_Hlo  <- prob_window(gg, max(0.5, H - 1.96 * seH), sig * rr)
    gg$p_Hhi  <- prob_window(gg, min(0.95, H + 1.96 * seH), sig * rr)
    app[[paste(nm, mode)]] <- list(series = nm, mode = mode, H = H, seH = seH, sigma = sig,
                                   windows = nrow(w), grid = gg)
  }
  cat(sprintf("%s: H = %.3f (se %.3f), sigma = %.4g, test windows = %d\n",
              nm, H, seH, sig, nrow(w)))
}
saveRDS(app, "results/application.rds")

## Scores: Brier and logarithmic, fBB with estimated H against the Brownian bridge.
## Paired t-test on the monthly sums of the Brier differences (windows starting in
## the same month and the four budgets are dependent; the months are treated as
## independent).
score <- function(p, y) c(brier = mean((p - y)^2),
                          log = -mean(y * log(pmax(p, 1e-6)) + (1 - y) * log(pmax(1 - p, 1e-6))))
summ <- do.call(rbind, lapply(names(app), function(key) {
  g <- app[[key]]$grid
  do.call(rbind, lapply(c(ks, NA), function(k) {
    gg <- if (is.na(k)) g else g[g$k == k, ]
    sH <- score(gg$p_H, gg$hit); s5 <- score(gg$p_half, gg$hit)
    dB <- tapply((gg$p_half - gg$hit)^2 - (gg$p_H - gg$hit)^2, gg$month, sum)
    data.frame(series = app[[key]]$series, sigma = app[[key]]$mode,
               k = ifelse(is.na(k), "all", k), n = nrow(gg),
               observed = mean(gg$hit), mean_p_H = mean(gg$p_H), mean_p_half = mean(gg$p_half),
               brier_H = sH["brier"], brier_half = s5["brier"],
               log_H = sH["log"], log_half = s5["log"],
               p_value = stats::t.test(dB)$p.value,
               brier_Hlo = score(gg$p_Hlo, gg$hit)["brier"],
               brier_Hhi = score(gg$p_Hhi, gg$hit)["brier"])
  }))
}))
write.csv(summ, "results/application_scores.csv", row.names = FALSE)
print(summ, digits = 4)

calib <- do.call(rbind, lapply(names(app), function(key) {
  g <- app[[key]]$grid; br <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
  do.call(rbind, lapply(c("p_H", "p_half"), function(col) {
    bin <- cut(g[[col]], br, include.lowest = TRUE)
    data.frame(series = app[[key]]$series, sigma = app[[key]]$mode, model = col,
               bin = levels(bin), n = as.vector(table(bin)),
               mean_p = as.vector(tapply(g[[col]], bin, mean)),
               observed = as.vector(tapply(g$hit, bin, mean)))
  }))
}))
write.csv(calib, "results/application_calibration.csv", row.names = FALSE)
print(calib, digits = 3)

## Figure: calibration of the two models, both markets, rolling scale.
cb <- calib[calib$sigma == "rolling" & calib$n >= 20, ]
cb$Model  <- ifelse(cb$model == "p_H", "Fractional Brownian bridge", "Brownian bridge (H = 0.5)")
lab <- c(SP500   = sprintf("S&P 500, H = %.2f", app[["SP500 rolling"]]$H),
         STOXX50 = sprintf("Euro Stoxx 50, H = %.2f", app[["STOXX50 rolling"]]$H))
cb$Market <- factor(lab[cb$series], levels = lab)
fig2 <- ggplot(cb, aes(mean_p, observed, colour = Model, shape = Model, linetype = Model)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dotted") +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  facet_wrap(~ Market) +
  scale_colour_manual(values = c("firebrick", "navy")) +
  scale_shape_manual(values = c(19, 1)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), labels = function(x) sprintf("%.1f", x)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability", y = "Observed frequency",
       colour = NULL, shape = NULL, linetype = NULL) +
  theme_fig()
dir.create("../figures", showWarnings = FALSE)
ggsave("../figures/fig_application.pdf", fig2, width = 9, height = 5)
plan(sequential)
