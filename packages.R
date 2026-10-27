## =============================================================================
## packages.R
## Packages used by the code of the article (Babaci & Guidoum).
## Missing packages are installed from CRAN, then all are loaded.
## =============================================================================

pkgs <- c("future",        # parallel back-ends (sequential, multisession, cluster)
          "future.apply",  # parallel lapply with reproducible random streams
          "quantmod",      # download of daily market data (application)
          "xts",           # time series of prices (application)
          "ggplot2")       # figures

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) install.packages(missing)

library(parallel)          # base R: L'Ecuyer-CMRG streams (nextRNGStream)
invisible(lapply(pkgs, library, character.only = TRUE))
