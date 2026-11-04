# Paths, site/source coordinates and constants.
# Raw-data locations can be overridden with environment variables EMIS_RAW_CCND, EMIS_RAW_BIRCH, EMIS_RAW_DEN.
suppressPackageStartupMessages({ library(data.table) })

find_repo <- function() {
  d <- normalizePath(getwd())
  while (!file.exists(file.path(d, "R", "config.R"))) {
    p <- dirname(d); if (p == d) stop("run from inside the sources_emissionsrates repo"); d <- p
  }
  d
}
REPO      <- if (exists(".repo")) .repo else find_repo()
PROJECT   <- dirname(REPO)                                   # h2s_wwtp_paper/
envpath   <- function(var, default) { v <- Sys.getenv(var); if (nzchar(v)) v else default }
RAW_CCND  <- envpath("EMIS_RAW_CCND",  file.path(PROJECT, "ccnd_network"))
RAW_BIRCH <- envpath("EMIS_RAW_BIRCH", file.path(PROJECT, "cdphe_birch_st_ADCO", "mtgoat"))
RAW_DEN   <- envpath("EMIS_RAW_DEN",   file.path(PROJECT, "noaa_isd"))
INTERIM   <- file.path(REPO, "data", "interim")             # cache (not committed)
PROCESSED <- file.path(REPO, "data", "processed")           # compact inputs (committed)
RESULTS   <- file.path(REPO, "results")
for (p in c(INTERIM, PROCESSED, RESULTS)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

TZ_MST <- "Etc/GMT+7"        # fixed Mountain Standard Time (UTC-7), used for every timestamp in the pipeline

SITES <- list(CM1 = c(39.80164, -104.90882), CM2 = c(39.79630, -104.95727), CM3 = c(39.82736, -104.90193),
              CM4 = c(39.82893, -104.93499), CM5 = c(39.81365, -104.92191), CM6 = c(39.78436, -104.95663),
              CM7 = c(39.80888, -104.91545), CM8 = c(39.81560, -104.94503), CM9 = c(39.78455, -104.96264),
              CM10 = c(39.820268, -104.936616))
SITE_LIST <- names(SITES)
BIRCH     <- c(39.82811, -104.93644)                         # CDPHE MT GOAT Adams County Birch Street
ALL_SITES <- c(SITES, list(BIRCH = BIRCH))
# sources: lat, lon, release height (m). Hite polygons from OSM ways 851439381 (main, 52 ha) and 1482663200 (north, 14 ha)
SOURCES <- list(HiteMain  = c(39.8069,  -104.9549,  3),  HiteNorth = c(39.8140, -104.9561, 3),
                HitePoint = c(39.8100,  -104.9556,  3),  Suncor    = c(39.803333, -104.945556, 20),
                P66       = c(39.79668, -104.94236, 5))
HITE_SPLIT <- list(HiteMain = 52/66, HiteNorth = 14/66)

PPB_PER_UGM3     <- 1/1.19   # H2S at ~83.5 kPa, 15 C (ug/m3 = 1.19 x ppb)
TPY_PER_GS       <- 31.536   # t/yr per g/s
CCND_GAIN        <- 0.70     # CM4 enhancement gain vs Birch reference (Q4 2025); applied to all CCND sites
CM8_SPEED_SCALE  <- 1/0.37   # CM8 5-min vector-mean speed -> Birch-equivalent (Q4 2025)
BIRCH_CLOCK_CUT  <- as.POSIXct("2025-10-24 01:45:26", tz = TZ_MST)   # Birch H2S stamps before this are 6 h early
TRANSPORT_MIN_WS <- 0.8      # m/s (regional wind); below = stagnation

src_files <- function() for (f in c("utils_ts","results","gaussian","ccnd_io","birch_io","wind","validate","detect","invert")) {
  p <- file.path(REPO, "R", paste0(f, ".R")); if (file.exists(p)) source(p) }
