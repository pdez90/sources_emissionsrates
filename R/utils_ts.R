# time-series helpers (regular 5-min / 1-min grids, pandas-like semantics)
roll_mean <- function(x, n, minp = n) {                      # trailing mean with a minimum count of valid values
  s <- frollsum(ifelse(is.finite(x), x, 0), n); k <- frollsum(as.numeric(is.finite(x)), n)
  out <- s/k; out[is.na(k) | k < minp] <- NA; out
}
roll_q10_centered <- function(x, n = 72, minp = n %/% 2, p = 0.10) {   # centred rolling quantile (baseline)
  # pandas rolling(n, center=True) semantics: window for i is [i - n%/%2, i + n - n%/%2 - 1]
  y <- x; y[!is.finite(y)] <- NA; L <- length(y); h <- n %/% 2
  if (all(is.na(y))) return(rep(NA_real_, L))
  lo <- pmax(1, seq_len(L) - h); hi <- pmin(L, seq_len(L) + n - h - 1)
  cs <- c(0, cumsum(!is.na(y))); k <- cs[hi + 1] - cs[lo]
  q <- rep(NA_real_, L)
  R <- as.numeric(caTools::runquantile(y, n, p, endrule = "NA", align = "center"))   # caTools window for j: [j-h+1, j+h]
  int <- (h + 1):(L - h); int <- int[int >= 2 & int <= L]
  q[int] <- R[int - 1]
  for (i in c(seq_len(min(L, h + 1)), seq(max(1, L - h), L))) { w <- y[lo[i]:hi[i]]; if (any(!is.na(w))) q[i] <- quantile(w, p, na.rm = TRUE, names = FALSE, type = 7) }
  q[k < minp] <- NA; q
}
enhancement <- function(M, n = 72) { E <- M; for (j in seq_len(ncol(M))) E[, j] <- M[, j] - roll_q10_centered(M[, j], n); E }
robust_sigma <- function(x) { x <- x[is.finite(x)]; 1.4826*median(abs(x - median(x))) }
bin_end <- function(t, sec) as.POSIXct(ceiling(as.numeric(t)/sec)*sec, origin = "1970-01-01", tz = TZ_MST)
mst_hour <- function(t, back_min = 30) as.integer(format(t - back_min*60, "%H", tz = TZ_MST))
day_of   <- function(t) as.integer(floor(as.numeric(t)/86400 - 7/24) ) # MST calendar day index
