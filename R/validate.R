# Birch (reference analyzer) vs CCND: clock, timing, gain/offset, and wind comparisons (Q4 2025 overlap).
series_at <- function(t_src, x, t_dst) x[match(as.numeric(t_dst), as.numeric(t_src))]
# mimic CCND publication: trailing 12-pt mean of 5-min values floored at 0
birch_1h_ccnd_style <- function(b5) roll_mean(pmax(b5$h2s, 0), 12, 12)
best_lag <- function(tx, x, ty, y, lags_min) {
  best <- c(NA, NA)
  for (L in lags_min) {
    yy <- series_at(ty, y, tx + L*60); m <- is.finite(x) & is.finite(yy)
    if (sum(m) > 200) { r <- cor(x[m], yy[m]); if (is.na(best[2]) || r > best[2]) best <- c(L, r) }
  }
  best
}
lag_table <- function(tb, b1h, tc, cm4, start = "2025-10-01", end = "2026-01-01", days = 3) {
  st <- seq(as.POSIXct(start, tz = TZ_MST), as.POSIXct(end, tz = TZ_MST), by = days*86400)
  rbindlist(lapply(st, function(a) { m <- tb >= a & tb <= a + days*86400
    r <- best_lag(tb[m], b1h[m], tc, cm4, seq(-480, 480, 5)); data.table(start = format(a), lag_min = r[1], r = r[2]) }))
}
deming <- function(x, y) {
  mx <- mean(x); my <- mean(y); sxx <- mean((x - mx)^2); syy <- mean((y - my)^2); sxy <- mean((x - mx)*(y - my))
  b <- (syy - sxx + sqrt((syy - sxx)^2 + 4*sxy^2))/(2*sxy); c(int = my - b*mx, slope = b)
}
gain_fit <- function(tb, b1h, tc, cm4, lag = -10, nboot = 500, seed = 1,
                     windows = list(c("2025-10-01", "2025-11-24"), c("2025-12-12", "2026-01-01"))) {
  inw <- Reduce(`|`, lapply(windows, function(w) tb >= as.POSIXct(w[1], tz = TZ_MST) & tb <= as.POSIXct(paste(w[2], "23:59:59"), tz = TZ_MST)))
  x <- b1h[inw]; y <- series_at(tc, cm4, tb[inw] + lag*60); dd <- day_of(tb[inw]); ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]; dd <- dd[ok]
  dm <- deming(x, y); ols <- coef(lm(y ~ x)); hi <- x > 5; sh <- coef(lm(y[hi] ~ x[hi]))[2]
  ud <- unique(dd); idx <- split(seq_along(dd), dd); set.seed(seed)
  bs <- replicate(nboot, { k <- unlist(idx[as.character(sample(ud, length(ud), TRUE))], use.names = FALSE); deming(x[k], y[k])[2] })
  list(n = length(x), r = cor(x, y), deming_slope = unname(dm[2]), deming_int = unname(dm[1]), ols_slope = unname(ols[2]),
       ols_int = unname(ols[1]), slope_birch_gt5 = unname(sh), boot_lo = unname(quantile(bs, 0.025)), boot_hi = unname(quantile(bs, 0.975)))
}
wind_compare <- function(P, b5, sites = c("CM4", "CM8", "CM10", "CM2")) {
  bsp <- sqrt(b5$u^2 + b5$v^2); bd <- wind_from(b5$u, b5$v)
  rbindlist(lapply(sites, function(s) {
    cu <- series_at(P$t, P$u[, s], b5$t); cv <- series_at(P$t, P$v[, s], b5$t); csp <- sqrt(cu^2 + cv^2)
    m <- is.finite(bsp) & is.finite(csp); mm <- m & bsp > 1; dd <- angdiff(wind_from(cu, cv), bd)[mm]
    data.table(site = s, speed_ratio = median(csp[m])/median(bsp[m]), dir_offset_med = median(dd), dir_offset_iqr = unname(diff(quantile(dd, c(.25, .75)))))
  }))
}
# lag of Birch H2S vs its own S-sector (170-245 deg) wind indicator, 1-21 Oct 2025, before and after the clock fix
clock_check <- function(tb, h_raw_1h, h_fix_1h, b5) {
  wd <- wind_from(b5$u, b5$v); sec <- roll_mean(ifelse(is.finite(b5$u), as.numeric(wd >= 170 & wd <= 245), NA), 12, 6)
  m <- tb >= as.POSIXct("2025-10-01", tz = TZ_MST) & tb <= as.POSIXct("2025-10-21 23:59", tz = TZ_MST)
  list(raw = best_lag(tb[m], h_raw_1h[m], b5$t, sec, seq(-720, 720, 15)), fixed = best_lag(tb[m], h_fix_1h[m], b5$t, sec, seq(-720, 720, 15)))
}
