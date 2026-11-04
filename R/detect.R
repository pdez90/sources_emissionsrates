# Network plume detection and source attribution (Q4 2025; CCND + Birch), after Patel et al. (2025, 2026).
# Enhancement = value - centred 6-h rolling 10th percentile.  Site plume = H2S enhancement > max(3 sigma, MINPPB) for
# >= 15 min.  Network event = site plumes merged across sites within 45 min.  Geometry likelihood per source uses every
# site with data: p_hit(dA, d) = 0.85 exp(-dA^2 / 2*20^2) min(1, 2 km/d) + 0.03 for enhanced sites, 1 - p_hit for
# silent ones, plus an 'other' class (p = 0.15).  The chemistry likelihood ratio (refinery vs Hite) is trained on
# geometry-clear site plumes.
CAND <- c("HiteMain", "HiteNorth", "Suncor", "P66")
SPECIES <- c("so2", "no2", "co", "tvoc", "benz", "tol", "xyl")
phit <- function(dA, d) 0.85*exp(-dA^2/(2*20^2))*min(1, 2000/d) + 0.03
panel_q4 <- function(P, b5, start = "2025-10-01 00:05:00", end = "2026-01-01 00:00:00") {
  tt <- seq(as.POSIXct(start, tz = TZ_MST), as.POSIXct(end, tz = TZ_MST), by = 300); i <- match(as.numeric(tt), as.numeric(P$t))
  X <- list(); for (k in c("h2s", SPECIES)) if (!is.null(P[[k]])) X[[k]] <- P[[k]][i, , drop = FALSE]
  bh <- roll_mean(series_at(b5$t, b5$h2s, tt), 12, 9)                   # Birch, CCND-style trailing hour (not floored)
  X$h2s <- cbind(X$h2s, BIRCH = bh)
  list(t = tt, X = X)
}
detect_events <- function(tt, X, U, V, min_ppb = 2, min_len = 3, pad_min = 45) {
  E <- lapply(X, enhancement); eh <- E$h2s; sites <- colnames(eh)
  sig <- sapply(sites, function(s) robust_sigma(eh[, s]))
  wf <- wind_from(U, V); ws <- sqrt(U^2 + V^2)
  brg <- lapply(setNames(sites, sites), function(s) lapply(setNames(CAND, CAND), function(c) bearing(ALL_SITES[[s]], SOURCES[[c]][1:2])))
  iv <- list()
  for (s in sites) {
    if (!is.finite(sig[s])) next
    f <- !is.na(eh[, s]) & eh[, s] > max(3*sig[s], min_ppb); r <- rle(f); en <- cumsum(r$lengths); st <- en - r$lengths + 1
    k <- r$values & r$lengths >= min_len
    if (any(k)) iv[[s]] <- data.table(a = st[k], b = en[k], site = s)
  }
  iv <- rbindlist(iv)[order(a, b)]; pad <- pad_min/5
  grp <- integer(nrow(iv)); g <- 0; cur_end <- -Inf
  for (r in seq_len(nrow(iv))) { if (iv$a[r] <= cur_end + pad) { cur_end <- max(cur_end, iv$b[r]) } else { g <- g + 1; cur_end <- iv$b[r] }; grp[r] <- g }
  iv[, ev := grp - 1L]
  rows <- list(); srows <- list(); n <- length(tt)
  for (e in unique(iv$ev)) {
    ie <- iv[ev == e]; t0 <- min(ie$a); t1 <- max(ie$b); tc <- t0 + (t1 - t0) %/% 2
    hit <- list()
    for (r in seq_len(nrow(ie))) { s <- ie$site[r]; seg <- eh[ie$a[r]:ie$b[r], s]; k <- which.max(seg)
      if (is.null(hit[[s]]) || seg[k] > hit[[s]][2]) hit[[s]] <- c(ie$a[r] + k - 1, seg[k]) }
    L <- setNames(numeric(length(CAND)), CAND); Lo <- 0; w0 <- max(1, t0 - pad); w1 <- min(n, t1 + pad)
    for (s in sites) {
      if (!any(is.finite(eh[w0:w1, s]))) next
      ti <- if (!is.null(hit[[s]])) hit[[s]][1] else tc; w <- wf[ti]; if (!is.finite(w)) next
      for (c in CAND) { p <- phit(abs(angdiff(w, brg[[s]][[c]]["bearing"])), brg[[s]][[c]]["dist"]); L[c] <- L[c] + log(if (!is.null(hit[[s]])) p else 1 - p) }
      Lo <- Lo + log(if (!is.null(hit[[s]])) 0.15 else 0.85)
      if (!is.null(hit[[s]])) {
        r <- list(ev = e, site = s, tpeak = tt[ti], dh2s = hit[[s]][2], wfrom = w, wspd = ws[ti])
        for (sp in SPECIES) if (!is.null(E[[sp]]) && s %in% colnames(E[[sp]])) { j <- max(1, ti - 3):min(n, ti + 3); x <- E[[sp]][j, s]
          r[[paste0("d", sp)]] <- if (any(is.finite(x))) mean(x, na.rm = TRUE) else NA_real_ }
        for (c in CAND) r[[paste0("dang_", c)]] <- abs(angdiff(w, brg[[s]][[c]]["bearing"]))
        srows[[length(srows) + 1]] <- as.data.table(r)
      }
    }
    ll <- c(L, other = Lo); post <- exp(ll - max(ll)); post <- post/sum(post)
    hs <- names(hit)[order(sapply(hit, `[`, 1))]
    R <- data.table(ev = e, t0 = tt[t0], t1 = tt[t1], nsites = length(hit), sites = paste(hs, collapse = ","), maxdh2s = max(sapply(hit, `[`, 2)))
    for (c in names(post)) R[[paste0("p_", c)]] <- post[[c]]
    rows[[length(rows) + 1]] <- R
  }
  EV <- rbindlist(rows); SR <- rbindlist(srows, fill = TRUE)
  EV[, `:=`(p_Hite = p_HiteMain + p_HiteNorth, p_Ref = p_Suncor + p_P66)]
  EV[, geom_attrib := fifelse(p_Hite > 0.8, "Hite", fifelse(p_Ref > 0.8, "Refinery", fifelse(p_other > 0.5, "other/local", "ambiguous")))]
  list(EV = EV, SR = SR, sig = sig)
}
chemistry <- function(EV, SR) {
  SR <- copy(SR); SR[, aH := pmin(dang_HiteMain, dang_HiteNorth)][, aR := pmin(dang_Suncor, dang_P66)]
  SR[, geo := fifelse(aH < 15 & aR > 30, "Hite-only", fifelse(aR < 15 & aH > 30, "Refinery-only", "other"))]
  for (sp in SPECIES) if (paste0("d", sp) %in% names(SR)) SR[[paste0("r_", sp)]] <- SR[[paste0("d", sp)]]/SR$dh2s
  tr <- SR[wspd >= TRANSPORT_MIN_WS & geo %in% c("Hite-only", "Refinery-only")][, co_b := as.numeric(!is.na(dco) & dco > 0.05)]
  feats <- list(co = list("co_b", c(-0.5, 0.5, 1.5)), benz = list("r_benz", c(-1e9, 0.01, 0.03, 0.07, 1e9)), no2 = list("r_no2", c(-1e9, 0.5, 1.5, 1e9)))
  hcount <- function(x, b) { x <- x[is.finite(x)]; as.numeric(table(cut(x, b, right = FALSE, include.lowest = TRUE))) }  # numpy.histogram bins
  LR <- lapply(feats, function(f) { ph <- hcount(tr[geo == "Hite-only"][[f[[1]]]], f[[2]]) + 1; pr <- hcount(tr[geo == "Refinery-only"][[f[[1]]]], f[[2]]) + 1
    list(b = f[[2]], l = log((pr/sum(pr))/(ph/sum(ph)))) })
  summ <- tr[, .(r_so2 = median(r_so2, na.rm = TRUE), r_benz = median(r_benz, na.rm = TRUE), r_tvoc = median(r_tvoc, na.rm = TRUE),
                 r_co = median(r_co, na.rm = TRUE), r_no2 = median(r_no2, na.rm = TRUE), frac_co_enh = mean(!is.na(dco) & dco > 0.05)), by = geo]
  f <- SR[, .(co_flag = mean(!is.na(dco) & dco > 0.05), benz_r = median(r_benz, na.rm = TRUE), no2_r = median(r_no2, na.rm = TRUE), wspd = median(wspd, na.rm = TRUE)), by = ev]
  for (k in c("co_flag", "benz_r", "no2_r", "wspd")) set(f, which(is.nan(f[[k]])), k, NA)
  EV <- merge(EV, f, by = "ev", all.x = TRUE)
  lr1 <- function(k, v) { if (!is.finite(v)) return(0); b <- LR[[k]]$b; l <- LR[[k]]$l; l[min(max(findInterval(v, b), 1), length(l))] }
  EV[, chem_llr := mapply(function(c, bz, no) lr1("co", if (is.finite(c)) as.numeric(c >= 0.5) else NA) + lr1("benz", bz) + lr1("no2", no), co_flag, benz_r, no2_r)]
  EV[, regime := fifelse((!is.na(wspd) & wspd < TRANSPORT_MIN_WS) | (nsites >= 6 & as.numeric(t1 - t0, units = "secs") > 5*3600), "stagnation", "transport")]
  EV[, comb_llr := log((p_Ref + 1e-6)/(p_Hite + 1e-6)) + chem_llr]
  EV[, p_Hite_comb := (p_Hite + p_Ref)/(1 + exp(comb_llr))][, p_Ref_comb := (p_Hite + p_Ref) - p_Hite_comb]
  EV[, final := fifelse(regime == "stagnation", "stagnation", fifelse(p_Hite_comb > 0.8, "Hite", fifelse(p_Ref_comb > 0.8, "Refinery", fifelse(p_other > 0.5, "other/local", "ambiguous"))))]
  gc <- EV[regime == "transport" & geom_attrib %in% c("Hite", "Refinery")]
  list(EV = EV, SR = SR, summ = summ, LR = LR, agree = if (nrow(gc)) mean((gc$chem_llr > 0) == (gc$geom_attrib == "Refinery")) else NA)
}
