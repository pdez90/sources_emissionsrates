# Emission-rate estimators.
# zeros-retained ratio (primary): Q = sum(obs) / sum(pred per unit Q) over every time the site is downwind of the target
#   source (pred > 0.05 ppb per g/s) and the target dominates the modelled footprint (others < dom x target), including
#   downwind times with no detected plume.  Optionally masks times inside events attributed to other sources.
# binned NNLS: joint Hite/Suncor/P66 fit on sums in (wind direction 10 deg x speed class x day/night x site) bins.
# plume-conditioned (Patel et al. 2026): per-plume centerline Q from the peak enhancement; biased high when emissions
#   or detections are intermittent (selection on concentration).
PARTS <- list(Hite = HITE_SPLIT, Suncor = list(Suncor = 1), P66 = list(P66 = 1))
footprints <- function(sites, U, V, hr, parts = PARTS)
  setNames(lapply(sites, function(s) lapply(parts, function(p) footprint(s, U, V, hr, p))), sites)
target_split <- function(F, target) {
  if (target == "Refinery") return(list(pt = F$Suncor + F$P66, oth = F$Hite))
  list(pt = F[[target]], oth = Reduce(`+`, F[setdiff(names(F), target)]))
}
ratio_by_site <- function(enh, FP, target, base, day, gain = CCND_GAIN, dom = 0.3, site_mask = NULL) {
  q <- c(); daily <- list()
  for (s in names(FP)) {
    ts <- target_split(FP[[s]], target); pt <- ts$pt; oth <- ts$oth
    m <- base & (ifelse(is.na(oth), 0, oth) < dom*ifelse(is.na(pt), 0, pt))
    if (!is.null(site_mask[[s]])) m <- m & site_mask[[s]]
    obs <- enh[, s]/(if (s == "BIRCH") 1 else gain)
    dw <- m & is.finite(obs) & is.finite(pt) & pt > 0.05
    if (sum(dw) < 20) next
    d <- data.table(o = obs[dw], p = pt[dw], d = day[dw])[, .(o = sum(o), p = sum(p)), by = d]
    q[s] <- sum(d$o)/sum(d$p)*TPY_PER_GS; daily[[s]] <- d
  }
  list(q = q, daily = daily)
}
boot_site_median <- function(daily, nb = 500, seed = 0) {           # day-block bootstrap of the median over sites
  if (!length(daily)) return(c(NA, NA, NA))
  days <- sort(unique(unlist(lapply(daily, `[[`, "d")))); O <- matrix(0, length(days), length(daily)); Pm <- O
  for (j in seq_along(daily)) { i <- match(daily[[j]]$d, days); O[i, j] <- daily[[j]]$o; Pm[i, j] <- daily[[j]]$p }
  set.seed(seed)
  b <- replicate(nb, { k <- sample.int(length(days), length(days), TRUE); so <- colSums(O[k, , drop = FALSE]); sp <- colSums(Pm[k, , drop = FALSE])
    median(so[sp > 0]/sp[sp > 0])*TPY_PER_GS })
  unname(quantile(b, c(0.025, 0.5, 0.975)))
}
site_spread <- function(q, excl = "CM9") { q <- q[setdiff(names(q), excl)]; q <- q[q > 0]; if (length(q) > 1) sd(log(q))*sqrt((length(q) - 1)/length(q)) else NA }
binned_nnls <- function(enh, FP, base, wf, sp, hr, gain = CCND_GAIN, site_mask = NULL) {
  db <- wf %/% 10; sb <- findInterval(sp, c(1.5, 3)); nb <- as.integer(hr < 7 | hr >= 19); G <- list()
  for (s in names(FP)) {
    F <- FP[[s]]; m <- base & is.finite(enh[, s]) & is.finite(F$Hite); if (!is.null(site_mask[[s]])) m <- m & site_mask[[s]]
    d <- data.table(o = enh[m, s]/(if (s == "BIRCH") 1 else gain), H = F$Hite[m], S = F$Suncor[m], P = F$P66[m], b = db[m], s_ = sb[m], n = nb[m])
    G[[s]] <- d[, .(o = sum(o), H = sum(H), S = sum(S), P = sum(P), c = .N), by = .(b, s_, n)]
  }
  G <- rbindlist(G)[c >= 12]; w <- 1/sqrt(G$c)
  fit <- nnls::nnls(cbind(G$H, G$S, G$P)*w, G$o*w)
  setNames(fit$x*TPY_PER_GS, c("Hite", "Suncor", "P66"))
}
plume_conditioned <- function(EV, SR, b1, target = "Hite", gain = CCND_GAIN, max_ang = 30, min_ws = 0.8) {
  cands <- if (target == "Hite") c("HiteMain", "HiteNorth") else c("Suncor", "P66")
  b1 <- b1[order(t)]; tg <- seq(min(b1$t), max(b1$t), by = 60); h <- series_at(b1$t, b1$h2s, tg)
  e1 <- h - { y <- h; L <- length(y); cs <- c(0, cumsum(!is.na(y)))
              q <- as.numeric(caTools::runquantile(y, 120, 0.1, endrule = "quantile", align = "center")); lo <- pmax(1, seq_len(L) - 60); hi <- pmin(L, seq_len(L) + 59)
              q[(cs[hi + 1] - cs[lo]) < 30] <- NA; q }
  rows <- list()
  for (e in EV[final == target]$ev) for (r in which(SR$ev == e)) {
    a <- unlist(SR[r, paste0("dang_", cands), with = FALSE]); names(a) <- cands; c <- names(which.min(a))
    if (a[c] > max_ang || !(isTRUE(SR$wspd[r] >= min_ws))) next
    tp <- SR$tpeak[r]; site <- SR$site[r]; dist <- bearing(ALL_SITES[[site]], SOURCES[[c]][1:2])["dist"]; hr <- mst_hour(tp)
    if (site == "BIRCH") { k <- which(tg >= tp - 75*60 & tg <= tp + 15*60); dc <- suppressWarnings(max(roll_mean(e1[k], 10, 7), na.rm = TRUE)) } else dc <- SR$dh2s[r]/gain
    fac <- centerline_factor(dist, SR$wspd[r], hr, H = SOURCES[[c]][3])
    rows[[length(rows) + 1]] <- data.table(ev = e, t = tp, site = site, source = c, dist_m = unname(dist), ws = SR$wspd[r], ang = unname(a[c]), dc_ppb = dc, Q_tpy = dc*fac*TPY_PER_GS)
  }
  Q <- rbindlist(rows); if (!nrow(Q)) return(list(Q = Q, s = list()))
  pe <- Q[, .(q = median(Q_tpy)), by = ev]$q; set.seed(0); bs <- replicate(2000, mean(sample(pe, length(pe), TRUE)))
  list(Q = Q, s = list(n_obs = nrow(Q), n_events = length(pe), median = median(pe), mean = mean(pe), boot_lo = unname(quantile(bs, 0.025)), boot_hi = unname(quantile(bs, 0.975))))
}
