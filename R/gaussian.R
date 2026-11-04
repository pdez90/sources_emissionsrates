# Gaussian plume kernel: Briggs urban dispersion, Pasquill class from wind speed and hour, ground reflection.
xy <- function(lat, lon, lat0 = 39.81, lon0 = -104.94)
  cbind((lon - lon0) * 111320 * cos(lat0 * pi/180), (lat - lat0) * 110540)
# bearing (deg from N) and distance (m) from point `from` to point `to` (lat, lon)
bearing <- function(from, to) {
  a <- xy(from[1], from[2]); b <- xy(to[1], to[2]); dx <- b[1] - a[1]; dy <- b[2] - a[2]
  c(bearing = (atan2(dx, dy) * 180/pi) %% 360, dist = sqrt(dx^2 + dy^2))
}
wind_from <- function(U, V) (atan2(-U, -V) * 180/pi) %% 360        # U,V = vector towards
angdiff   <- function(a, b) ((a - b + 180) %% 360) - 180
sig <- function(x, cls) {                                        # cls 1 = A-B, 2 = C, 3 = D, 4 = E-F
  f  <- sqrt(1/(1 + 0.0004 * x))
  sy <- ifelse(cls <= 1, 0.32*x*f, ifelse(cls == 2, 0.22*x*f, ifelse(cls == 3, 0.16*x*f, 0.11*x*f)))
  sz <- ifelse(cls <= 1, 0.24*x*sqrt(1 + 0.001*x), ifelse(cls == 2, 0.20*x,
         ifelse(cls == 3, 0.14*x/sqrt(1 + 0.0003*x), 0.08*x/sqrt(1 + 0.0015*x))))
  list(sy = sy, sz = sz)
}
stab <- function(sp, hr) {
  day <- hr >= 9 & hr <= 16
  ifelse(day, ifelse(sp < 2, 1, ifelse(sp < 5, 2, 3)), ifelse(sp < 3, 4, 3))
}
# concentration per unit emission (ppb per g/s) at `site` for wind U,V (m/s) and MST hour `hr`.
# parts = named list source -> fraction of the emission
footprint <- function(site, U, V, hr, parts, zr = 3, umin = 0.5) {
  n <- sqrt(U^2 + V^2); sp <- pmax(n, umin); wf <- wind_from(U, V); cls <- stab(sp, hr); out <- numeric(length(U))
  for (src in names(parts)) {
    s <- SOURCES[[src]]; bd <- bearing(ALL_SITES[[site]], s[1:2]); H <- s[3]
    dth <- angdiff(wf, bd["bearing"]) * pi/180
    xd <- bd["dist"] * cos(dth); yd <- bd["dist"] * sin(dth); xd[xd <= 10] <- NA
    g <- sig(xd, cls)
    cc <- parts[[src]]/(2*pi*sp*g$sy*g$sz) * exp(-yd^2/(2*g$sy^2)) * (exp(-(zr - H)^2/(2*g$sz^2)) + exp(-(zr + H)^2/(2*g$sz^2)))
    cc[!is.finite(cc)] <- 0
    out <- out + cc * 1e6 * PPB_PER_UGM3
  }
  out[!is.finite(U) | !is.finite(V)] <- NA
  out
}
# g/s per ppb for a centerline (y = 0) observation at distance `dist` (Patel et al. 2026 per-plume method)
centerline_factor <- function(dist, ws, hr, H = 3, zr = 3) {
  ws <- max(ws, 0.5); g <- sig(dist, stab(ws, hr))
  vert <- exp(-(zr - H)^2/(2*g$sz^2)) + exp(-(zr + H)^2/(2*g$sz^2))
  2*pi*g$sy*g$sz*ws/vert/1e6/PPB_PER_UGM3
}
