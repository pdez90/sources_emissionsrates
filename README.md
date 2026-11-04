# sources_emissionsrates

Source attribution and emission-rate estimates for hydrogen sulfide (H2S) in Commerce City / north Denver, Colorado. The study area covers the Metro Water Recovery **Robert W. Hite Treatment Facility** and the **Suncor** and **Phillips 66** refinery complex. The estimates come from two networks:

- the Suncor-funded Commerce City–North Denver (**CCND**) community sensor network, 10 sites;
- the CDPHE **MT GOAT Adams County Birch Street** reference site (Teledyne T101, 1 Hz H2S).

All code is in **R**.

> **Status: preliminary research code.** The results below are working estimates for a manuscript in preparation. They are not peer reviewed.

## Quick start

```bash
Rscript install_packages.R        # one time: data.table, R.utils, readxl, jsonlite, caTools, nnls
./run_all.sh --from-processed     # ~2 min, uses the committed data/processed files
./run_all.sh                      # full rebuild from raw files (~15 min); see "Raw data"
```

`scripts/99_vet.R` prints PASS/FAIL for every check and exits non-zero if any check fails. Each step can also be run on its own, for example `Rscript scripts/30_invert.R`, from any directory.

## Layout

```
R/                    library (sourced by every script)
  config.R            paths, site/source coordinates, constants (gain, wind scale, Birch clock fix)
  utils_ts.R          rolling means / centred rolling quantile (pandas semantics), binning
  ccnd_io.R           CCND workbook extraction (readxl) and harmonization -> 5-min MST panel
  birch_io.R          MT GOAT AQDx parsing, QC flags, Birch H2S clock fix -> 1-min series
  wind.R              NOAA ISD (DEN) parsing, CCND vane QA, regional transport wind
  gaussian.R          Gaussian plume kernel (Briggs urban), footprints, centerline factor
  validate.R          Birch vs CCND: timing, gain/offset (Deming + day bootstrap), wind comparison
  detect.R            network plume detection, geometry likelihood, chemistry likelihood ratio
  invert.R            zeros-retained ratio estimator, binned NNLS, plume-conditioned estimator, bootstrap
  results.R           key_numbers.json, CSV writer
scripts/              01-04 build data, 10 validate, 20 detect, 30 invert, 99 vet
data/processed/       compact inputs (committed); every timestamp is fixed MST (UTC-7), interval end
  ccnd_5min_h2s.csv.gz    t_mst, quarter, h2s_CM1..CM10 (ppb, as published: trailing 1-h means, floored at 0)
  ccnd_5min_wind.csv.gz   u_, v_, ws_CM1..CM10 (m/s; u,v = vector towards)
  ccnd_5min_gases.csv.gz  no2, so2, co, tvoc, benz, tol, ebz, xyl x CM1..CM10 (as published)
  birch_1min.csv.gz       Birch H2S (ppb, clock-fixed, QC-filtered, negatives kept), u, v, ws (m/s), 2025-10 .. 2026-06
  den_hourly.csv.gz       Denver Intl. Airport hourly wind
results/              outputs (CSV + key_numbers.json)
tests/expected.json   regression targets for 99_vet.R
```

## Raw data (not committed; see DATA.md)

| Data | Default location (relative to the repo's parent folder) | Source |
|---|---|---|
| CCND quarterly workbooks | `../ccnd_network/` | CCND program website, sensor reports |
| MT GOAT Birch Street zips (`MTGOAT_AdamsCountyBirchStreetBeta_YYYYQ#.zip`) | `../cdphe_birch_st_ADCO/mtgoat/` | CDPHE air toxics repository, HB21-1189 stationary monitoring |
| NOAA ISD DEN (`725650-03017-YYYY`) | `../noaa_isd/` | NOAA NCEI ISD |

Override these locations with `EMIS_RAW_CCND`, `EMIS_RAW_BIRCH` and `EMIS_RAW_DEN`.

## Method summary

1. **Harmonize.**
   - CCND gases are trailing 1-h means of 5-min data, floored at 0, in local time. They are converted to fixed MST.
   - The CCND wind files carry several time-stamp errors, which are fixed:
     - Q4 2023: the Date and Time columns conflict;
     - Q1 2024: stamped in UTC;
     - Q1 2025: AM/PM lost;
     - Q3 2021: stamps collapse to midnight after 2021-09-04 22:02, and those rows are dropped.
   - The Birch H2S stamps before 2025-10-24 01:45 MST are 6 h early: there is a 6.0 h gap at that point, and before it H2S lags Birch's own wind by 330 min. They are shifted.
2. **Validate against the reference.** CM4 sits 150 m from Birch.
   - CM4 lags Birch by about 10 min.
   - Enhancement gain is **0.72** (Deming; day bootstrap 0.63–0.82), with a +2.6 ppb offset.
   - CCND wind speeds are about 0.37–0.48× Birch's.
   - Several CCND vanes drift or jump by 30–90° between quarters (`results/wind_qa_vs_den.csv`). CM8 is stable: −3 to −17° against DEN, and +6° against Birch.
   - **The regional wind is therefore Birch, or CM8 scaled by 1/0.37.**
3. **Detect and attribute plumes** (Q4 2025), following Patel et al. 2025 and 2026.
   - Enhancement is measured above a 6-h centred rolling 10th-percentile baseline. A plume must exceed max(3σ, 2 ppb) for at least 15 min. Plumes are merged across sites within 45 min.
   - The geometry likelihood counts both the sites that saw the plume and the sites that were downwind but silent.
   - A chemistry likelihood ratio (benzene/H2S, CO co-enhancement, NO2/H2S) is trained on plumes whose geometry is unambiguous.
   - Stagnation events (wind below 0.8 m/s, or seen at six or more sites for more than 5 h) are flagged as not plume-resolvable.
4. **Estimate emissions.**
   - *Primary: the zeros-retained ratio.* Q = Σ observed enhancement / Σ Gaussian footprint, summed over every time the site is downwind of the target source and the target dominates the modelled footprint (other sources < 0.3 × target). Downwind times with no detected plume are included as zeros.
   - CCND enhancements are divided by a gain of 0.70. The result is the median over sites (CM9 excluded), with a day-block bootstrap.
   - *Sensitivities:* geometry (Hite main and north polygons, point), gain, the dominance filter, regime, and a joint binned NNLS.
   - *Plume-conditioned* (Patel et al. 2026 centerline): reported as an **upper bound**. Only about 9% of Hite-downwind transport time contains a detected Hite plume, so selecting on plumes biases this estimate high, and it is unstable (heavy tail).

## Key results (t/yr; preliminary)

| Estimate | Hite WWTP | Suncor + P66 |
|---|---|---|
| 2025, zeros-retained, transport (CM8-scaled wind) | **28** (25–31) | **10.5** (9–12) |
| 2025, all hours | 22 | 8 |
| Q4 2025, zeros-retained, transport (Birch wind) | 29 (21–35) | 13 (10–21) |
| Q4 2025, plume-conditioned per-event median (upper bound; unstable) | ~110–170 | ~30 |
| 2024 (diagnostic only; different CCND processing regime) | 6 | 1.5 |
| Reported (CDPHE HB22-1244 TAC, CY2024) | 1.9 (ENGIE cogeneration only) | Suncor 2.3 |

Notes:

- The area-split geometry (main 79% / north 21%) fits far better than north-only (site spread 0.22 vs 0.63).
- CM9 sees a broad NNW sector that is not explained by Hite.

## Main caveats

- The gain comes from one collocated sensor over one quarter, and it is applied network-wide.
- The regional-wind assumption is weakest at night and in stagnation.
- The rolling-percentile background can absorb part of broad plumes.
- Suncor and P66 lie in nearly the same direction from most sites.
- 2024 should not be pooled with 2025 until the processing change in the CCND data is explained.

## Provenance

This R code is a port of a Python prototype. It was checked against that prototype on the same raw files:

| Component | Agreement with the prototype |
|---|---|
| CCND H2S/gas panel | identical |
| Wind panel | identical except 10 values |
| Birch 1-min series | identical |
| 2024–2025 estimates | identical |
| Q4 2025 plume events | the same 290 events; 92% share the same geometry class |

The remaining class differences come from one deliberate change. The R version evaluates wind at the event's central 5-min step; the prototype dropped sites when the midpoint fell between steps.

## References
- Patel, M. Y. et al. (2025) Plume detection and emissions quantification potential using a dense sensor network. *ACS EST Air* 2, 1099–1106.
- Patel, M. Y. et al. (2026) Using network observations to constrain CO and CO2 emissions from an oil refinery. *JGR Atmospheres*.
