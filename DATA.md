# Raw data provenance

## CCND community network (CM1–CM10)
- Quarterly Excel workbooks, Q3 2021 – Q4 2025, from the CCND program site's sensor-reports section.
- Gas sheets (H2S, NO2, SO2, CO, tVOC, BTEX) are 1-h rolling averages at 5-min stamps. Sheet names vary by year: "H2S", "H2S 1-hr", "H2S 1-Hour", "H2S-1hr".
- The WS (mph) and WD (deg) sheets are 1-min data.
- The Info sheet says "All data is collected in MST or MDT, as applicable".
- Instrument: SensIT RAMP with an ENVEA Cairsens H2S/CH4S cell. Per the QAPP: detection limit 10 ppb, quantification limit 20 ppb, uncertainty ±30%.
- Issues fixed in `R/ccnd_io.R`:
  - Q4 2023 wind Date/Time conflict;
  - Q1 2024 wind in UTC;
  - Q1 2025 wind AM/PM lost;
  - Q3 2021 wind stamps collapsing to midnight from 2021-09-04 22:02 (those rows dropped).
- Issues not fixed (data excluded or flagged):
  - CM8 Q1 2025;
  - a processing step change between 2024 and 2025;
  - no SO2 sheet in Q3–Q4 2023.

## CDPHE MT GOAT, Adams County Birch Street (39.82811, −104.93644)
- Source: the CDPHE air toxics repository, "Stationary Community Air Toxics Monitoring for HB21-1189 — MT GOAT Adams County Birch Street Site". The quarterly data packets are Google Drive zips (2025Q4, 2026Q1, 2026Q2).
- Format: AQDx v3 weekly CSVs at 1 s, fixed MST.
- H2S: Teledyne T101 UV fluorescence, MDL 0.4 ppb. Zero every 6 h (10 min); bump test at ~600 ppb daily.
- Met: MaxiMet GMX501. Wind speed is in mph (unit code 012).
- Clock issue: H2S stamps before 2025-10-24 01:45:26 MST are 6 h early (fixed in `R/birch_io.R`).

## NOAA ISD, Denver International Airport (USAF 725650, WBAN 03017)
- Standard ISD fixed-width files, 2024–2025, FM-15 reports.
- The 2025 file ends 2025-08-26, when it was downloaded.
