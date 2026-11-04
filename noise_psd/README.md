# noise_psd/ — Documentation

## I. Overview and Scientific Purpose

This sub-package computes the **Power Spectral Density (PSD) of DAS noise** for the Kefalonia deployment (April–June 2024), then overlays meteorological data (ERA5 wind speed and significant wave height) to diagnose noise source mechanisms. Outputs are used by the `check_minmag` package to assess the minimum magnitude of completeness (Mc) and validate detection thresholds.

**Position in the package**:

```
check_minmag_pkg/
└── noise_psd/           ← this folder
    ├── psd_daily_h5/    (NOT bundled — must recompute from DAS data)
    └── ...scripts...

Products/
└── noise_trends_f8-20Hz_p5.h5   ← final product, already bundled
```

**Script pipeline** (Script 1 must run first; Scripts 2/3/4 are independent consumers):

```
Script 1 → psd_daily_h5/   → Script 2 (spectrogram + meteorology overlay plot)
                            → Script 3 (PPSD 2D histogram)
                            → Script 4 (noise trend analysis) → products/noise_trends_f8-20Hz_p5.h5
```

---

## II. Script-by-Script Description

### Script 1: `1_calculate_psd_h5.py` — DAS → Daily PSD h5 Files

- **Purpose**: Read continuous DAS data day-by-day, divide into 40-minute chunks, compute Welch PSD for each chunk, and save as daily h5 files.
- **Inputs**:
  - Continuous DAS data: `DAS_DATA_DIR = os.path.join(_PKG_ROOT, 'IN_DAS_data') + os.sep` (explicit top-level variable; package-relative default — place or symlink the 80 Hz archive `DATA_DEC_80Hz_16m_SORTED/` there, or change this line before running).
  - Key parameters: `fs=80`, `nperseg_sec=200` (Welch segment length 200 s → nperseg=16000 samples), `chunk_sec=2400` (chunk = 40 min = 2400 s), 55 channels (`das_dist_start=0`, `das_dist_end=900`, every channel ≈ 16.34 m spacing → ≈ 0–14.7 km of fiber).
- **Outputs**: `psd_daily_h5/<YYYYMMDD>.h5`, one file per day. Structure of each h5:
  - `frequency`: dataset, shape `(8001,)`, the frequency axis 0–40 Hz (Δf = 0.005 Hz, from 80 Hz / nperseg=16000).
  - `ch_X.Xkm/`: one group per channel (named by distance in km); each group contains datasets named by UTC timestamp `HH:MM:SS` (40-minute interval); each dataset is shape `(8001,)` float64, PSD in units `(strain-rate)²/Hz`.
- **Dependencies**: `functions/load_das_h5_PE2024.py`'s `load_das_custom`; standard library + `numpy` + `scipy.signal.welch`.
- **Note**: `psd_daily_h5/` is **not bundled** in this package (tens of GB); must be recomputed from raw DAS data. This is the only script that requires DAS access; Scripts 2/3/4 only need `psd_daily_h5/`.

---

### Script 2: `2_plot_results_with_wind_wave.py` — Spectrogram + Meteorological Overlay

- **Purpose**: Read `psd_daily_h5/` and ERA5 meteorological data; plot a multi-panel spectrogram time series with overlaid wind speed and significant wave height.
- **Inputs**:
  - `psd_daily_h5/` (Script 1 output).
  - `meteo/load_kefa_meteo.py`'s `load_kefa_meteo_data(t_start, t_end)` → `(atmos_df, wave_df)`.
- **Key parameters**: channel selection `channels_to_plot` (subset of the 55 channels), frequency band for the overview panel (e.g., 8–20 Hz), time range for the full deployment period.
- **Outputs**: PNG spectrogram plots with meteorological time series overlay; not written to `products/`.
- **Dependencies**: `functions/plot.py`'s `Ticks_xaxis_day_or_hour`; `meteo/load_kefa_meteo.py`.

---

### Script 3: `3_plot_ppsd.py` — Probabilistic PSD (PPSD)

- **Purpose**: Compute a 2D histogram (PPSD) of PSD across all time windows and plot a frequency–dB heatmap.
- **Inputs**: `psd_daily_h5/` (Script 1 output).
- **Key parameters**: dB bin edges (e.g., −220 to −100 dB re 1 (strain-rate)²/Hz), frequency bin edges.
- **Outputs**: PPSD PNG figures.
- **Dependencies**: standard library + `numpy` + `matplotlib`. Does **not** import `functions/`.
- **Note**: This script does **not** produce any file consumed by `products/` or by the MF release.

---

### Script 4: `4_analyse_noise_trends.py` — ★ Noise Trend Analysis (Produces Interface Product)

- **Purpose**: Compute the 5th-percentile noise amplitude in the 8–20 Hz band per time window per channel across the deployment period; save as the frozen product h5 for MF interface use.
- **Inputs**: `psd_daily_h5/` (Script 1 output).
- **Key parameters**:
  - `freq_band = (8, 20)` Hz (8–20 Hz band, indices into the 8001-point frequency axis).
  - `percentile = 5` (5th percentile over time).
  - `channels_to_average_noise_level` = 28 channels corresponding to 8.3–14 km distance (sub-selection for the lower panel trend curve).
  - `output_file = 'products/noise_trends_f8-20Hz_p5.h5'`.
- **Core math**:
  1. For each time window (chunk) and channel, sum the PSD within the 8–20 Hz band indices, take square root → noise amplitude in strain-rate RMS units.
  2. Collect the amplitude time series for all chunks; compute `np.percentile(..., 5)` over the time axis → per-channel P5 amplitude.
  3. Write `noise_amplitude` array and axis labels to the output h5.
- **Outputs**: `{PKG_ROOT}/products/noise_trends_f8-20Hz_p5.h5`, containing:
  - `noise_amplitude`: shape `(2135, 55)`, float64. Each row = one time window; each column = one of the 55 channels. Values are RMS noise amplitude in the 8–20 Hz band.
  - `channels_km`: shape `(55,)`, channel distances in km.
  - `timestamps_utc_str`: shape `(2135,)`, ISO format UTC strings for each time window.
- **This file is already bundled as `products/noise_trends_f8-20Hz_p5.h5`** (967 KB). Re-running this script will overwrite it (requires `psd_daily_h5/`).
- **Dependencies**: `functions/plot.py`'s tick formatters; `meteo/load_kefa_meteo.py` (for the overlay meteorological time series in the output figure panel).

---

## III. Frozen Interface Product

`products/noise_trends_f8-20Hz_p5.h5` is the **only file from this sub-package consumed by the companion MF workflow repository `WANGYuhan2002/yuhan-mfdas` (`yuhan_dasmf_v1.0/` in `matched_filter_release.zip` on Zenodo)**. Its datasets and shapes are a frozen contract:

| Dataset | Shape | Description |
|---|---|---|
| `noise_amplitude` | `(2135, 55)` | 5th-percentile noise RMS amplitude in 8–20 Hz band per time window per channel |
| `channels_km` | `(55,)` | Channel distances in km along the fiber |
| `timestamps_utc_str` | `(2135,)` | UTC timestamps of each time window (ISO string) |

In `WANGYuhan2002/yuhan-mfdas`, `Run_4_assign_overlaps_with_kurtosis/5_Show_withMag_noise_Mc_bvalue_slip_momentum.py` reads this file from `IN_noise_meteo/noise_trends_f8-20Hz_p5.h5`. Note that the GitHub version of `WANGYuhan2002/yuhan-mfdas` only ships `load_kefa_meteo.py` in `IN_noise_meteo/`; the data files `noise_trends_f8-20Hz_p5.h5` and `kefa_*_data.csv` live here in `WANGYuhan2002/yuhan_das_spectra` (and in the Zenodo archive).

---

## IV. Important Notes

### Resolution Incompatibility with eq_spectra/

`eq_spectra/2_plot_psd_by_eq.py` requires a PPSD background computed at **10-minute chunks / 40-second nperseg** (subdirectory `chunk0:10:00_nperseg40s/`), but Script 1 in this folder produces 40-minute chunks / 200-second nperseg (subdirectory `chunk0:40:00_nperseg200s/`). These are **incompatible** — `2_plot_psd_by_eq.py` will skip the PPSD background panel if the correct resolution directory is not found. This is a known issue requiring author action (either recompute a separate 10-minute resolution PSD, or accept the missing background).

### `psd_daily_h5/` Not Bundled

Script 1 must access raw DAS data from the mounted drive to produce `psd_daily_h5/`. This intermediate directory (tens of GB) is not bundled. Anyone needing to re-run the full pipeline from scratch must have access to the DAS data and set `DAS_DATA_DIR` correctly in Script 1 (or place the data at its default `IN_DAS_data/`). Without DAS data Script 1 does not fail: it writes 60 empty daily h5 files and exits with status 0. The one-day sample `2024-04-24.h5` (69.3 MB) is in the Zenodo archive (`matched_filter_release.zip` → `check_minmag_pkg/noise_psd/example_outputs/sample_intermediate_psd_daily/2024-04-24.h5`), not in the GitHub repository (see `example_outputs/sample_intermediate_psd_daily/`).
