# noise_psd / example_outputs — Representative result figures

> These are **example figures** produced by the `noise_psd/` pipeline (called "Pipeline A",
> the noise-PSD main line, in the package README), included only so you can see at a glance
> "what the code actually produces".
> **Not the full set**: only a few representative channels' PNGs are kept here (~21 MB).
> The full set can be recomputed from external large DAS by the `1_/2_/3_/4_*.py` scripts one
> level up.
> This folder holds **no intermediate big data** (`psd_daily_h5/` per-day PSD, tens of GB) — figures only.
> Example figures are read-only copies from the source project (source unmodified):
> `check_minmag_validity/Run_4_amp_spec_final/chunk0:40:00_nperseg200s/`
> Written: 2026-06-15.

---

## Data flow at a glance (framework before looking at the figures)

```
External large DAS  ──[1_calculate_psd_h5.py]──>  psd_daily_h5/  (intermediate, not bundled, not in this folder)
                                              │
                 ┌────────────────────────────┼────────────────────────────┐
                 ▼                            ▼                            ▼
        [2_plot_results_with_      [3_plot_ppsd.py]            [4_analyse_noise_trends.py]
         wind_wave.py]              PPSD probabilistic spectra  noise-trend heatmap + product h5
         spectrogram + meteo       (per-channel + summary)     (★ produces noise_trends h5)
```

`1_` is the only data-producing script and **draws nothing**, so this folder has no figure for it;
`2_/3_/4_` each independently consume `psd_daily_h5/`, and the figures here are the outputs of those three.

---

## Which script each figure here comes from

### From `4_analyse_noise_trends.py` (noise-trend analysis; ★ visualises the in-package product noise_trends h5)

| File | Description |
|---|---|
| `noise_heatmap_f8-20Hz_p5.png` | **Two-panel figure.** Top: channel (0.5–15.2 km) × time heatmap of the "noise amplitude proxy" (P5 percentile over the 8–20 Hz band, then square-rooted); bottom: noise-amplitude time series averaged over the 9–14 km segment (`channels_to_average=[8.3…14]`) along channels (log y-axis). **This frequency-band / percentile combination is exactly the visualisation corresponding to the product `products/noise_trends_f8-20Hz_p5.h5` referenced by the MF release.** |
| `noise_heatmap_f0.02-1Hz_p5.png` | Same two-panel layout, but for the low band 0.02–1 Hz (ocean-wave secondary microseism band), for comparison with the high band. |
| `noise_subplots_f8-20Hz_p5.png` | **Per-channel subplots.** Noise-amplitude time series of all 55 channels arranged in a grid, to inspect trends channel by channel. The file is large (~11 MB), so only the high-band one is kept as a representative; the low-band `noise_subplots_f0.02-1Hz_p5.png` can likewise be recomputed by `4_`. |

### From `2_plot_results_with_wind_wave.py` (spectrogram + amplitude time series + overlaid wind/wave meteo)

| File | Description |
|---|---|
| `sg_meteo_channel_9.0.png` | Single-channel **3-in-1 integrated figure**: ① time–frequency spectrogram (ASD, LogNorm); ② band-limited (0.02–1 Hz, 4–35 Hz) ASD time series + significant wave height SWH on the right axis; ③ wind speed + mean wave period MWP on the right axis. Meteo data come from the in-package `meteo/`. Three representative channels (9.0 / 11.0 / 13.0 km) from the 9–14 km analysis segment are shown. |
| `sg_meteo_channel_11.0.png` | Same, channel 11.0 km. |
| `sg_meteo_channel_13.0.png` | Same, channel 13.0 km. |

> The full set is one figure per channel for all 55 channels (~2 MB each in the source); only 3 representatives are kept here.

### From `3_plot_ppsd.py` (Probabilistic Power Spectral Density, McNamara–Buland style)

| File | Description |
|---|---|
| `ppsd_ch5p0km_20240424-20240622.png`  | Single-channel **PPSD**: a 2D histogram in the "frequency × power (dB)" plane over all PSD spectra within the channel's full time window, normalised to a probability density per frequency. Channel 5.0 km. |
| `ppsd_ch9p0km_20240424-20240622.png`  | Same, channel 9.0 km. |
| `ppsd_ch11p0km_20240424-20240622.png` | Same, channel 11.0 km. |
| `ppsd_ch13p0km_20240424-20240622.png` | Same, channel 13.0 km. |
| `ppsd_summary_50percentile_20240424-20240622.png` | **Summary figure**: per-channel P50 median percentile of PSD over time, stacked into a [channel × frequency] matrix heatmap — the median noise level of all 55 channels in one figure. |

> The full set is one PPSD per channel for all 55 channels (~0.47 MB each) + 1 summary; only 4 + the summary are kept here.

---

## How to recompute the full set

With external large DAS mounted, under `noise_psd/` run in order: first `1_calculate_psd_h5.py`
to generate `psd_daily_h5/` (intermediate products), then `2_/3_/4_*.py` each to produce the
figures for all channels (and `4_`'s product h5). Parameters (channel list, frequency band, time
window, etc.) are in the Config block at the top of each script and in the pipeline's `../README.md`.
