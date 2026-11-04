# Intermediate-data sample: single-day PSD h5

> **GitHub version:** the h5 file described here is **not included** in this repository (only this note is). The 69.3 MB daily PSD file `2024-04-24.h5` (72,712,352 bytes) is in the Zenodo archive [10.5281/zenodo.20715241](https://doi.org/10.5281/zenodo.20715241) (`matched_filter_release.zip` → `check_minmag_pkg/noise_psd/example_outputs/sample_intermediate_psd_daily/2024-04-24.h5`). Copy it here, or into `noise_psd/chunk0:40:00_nperseg200s/psd_daily_h5/` to use it as input for scripts 2_/3_/4_.

`2024-04-24.h5` (~70 MB) is a **sample intermediate product** of
`noise_psd/1_calculate_psd_h5.py` — the per-channel Welch PSD for a single day.

- **Produced by**: `1_calculate_psd_h5.py` (DAS → one h5 per day).
- **Structure**: `frequency` (length 8001, 0–40 Hz, Δf = 0.005 Hz) + one group per channel
  `ch_X.Xkm/`, whose dataset names are the start timestamp of each 40-min chunk
  (`YYYYmmddTHHMMSS`) and whose values are that chunk's PSD `Pxx` for that channel
  (unit `(strain/s)²/Hz`).
- **Consumed by**: `2_plot_results_with_wind_wave.py` / `3_plot_ppsd.py` / `4_analyse_noise_trends.py`.

> Only **one day** is included here as a sample. The full set (2024-04-24 → 06-22, 60 days,
> ~6.9 GB) is not bundled and can be recomputed by `1_calculate_psd_h5.py` on external
> continuous DAS.
