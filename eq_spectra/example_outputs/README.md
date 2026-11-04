# eq_spectra / example_outputs — Representative result figures

> These are **example figures** produced by the `eq_spectra/` pipeline (called "Pipeline B",
> earthquake spectra by magnitude, in the package README), included only so you can see at a
> glance "what the code actually produces".
> **Not the full set**: only 9 representative channels' per-event spectral PNGs are kept here (~11 MB).
> The full set can be recomputed from external large DAS by `1_calculate_psd_by_eq*.py` +
> `2_plot_psd_by_eq.py` one level up.
> This folder holds **no intermediate big data** (per-event `psd_data_channel_*.h5`, continuous
> PSD daily h5) — figures only.
> Written: 2026-06-15.

---

## These figures are the output of `2_plot_psd_by_eq.py`

`2_plot_psd_by_eq.py` draws one combined spectral figure **per channel**, overlaying two kinds of data:

1. **Noise background (PPSD)**: the probability-density distribution of all continuous PSD over a
   long time window for that channel (blue 2D histogram, frequency × power), representing the
   ambient noise field.
2. **Earthquake spectra by magnitude (foreground spectral lines)**: per-event PSD for the same
   channel, binned by magnitude (M 0.0–0.2, …, 4.5–5.0); each bin's **median spectrum** is drawn
   as a coloured line (colour graded by magnitude via the plasma colormap), annotated with the
   magnitude range and the bin's event count.

So each figure answers: **"on this channel, how far above the ambient noise floor do earthquake
spectra of different magnitudes rise"** — the larger the magnitude, the higher the spectral line
and the more it emerges from the noise background. The red semi-transparent box is the
"amplitude extraction band" (8–20 Hz), used to take amplitudes from the spectra downstream.

> Note: the `M-1p0_5p0` in each figure's title/filename is the **overall magnitude selection
> window (−1.0, 5.0)**, not a single magnitude; each figure actually overlays **several
> median spectral lines of different magnitudes** (see the in-figure legend). The PNGs are
> distinguished by **channel position** (8.3–14.0 km along the cable), not by magnitude.

---

## Which variant they come from

The package defaults to `template_0423_80Hz` + `eq_list_MF_v3_40_RAWCCsMag_5-20Hz`
(see `self.Template` and `self.MagList_subdir` in the Config at the top of `2_plot_psd_by_eq.py`).
The figures here are taken from the source project's **main variant** (RAWCCs magnitude,
consistent with the package default; read-only copy, source unmodified):

`check_minmag_validity/Run_5_eq_spec_by_magnitude/`
`eq_list_MF_v3_40_RAWCCsMag_5-20Hz/template_0423_80Hz/psd_plots_per_channel/`

> The main variant covers 13 channels over 8.3–14.0 km with consistent quality — representative
> enough — so the `Run_5_eq_spec_by_magnitude_more_channels/` variant (more channels, 56 figures)
> was not taken.

---

## File list in this folder (9 representative channels, evenly sampled along the cable)

| File | Channel position |
|---|---|
| `psd_analysis_ch8p3km_20240424-20240622_M-1p0_5p0.png`  | 8.3 km |
| `psd_analysis_ch8p6km_20240424-20240622_M-1p0_5p0.png`  | 8.6 km |
| `psd_analysis_ch9p0km_20240424-20240622_M-1p0_5p0.png`  | 9.0 km |
| `psd_analysis_ch10p0km_20240424-20240622_M-1p0_5p0.png` | 10.0 km |
| `psd_analysis_ch11p0km_20240424-20240622_M-1p0_5p0.png` | 11.0 km |
| `psd_analysis_ch11p5km_20240424-20240622_M-1p0_5p0.png` | 11.5 km |
| `psd_analysis_ch12p0km_20240424-20240622_M-1p0_5p0.png` | 12.0 km |
| `psd_analysis_ch13p0km_20240424-20240622_M-1p0_5p0.png` | 13.0 km |
| `psd_analysis_ch14p0km_20240424-20240622_M-1p0_5p0.png` | 14.0 km |

---

## How to recompute the full set

With external large DAS mounted, under `eq_spectra/`: first run `1_calculate_psd_by_eq.py`
(or the parallel version `1_calculate_psd_by_eq_Parallel.py`) to compute per-event PSD,
producing `eq_list/<MagList>/<Template>/psd_data_channel_*.h5` (intermediate, not bundled);
the continuous PSD daily h5 needed for the noise background is recomputed by
`noise_psd/1_calculate_psd_h5.py`; once both are ready, run `2_plot_psd_by_eq.py` to produce
figures for all channels in the `CHANNELS_FOR_PPSD` list. Parameters are in the Config at the
top of `2_plot_psd_by_eq.py` and in the pipeline's `../README.md`.
