# Intermediate-data sample: per-event PSD h5 for a single channel

> **GitHub version:** the h5 file described here (`psd_data_channel_11.0km.h5`) is **not included** in this repository (its HDF5 metadata holds a machine-specific path); recompute it with `eq_spectra/1_calculate_psd_by_eq.py` as described below. Neither of the two sample intermediate h5 files of this package is shipped in the GitHub version; the 69.3 MB daily PSD sample of `noise_psd/` is in the Zenodo archive [10.5281/zenodo.20715241](https://doi.org/10.5281/zenodo.20715241) (`matched_filter_release.zip` → `check_minmag_pkg/noise_psd/example_outputs/sample_intermediate_psd_daily/2024-04-24.h5`).

`psd_data_channel_11.0km.h5` (~980 KB) is a **sample intermediate product** of
`eq_spectra/1_calculate_psd_by_eq.py` — the PSD of all earthquake events on one channel
(11.0 km, template_8).

- **Produced by**: `1_calculate_psd_by_eq.py` (reads the MF event catalogue
  `Events_classified_Real_Mag.csv` + external DAS, computing a Welch PSD per event).
- **Structure**: one dataset per earthquake event; each dataset's attrs include that event's
  `magnitude` etc.; files are split by channel (`psd_data_channel_<km>.h5`).
- **Consumed by**: `2_plot_psd_by_eq.py` (overlays magnitude-binned median earthquake spectra on
  the noise background).

> Only **one channel** is included here as a sample. The full set (each template × each channel,
> ~1.7 GB) is not bundled and can be recomputed by `1_calculate_psd_by_eq.py` on external DAS +
> the MF event catalogue.
