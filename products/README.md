# products/ — Documentation

## I. Contents and Purpose

This folder contains a single frozen product:

| File | Size | Description |
|---|---|---|
| `noise_trends_f8-20Hz_p5.h5` | ~948 KB (967448 bytes) | DAS noise amplitude time series product, produced by `noise_psd/4_analyse_noise_trends.py` |

This is the final numerical output of the `noise_psd/` noise analysis pipeline: within the 8–20 Hz frequency band, for each time step and each channel, the 5th percentile (P5) of the PSD is taken and then square-rooted to yield a "noise amplitude" proxy. It is also the **frozen interface contract product between this package and the MF release** — the MF uses it to link detection capability (minimum detectable magnitude Mc) to actual noise levels.

## II. Item-by-Item Description (Fields / Data Flow)

### Data Flow (How It Is Produced)

```
noise_psd/1_calculate_psd_h5.py  →  daily PSD h5 files (chunk=40min / nperseg=200s, 55 channels)
   │  noise_psd/4_analyse_noise_trends.py
   │    · take P5 (5th percentile) of PSD within the 8–20 Hz band
   │    · square-root to obtain amplitude proxy:  sqrt(P5 of PSD)
   ▼
products/noise_trends_f8-20Hz_p5.h5   ← this product (frozen)
```

Filename convention: `f8-20Hz` = frequency band 8–20 Hz; `p5` = 5th percentile (determined by `frequency_band=(8,20)` and `percentile_to_calculate=5` in the generating script's `Config`).

### h5 Contents (Field by Field, Verified by Inspection)

**Datasets**

| Name | Shape | dtype | Meaning |
|---|---|---|---|
| `noise_amplitude` | `(2135, 55)` | `float64` | `[time step, channel]`; each cell = `sqrt(P5 of PSD)` for that time step and channel in 8–20 Hz |
| `timestamps_utc_str` | `(2135,)` | `\|S19` | UTC start-time string for each time step; format `%Y-%m-%dT%H:%M:%S` (e.g. `2024-04-24T10:20:00` … `2024-06-22T17:00:00`) |
| `channels_km` | `(55,)` | `float64` | Channel position along fiber (km), 0.5 … 15.2 km |

**Attributes**

| Name | Value |
|---|---|
| `amplitude_unit` | `(strain/s)/sqrt(Hz)` |
| `amplitude_description` | `sqrt(P5 of PSD)` |
| `frequency_band_hz` | `(8, 20)` |
| `time_step_duration` | `0:40:00` (40 minutes per time step, consistent with the PSD computation chunk) |
| `description` | `Noise amplitude trend analysis results.` |
| `creation_date_utc` | Creation timestamp (e.g. `2025-07-06T19:37:52...`) |

Reading convention: `timestamps_utc_str` is stored as byte strings; decode with `.decode('utf-8')` then parse with `%Y-%m-%dT%H:%M:%S`. The two dimensions of `noise_amplitude` correspond to `timestamps_utc_str` (rows / time steps) and `channels_km` (columns / channels).

## III. Integration with the Rest of This Package and with MF

- **Within this package**: this product is written by `noise_psd/4_analyse_noise_trends.py`'s `save_results_to_h5()` (the same script also produces heatmaps and per-channel subplots from it). To recompute: first run `noise_psd/1_` to generate daily PSD h5 files, then run `noise_psd/4_`.
- **Frozen interface contract with the MF release (`WANGYuhan2002/yuhan-mfdas`, folder `yuhan_dasmf_v1.0/` in `matched_filter_release.zip`)**: `Run_4_assign_overlaps_with_kurtosis/5_Show_withMag_noise_Mc_bvalue_slip_momentum.py` in `WANGYuhan2002/yuhan-mfdas` directly reads the three datasets (`noise_amplitude`, `channels_km`, `timestamps_utc_str`) from this product, and averages over the **9–14 km channels** (`channels_to_avg = [9,10,11,12,13,14]`) along the channel dimension (`noise_matrix[:, indices].mean(axis=1)`) to obtain a noise amplitude time series for plotting the magnitude–Mc–b-value–slip–seismic moment combined figure.
- In the Zenodo archive (`matched_filter_release.zip`), this product is also copied into `yuhan_dasmf_v1.0/IN_noise_meteo/noise_trends_f8-20Hz_p5.h5`; in the GitHub repositories, `WANGYuhan2002/yuhan-mfdas` only ships `load_kefa_meteo.py` in `IN_noise_meteo/`, while the data files `noise_trends_f8-20Hz_p5.h5` and `kefa_*_data.csv` live here in `WANGYuhan2002/yuhan_das_spectra` (and in the Zenodo archive).
- Therefore the three dataset names, shape semantics (`[time step, channel]`), time string format, and `channels_km` values form the **stable interface with MF** and must not be changed; otherwise the MF-side reading and 9–14 km averaging logic will break.
