# eq_spectra/ — Documentation

## I. Overview and Scientific Purpose

This sub-package computes the **per-event Power Spectral Density (PSD) spectra** of DAS waveforms for earthquakes detected by the matched-filter method, then plots them against a noise PSD background. The goal is to assess spectral content and signal-to-noise ratios of detected events — an important step in validating the minimum magnitude of completeness (Mc) for the MF catalog.

**Two parallel analysis lines** are provided:

| Line | Scripts | Channels | Event List (Default) |
|---|---|---|---|
| Main (serial) | `1_calculate_psd_by_eq.py` + `2_plot_psd_by_eq.py` | 13 channels (8.3–14 km) | `eq_list_MF_v3_40_RAWCCsMag_5-20Hz/` |
| Parallel variant | `1_calculate_psd_by_eq_Parallel.py` + `2_plot_psd_by_eq.py` | 28 channels (1–14 km) | Same |

Both step-1 scripts write to the same per-event h5 format; `2_plot_psd_by_eq.py` reads whichever step-1 output is present.

---

## II. Event Lists

Two event lists are bundled in `eq_spectra/`:

| Directory | Event count | Magnitude source | Frequency band |
|---|---|---|---|
| `eq_list_MF_v3_40_RAWCCsMag_5-20Hz/` | ~13,492 events | RAWCCs-based quasi-magnitude, calibrated in 5–20 Hz | **Default** |
| `eq_list_MF_v3_40_ampMag_10-35Hz/` | ~13,492 events | Peak-amplitude-based quasi-magnitude, calibrated in 10–35 Hz | Alternate |

Both event lists are tab-separated text files. The only columns read by the scripts are:
- `Rt_trig`: trigger time (string, format `YYYYMMDDTHHMMSS`)
- `Real_Mag`: magnitude (float)
- `Template`: template name (string, e.g. `template_10_80Hz`)

The default (`RAWCCsMag`) is used unless the user changes `eq_list_dir` in the script.

---

## III. Script-by-Script Description

### Script 1a: `1_calculate_psd_by_eq.py` — Per-Event PSD (Serial, 13 Channels)

- **Purpose**: For each event in the event list, read a 4-second DAS window centered on the trigger time, compute Welch PSD, and save as a per-event dataset in an h5 file. Serial (single-process).
- **Inputs**:
  - DAS data: `DAS_DATA_DIR` (package-relative default `IN_DAS_data/`, same as `noise_psd/1_`; place or symlink the 80 Hz archive there, or change the line).
  - Event list: `eq_list_MF_v3_40_RAWCCsMag_5-20Hz/` (default).
  - Key parameters: `fs=80`, `nperseg=4×80=320` samples (4-second Welch segment), `das_dist_start=8300` m, `das_dist_end=14000` m (13 channels at 8.3–14 km), time window = `[t_trig − 1.5 s, t_trig + 2.5 s]`.
- **Outputs**: per-template h5 file `eq_spectra/psd_by_eq_h5/<template_name>_13ch.h5`, containing:
  - `frequency`: shape `(161,)`, float64, frequency axis 0–40 Hz (Δf = 0.25 Hz = 80/320).
  - Per-event dataset `event_YYYYMMDDTHHMMSS`: shape `(161, 13)` float64, PSD per channel; attrs: `magnitude` (float), `trigger_time` (string).
- **Dependencies**: `functions/load_das_h5_PE2024.py`'s `load_das_custom`.

---

### Script 1b: `1_calculate_psd_by_eq_Parallel.py` — Per-Event PSD (Parallel, 28 Channels)

- **Purpose**: Same as Script 1a, but uses multiprocessing for speed and extends the channel range from 8.3–14 km (13 channels) to 1–14 km (28 channels).
- **Key differences from 1a**:
  - `das_dist_start=1000` m, `das_dist_end=14000` m → 28 channels.
  - Output filename suffix `_28ch.h5`.
  - Uses `multiprocessing.Pool` for parallel processing.
- **Inputs/Outputs**: same structure as 1a, but per-event dataset shape is `(161, 28)`.

---

### Script 2: `2_plot_psd_by_eq.py` — PSD Plotting with PPSD Background

- **Purpose**: Read per-event PSD h5 files (from Script 1a or 1b), bin events by magnitude, and plot median PSD curves per magnitude bin overlaid on a noise PPSD background.
- **Inputs**:
  - Per-event PSD h5 (from Script 1a or 1b).
  - **PPSD background**: a PPSD h5 computed by `noise_psd/3_plot_ppsd.py` at **10-minute chunks / 40-second nperseg** resolution → subdirectory `chunk0:10:00_nperseg40s/`.
- **Outputs**: PNG figures with PSD curves per template per magnitude bin overlaid on PPSD background.
- **Critical note**: See the incompatibility warning below.

---

## IV. Important Incompatibility: PPSD Background Resolution Mismatch

`2_plot_psd_by_eq.py` looks for a PPSD background in:
```
noise_psd/psd_daily_h5/chunk0:10:00_nperseg40s/
```
But the standard `noise_psd/` pipeline (Script 1) produces:
```
noise_psd/psd_daily_h5/chunk0:40:00_nperseg200s/
```
These are **incompatible** — the required 10-minute/40-second resolution directory **does not exist** if only the standard noise_psd pipeline has been run.

**Behavior**: `2_plot_psd_by_eq.py` will **silently skip the PPSD background panel** if the directory is not found; it still plots the per-event PSD curves without a background.

**Resolution**: The author needs to either:
1. Re-run `noise_psd/1_calculate_psd_h5.py` with `chunk_sec=600` (10 min) and `nperseg_sec=40` to produce the required resolution, OR
2. Accept the missing PPSD background and proceed with PSD-curves-only plots.

---

## V. Data Format: Per-Event PSD h5

```
eq_spectra/psd_by_eq_h5/<template_name>_<Nch>.h5
├── frequency                  (161,)  float64   0–40 Hz, Δf=0.25 Hz
├── event_20240425T102312      (161, N_channels)  float64   PSD per channel
│     attrs: magnitude=2.4, trigger_time="20240425T102312"
├── event_20240425T103541      ...
└── ...
```

- The dataset name encodes the trigger time: `event_YYYYMMDDTHHMMSS` (seconds precision).
- Shape: `(161, 13)` for 13-channel output; `(161, 28)` for 28-channel output.
- Per-event PSD h5 files are **not bundled** in this package; they must be recomputed from DAS data using Script 1a or 1b.

---

## VI. Execution Dependencies

```
DAS mounted drive  +  eq_list_MF_v3_40_RAWCCsMag_5-20Hz/
         │
         ▼
[Script 1a serial] or [Script 1b parallel]   → psd_by_eq_h5/<template>_Nch.h5
                                                      │
         noise_psd/psd_daily_h5/                      ▼
         chunk0:10:00_nperseg40s/  ────────→  [Script 2: plot]  → PNGs
         (PPSD background, OPTIONAL)
         (if absent: PPSD panel skipped)
```

Only DAS data is required to run Steps 1→2. The PPSD background for Step 2 requires an additional re-run of `noise_psd/1_calculate_psd_h5.py` at a different resolution.
