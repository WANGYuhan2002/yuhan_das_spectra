# functions/ — Documentation

## I. Contents and Purpose

This folder contains only two shared modules — all `functions` dependencies actually `import`-ed by the scripts selected for this package:

| Module | Purpose | Called By |
|---|---|---|
| `load_das_h5_PE2024.py` | DAS HDF5 custom loader (core function `load_das_custom`); reads continuous DAS data by time window | `noise_psd/1_calculate_psd_h5.py`, `eq_spectra/1_calculate_psd_by_eq.py`, `eq_spectra/1_calculate_psd_by_eq_Parallel.py`, `noise_psd/2_`, `noise_psd/4_` (the latter two import it but actual DAS reading is done in `1_`) |
| `plot.py` | Plotting utilities (`Ticks_xaxis_day_or_hour` / `Ticks_xaxis_minute_or_second` and other time-axis tick formatters) | `noise_psd/2_plot_results_with_wind_wave.py`, `noise_psd/4_analyse_noise_trends.py` |

**Why only these two**: Verified by `grep` against the `functions` symbols actually `import`-ed by the 7 selected scripts in this package — only `load_das_h5_PE2024` and `plot` are used, so only these two modules are included. Other modules from the original `check_minmag` project's `functions/` directory (`compute_correlation`, `trigger`, `classify`, `extract`, `extract_station_info`, `cluster`, `readers`, `das_mffi_io`, `append_cut_convert`, `stalta`, etc., and the entire `pydas_readers/` directory) are not imported by any of the selected scripts and have been excluded.

Verified by testing: `load_das_h5_PE2024.py` is self-contained — the entire file only imports from the standard library (`os`/`sys`/`glob`/`datetime`, etc.) + `h5py` + `numpy`, with **no dependency** on `pydas_readers/` or any other `functions` module. Tested: `import load_das_h5_PE2024` succeeds on its own and provides `load_das_custom` with no `pydas_readers` dependency (the planning document's conjecture of "depends on pydas_readers" is disproved by direct code reading). `plot.py` only imports `matplotlib` + `numpy`. Neither requires any transitive dependencies.

## II. Script-by-Script Description (I/O and Data Flow)

### `load_das_h5_PE2024.py`

- **Origin**: Derived from Pascal Edme (ETH Zürich)'s tool for reading Silixa iDASv2 / downsampled DAS HDF5 files, 2023 version.
- **Main functions**:
  - `load_headers_only(file)`: Reads only the header of a single h5 file (without reading large data blocks); returns a dict containing `fs` (sampling rate, derived from `1000/step`), `dx` (channel spacing), `nchan`, `npts`, `t0`/`t1` (inferred from filename using `%Y%m%d_%H%M%S.%f.h5` format), `d0`/`d1` (fiber distance start/end), `amp_scaling` (amplitude calibration, default 1), etc.
  - `load_file(file)`: Reads the `Traces` dataset from a single h5 file, divides by `amp_scaling` to convert to strain-rate units `(m/m)/s`; returns `(data, headers, axis)` where `data` has shape `[samples, channels]`.
  - **`load_das_custom(t_start, t_end, d_start=0, d_end=0, input_dir, nth_channel=1, ...)`** — the actual entry point used by scripts in this package:
    - Input: start/end `datetime`, optional fiber distance range `d_start`/`d_end`, data directory `input_dir`, optional channel stride `nth_channel`.
    - Behavior: locates files in `input_dir` by filename timestamp, stitches across files and days, extracts only the required time (and distance) interval, avoiding full block reads.
    - Output: `(data, headers, axis)` where `data` is a strain-rate array of shape `[samples, channels]`; `axis` contains `tt` (seconds) / `date_times` (absolute times) / `dd` (channel distances).
- **Call sites**:
  - `noise_psd/1_calculate_psd_h5.py:163` — `load_das_custom(t_start, t_end, input_dir=DAS_DATA_DIR)` to read continuous DAS and compute Welch PSD.
  - `eq_spectra/1_calculate_psd_by_eq.py:129` (and `_Parallel.py`) — same `load_das_custom` to read DAS around trigger times for per-event PSD.

### `plot.py`

- **Contents**: Matplotlib plotting utilities. The two time-axis tick formatters actually used by this package are:
  - `Ticks_xaxis_day_or_hour(ax, Day_tick, Hour_tick, rotation, font)`: Sets major/minor ticks and date format at "day" or "hour" intervals.
  - `Ticks_xaxis_minute_or_second(ax, Minute_tick, Second_tick, rotation, font)`: Sets ticks at "minute" or "second" intervals.
  - Also contains `waterfall_date`, `generate_gs_grid`, and other helpers (not directly called by current scripts in this package, but included with the module).
- **Call sites**: `noise_psd/2_` and `noise_psd/4_` use `from plot import Ticks_xaxis_day_or_hour, Ticks_xaxis_minute_or_second` during figure generation.

## III. Interface with Other Parts of This Package and MF

- This folder is the **low-level dependency for computation and plotting**, imported by scripts in `noise_psd/` and `eq_spectra/` after adding the path via `sys.path` (each script derives the path from `_PKG_ROOT/functions` relative to `__file__`, without absolute paths).
- This folder does **not** directly output any interface or product to the companion MF release (`WANGYuhan2002/yuhan-mfdas`, `yuhan_dasmf_v1.0/` in `matched_filter_release.zip` on Zenodo) — what MF consumes is `meteo/` (the meteorological loading module + CSVs) and `products/` (the frozen noise product h5), not `functions/`. `functions/` is only used by scripts in this package when "recomputing PSD / earthquake spectra from scratch"; MF uses the frozen products directly and does not touch this folder.
