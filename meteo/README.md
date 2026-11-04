# meteo/ — Documentation

## I. Contents and Purpose

This folder provides ERA5 (ECMWF Reanalysis v5) meteorological data for the Kefalonia DAS deployment area (April–June 2024) and a unified loading interface for those data.

> **Data source and attribution:** the two CSVs are point extractions from the ERA5 reanalysis (Hersbach et al., 2020, https://doi.org/10.1002/qj.3803), produced by ECMWF for the Copernicus Climate Change Service (C3S) and obtained from the Copernicus Climate Data Store. Contains modified Copernicus Climate Change Service information. Neither the European Commission nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains. See also `../NOTICE`.

| File/Folder | Purpose |
|---|---|
| `load_kefa_meteo.py` | **Frozen interface**: `load_kefa_meteo_data(t_start, t_end)` → `(atmos_df, wave_df)` |
| `kefa_atmos_data.csv` | Atmospheric data (wind speed, pressure, sea-surface temperature, wind direction) |
| `kefa_wave_data.csv` | Ocean-wave data (significant wave height SWH, mean wave period) |
| `acquire_era5/0_download_from_api.py` | CDS API download script |
| `acquire_era5/3_to_csv.py` | ERA5 NetCDF → CSV conversion script |

The companion MF workflow repository (`WANGYuhan2002/yuhan-mfdas`, folder `yuhan_dasmf_v1.0/` in `matched_filter_release.zip` on Zenodo) consumes only two columns from this folder: `atmos_df['wind_speed']` and `wave_df['swh']`.

---

## II. File Descriptions

### `load_kefa_meteo.py` — Frozen Loading Interface

**Entry point**: `load_kefa_meteo_data(t_start, t_end, data_dir=None)`

- **Arguments**:
  - `t_start`, `t_end`: `datetime.datetime` objects (naive UTC or localized), the requested time range.
  - `data_dir`: directory containing the CSV files. Defaults to `os.path.dirname(__file__)` (i.e., the `meteo/` folder itself), so no hard-coded absolute path is needed — works regardless of where the package is placed.
- **Return value**: `(atmos_df, wave_df)`, where:
  - `atmos_df`: `pandas.DataFrame`, columns `['valid_time', 'sst', 'msl', 'wind_speed', 'wind_direction']`; `valid_time` is a UTC `datetime64` index; rows are hourly ERA5 timestamps within `[t_start, t_end]`.
  - `wave_df`: `pandas.DataFrame`, columns `['valid_time', 'swh', 'mwp']`; `valid_time` is a UTC `datetime64` index; rows are hourly ERA5 timestamps within `[t_start, t_end]`.

**Frozen contract**: the function name (`load_kefa_meteo_data`) and positional argument order (`t_start`, `t_end`) must not be changed. The MF release uses this function by importing at a hardcoded path; changing the interface would break three downstream call sites.

**Portability fix (de-hardcoding)**: the original version had `data_dir` hard-coded to an absolute path. The current version uses `data_dir=None` defaulting to `os.path.dirname(__file__)`, making the package self-contained.

---

### `kefa_atmos_data.csv` — Atmospheric Data

- **Source**: ERA5 hourly surface variables, area = [38.6, 20.4, 38.1, 21.0] (lat_max, lon_min, lat_min, lon_max), nearest-grid-point extraction at lat=38.5, lon=20.5.
- **Columns**: `valid_time` (UTC hourly timestamps), `sst` (sea-surface temperature, K), `msl` (mean-sea-level pressure, Pa), `wind_speed` (m/s, derived from u10/v10 components), `wind_direction` (degrees from North).
- **Record count**: approximately 1423 records (covering April–June 2024, hourly).
- **MF usage**: only `wind_speed` is used.

---

### `kefa_wave_data.csv` — Wave Data

- **Source**: ERA5 hourly ocean-wave variables, same area/grid point as atmospheric data.
- **Columns**: `valid_time` (UTC hourly timestamps), `swh` (significant wave height, m), `mwp` (mean wave period, s).
- **MF usage**: only `swh` is used.

---

### `acquire_era5/0_download_from_api.py` — CDS API Download Script

- Requests ERA5 data via the `cdsapi` library (Copernicus Climate Data Store). As shipped it requests the *monthly-mean* products (`reanalysis-era5-single-levels-monthly-means`, `reanalysis-era5-wave-monthly-means`) and writes `kefalonia_atmospheric_data.nc` / `kefalonia_wave_data.nc`; `3_to_csv.py` instead reads hourly files `kefa_wind/data_stream-*.nc`, and the shipped CSVs are hourly. The two scripts do not chain as shipped (GitHub version note; code unchanged) — check the request before re-downloading.
- Area parameter: `area=[38.6, 20.4, 38.1, 21]` (lat_max, lon_min, lat_min, lon_max).
- **Not needed for re-running analysis**: the CSVs are already bundled in this folder. Only needed if ERA5 data needs to be refreshed or extended.
- Requires a CDS API account and `~/.cdsapirc` credentials file.

---

### `acquire_era5/3_to_csv.py` — NetCDF → CSV Conversion

- Reads downloaded ERA5 NetCDF files; extracts the nearest grid point (lat=38.5, lon=20.5); saves to `kefa_atmos_data.csv` and `kefa_wave_data.csv`.
- **Not needed for re-running analysis**: CSVs are already bundled. Only needed if re-downloading or extending ERA5 data.

---

## III. Interface with the MF Release (`WANGYuhan2002/yuhan-mfdas`)

Three planned call sites in `WANGYuhan2002/yuhan-mfdas` (`yuhan_dasmf_v1.0/` in `matched_filter_release.zip` on Zenodo):

| Call site | Script | Used columns |
|---|---|---|
| Output_1 plotting | `Output_1_SCCs_and_MAD/.../999_plot*.py` | `wind_speed`, `swh` |
| Run_4 plotting | `Run_4_assign_overlaps_with_kurtosis/5_Show*.py` | `wind_speed`, `swh` |
| Run_4 export | `Run_4_assign_overlaps_with_kurtosis/5_extract*.py` | `wind_speed`, `swh` |

> **GitHub version note**: In the shipped `WANGYuhan2002/yuhan-mfdas` (`yuhan_dasmf_v1.0/`) code, only `Output_1_SCCs_and_MAD/MAD_arrays/999_plot_MADarr_with_meteo_and_stalta.py` imports this loader, from its own copy in `IN_noise_meteo/` (`sys.path.append(f"{PROJECT_ROOT}/IN_noise_meteo")` + `from load_kefa_meteo import load_kefa_meteo_data`); the two `Run_4_assign_overlaps_with_kurtosis/5_*` scripts read only the noise h5 (`IN_noise_meteo/noise_trends_f8-20Hz_p5.h5`). The GitHub version of `WANGYuhan2002/yuhan-mfdas` only ships `load_kefa_meteo.py` in `IN_noise_meteo/`; the data files `noise_trends_f8-20Hz_p5.h5`, `kefa_atmos_data.csv`, and `kefa_wave_data.csv` live here in `WANGYuhan2002/yuhan_das_spectra` (`products/` and `meteo/`) and in the Zenodo archive `matched_filter_release.zip`.
