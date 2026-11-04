# -*- coding: utf-8 -*-
"""
calculate_psd_v3.py

VERSION: v3.0
DESCRIPTION:
This script performs time-frequency analysis on DAS data by calculating Power
Spectral Density (PSD) using Welch's method. This version implements a hybrid
storage solution, saving the results for each day into a separate, date-named
HDF5 file.

Key Features:
- Adopts a hybrid storage model: one HDF5 file per day of data.
- Centralized parameter management using a Config class.
- The frequency vector is stored once per daily HDF5 file.
- Calculation parameters are saved to a 'config.json' file for reproducibility.
- Automatically skips calculations for which data already exists in the correct
  daily HDF5 file, allowing for resumable processing.
"""
import numpy as np
import scipy.signal as signal
from datetime import datetime, timedelta
import os
import sys
import gc
import h5py
import json
from tqdm import tqdm

# --- [去硬编码] 包根推导 + functions 路径 (相对 __file__, 不用绝对路径/环境变量) ---
# 包结构: check_minmag_pkg/{functions, noise_psd, ...}; 本文件在 noise_psd/ 下。
_PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(_PKG_ROOT, "functions"))
import load_das_h5_PE2024 as load_das_h5_PE2024

# --- [显式可见] 外部大 DAS 数据目录 (PSD 重算才需要; 换机改这里) ---
# 原脚本把它埋在 Config.input_dir; 这里提到 Config 之外做成显式顶层变量。
DAS_DATA_DIR = os.path.join(_PKG_ROOT, 'IN_DAS_data') + os.sep  # [GitHub 版] 默认 = 包内 IN_DAS_data/（放入或软链接 80 Hz 归档 DATA_DEC_80Hz_16m_SORTED/）；换机改这一行

# ####################################################################################
# # --- 1. CONFIGURATION CLASS ---
# ####################################################################################

class Config:
    """
    A class to hold all script parameters and manage output paths.
    """
    def __init__(self):
        # --- Input and Output Base Paths ---
        # [去硬编码] 输出基目录改为包内 noise_psd/ 自身目录 (相对包根);
        #   原: f"{HOME_dir}Run_4_amp_spec_final/"
        # [去硬编码] input_dir 取顶层显式变量 DAS_DATA_DIR (换机改文件顶部那一行)。
        self.output_dir_base = os.path.join(_PKG_ROOT, "noise_psd")
        self.input_dir = DAS_DATA_DIR

        # --- Time Range and Chunking Settings ---
        # self.t_all_start = datetime(2024, 5, 2, 12, 0, 0)
        # self.t_all_end = datetime(2024, 5, 5, 12, 0, 0)
        
        self.t_all_start = datetime.strptime('2024-04-24 10:20:00', '%Y-%m-%d %H:%M:%S')
        self.t_all_end = datetime.strptime('2024-06-22 17:10:00', '%Y-%m-%d %H:%M:%S')
        

        # --- Channel Selection ---
        # self.dist_channels_km = np.array([2, 4, 6, 8, 9, 10, 11, 12, 13, 14])
        # self.dist_channels_km = np.array([8.3, 8.6, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14.0])
        self.dist_channels_km = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8, 8.3, 8.6,
                                          9, 9.2, 9.4, 9.6, 9.8, 10, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 12.2, 12.4, 12.6, 12.8,
                                          13, 13.2, 13.4, 13.6, 13.8, 14,
                                          14.2, 14.4, 14.6, 14.8, 15, 15.2])
                                          
        
        
        self.dx = 16.0 * 1.020951986312866
        self.channels_to_process_idx = [int(d * 1e3 / self.dx) for d in self.dist_channels_km]
        self.sampling_rate_hz = 80.0
        
        
        #? --- Core Calculation Parameters for Welch's Method ---
        self.chunk_duration = timedelta(minutes=40)  # 40 minutes in hours
        self.nperseg_seconds = 200.0

        # --- Automatically Generated Paths ---
        params_str = f"chunk{self.chunk_duration}_nperseg{int(self.nperseg_seconds)}s"
        self.run_output_dir = os.path.join(self.output_dir_base, params_str)
        # Directory to store the daily HDF5 files
        self.psd_data_dir = os.path.join(self.run_output_dir, "psd_daily_h5")

    def to_serializable_dict(self):
        """Converts object attributes to a JSON-serializable dictionary."""
        output = {}
        for key, value in self.__dict__.items():
            if isinstance(value, (datetime, timedelta)):
                output[key] = str(value)
            elif isinstance(value, np.ndarray):
                output[key] = value.tolist()
            else:
                try:
                    json.dumps(value)
                    output[key] = value
                except TypeError:
                    output[key] = f"non-serializable type: {type(value).__name__}"
        return output

    def save_to_json(self):
        """Saves the configuration to a JSON file in the run output directory."""
        os.makedirs(self.run_output_dir, exist_ok=True)
        path = os.path.join(self.run_output_dir, 'config.json')
        serializable_dict = self.to_serializable_dict()

        with open(path, 'w', encoding='utf-8') as f:
            json.dump(serializable_dict, f, indent=4, ensure_ascii=False)
        print(f"Configuration saved to: {path}")

# ####################################################################################
# # --- 2. CALCULATION SCRIPT ---
# ####################################################################################

def calculate_and_save_psd_daily_hdf5(config: Config):
    """
    Calculates PSD and saves results into daily HDF5 files.
    """
    print(f"--- STAGE 1: Calculating and Saving PSDs ---")
    print(f"Analyzing data from {config.t_all_start} to {config.t_all_end}")
    print(f"Daily HDF5 files will be saved in: {config.psd_data_dir}")

    t_array = np.arange(
        np.datetime64(config.t_all_start),
        np.datetime64(config.t_all_end),
        config.chunk_duration
    )
    
    nperseg = int(config.nperseg_seconds * config.sampling_rate_hz)
    freq_vector = np.fft.rfftfreq(nperseg, d=1./config.sampling_rate_hz)

    # --- Main Processing Loop ---
    for t_start_dt64 in tqdm(t_array, desc="Processing Chunks"):
        t_start = t_start_dt64.astype(datetime)
        t_end = t_start + config.chunk_duration
        
        # --- Determine correct daily HDF5 file path ---
        date_str = t_start.strftime('%Y-%m-%d')
        daily_h5_path = os.path.join(config.psd_data_dir, f"{date_str}.h5")
        timestamp_str = t_start.strftime('%Y%m%dT%H%M%S')

        try:
            # --- Open the correct daily file and process the chunk ---
            # 'a' mode: read/write if exists, create otherwise.
            with h5py.File(daily_h5_path, 'a') as hf:
                # 1. Initialize the file if it's newly created
                if 'frequency' not in hf:
                    print(f"\nCreating new daily file: {daily_h5_path}")
                    hf.create_dataset('frequency', data=freq_vector)
                    for chan_km in config.dist_channels_km:
                        hf.create_group(f"ch_{chan_km:.1f}km")

                # 2. Check for existing data to allow resuming
                first_chan_group = f"ch_{config.dist_channels_km[0]:.1f}km"
                if timestamp_str in hf[first_chan_group]:
                    continue

                # 3. Load data (only if not skipped)
                das_data, headers, _ = load_das_h5_PE2024.load_das_custom(
                    t_start, t_end, input_dir=config.input_dir
                )
                
                das_data = das_data[:, config.channels_to_process_idx].astype(np.float64)
                
                current_fs = headers.get('fs', config.sampling_rate_hz)
                nperseg_current = int(config.nperseg_seconds * current_fs)

                if das_data is None or das_data.shape[0] < nperseg_current:
                    print(f"\nWarning: Skipping chunk at {t_start} due to insufficient data.")
                    continue

                # 4. Calculate and save PSD for each channel
                for j, chan_km in enumerate(config.dist_channels_km):
                    trace = das_data[:, j]
                    _, Pxx = signal.welch(trace, fs=current_fs, nperseg=nperseg_current, scaling='density')
                    
                    group_name = f"ch_{chan_km:.1f}km"
                    hf[group_name].create_dataset(timestamp_str, data=Pxx, compression="gzip")
        
        except Exception as e:
            error_time = t_start.strftime('%Y-%m-%d %H:%M:%S')
            print(f"\nError processing chunk at {error_time}: {e}")
        
        # --- Memory Management ---
        finally:
            if 'das_data' in locals():
                del das_data
            if 'headers' in locals():
                del headers
            gc.collect()

    print(f"\n--- STAGE 1 Complete. PSD calculation finished. ---")


if __name__ == '__main__':
    # 1. Initialize configuration
    config = Config()

    # 2. Create directories and save the configuration file for this run
    # The main run directory holds the config
    os.makedirs(config.run_output_dir, exist_ok=True) 
    # The sub-directory holds the daily H5 files
    os.makedirs(config.psd_data_dir, exist_ok=True)
    config.save_to_json()

    # 3. Start the calculation process
    calculate_and_save_psd_daily_hdf5(config)