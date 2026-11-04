# -*- coding: utf-8 -*-
"""
plot_ppsd.py

VERSION: v1.0
DESCRIPTION:
This script generates Probabilistic Power Spectral Density (PPSD) plots for
specified DAS channels using a pre-computed daily HDF5 file structure.
It reads all PSD data within a given time window for each channel and calculates
the statistical distribution of power levels at each frequency, visualizing it
as a 2D histogram.

This script is adapted to use the data loading logic from 'plot_results.py'.
"""
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import os
import sys
import h5py
from tqdm import tqdm
import json

# ==============================================================================
# SECTION 0: CONFIGURATION CLASS
# Manages all paths and parameters for the PPSD analysis.
# ==============================================================================
# --- [去硬编码] 包根推导 (相对 __file__); 本脚本不 import functions, 仅需输出目录 ---
# 包结构: check_minmag_pkg/{noise_psd, ...}; 本文件在 noise_psd/ 下。
_PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class Config:
    def __init__(self):
        # --- 1. Base Paths and Data Location ---
        # [去硬编码] 输出基目录改为包内 noise_psd/ 自身目录;
        #   原: os.path.join(HOME_dir, "Run_4_amp_spec_final/")
        self.output_dir_base = os.path.join(_PKG_ROOT, "noise_psd")
        
        # Parameters from the PSD calculation script (must match)
        self.nperseg_seconds = 200
        self.chunk_duration = timedelta(minutes=40)
        
        # Path to the directory containing daily H5 files
        params_str = f"chunk{self.chunk_duration}_nperseg{int(self.nperseg_seconds)}s"
        self.psd_data_dir = os.path.join(self.output_dir_base, params_str, 'psd_daily_h5')

        # --- 2. PPSD Analysis Parameters ---
        # Define channels to analyze and plot individually
        self.channels_to_plot = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8, 8.3, 8.6,
                                          9, 9.2, 9.4, 9.6, 9.8, 10, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 12.2, 12.4, 12.6, 12.8,
                                          13, 13.2, 13.4, 13.6, 13.8, 14,
                                          14.2, 14.4, 14.6, 14.8, 15, 15.2])

        # Define the time window for the PPSD calculation
        self.time_window = (
            datetime.strptime('2024-04-24 10:20:00', '%Y-%m-%d %H:%M:%S'),
            datetime.strptime('2024-06-22 17:10:00', '%Y-%m-%d %H:%M:%S')
        )
        
        # --- 3. Plotting and Output Parameters ---
        self.time_window_str = f"{self.time_window[0].strftime('%Y%m%d')}-{self.time_window[1].strftime('%Y%m%d')}"
        
        # Create a dedicated output directory for PPSD plots
        self.plot_output_dir = os.path.join(os.path.dirname(self.psd_data_dir), 'ppsd_plots', self.time_window_str)

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

    def save_config(self):
        """Saves the configuration to a JSON file in the plot output directory."""
        os.makedirs(self.plot_output_dir, exist_ok=True)
        path = os.path.join(self.plot_output_dir, 'ppsd_config.json')
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(self.to_serializable_dict(), f, indent=4, ensure_ascii=False)
        print(f"\nConfiguration saved to: {path}")

# ==============================================================================
# SECTION 1: DATA LOADING FUNCTION (from plot_results.py)
# ==============================================================================

def load_psd_data_from_daily_h5(psd_data_dir, channel_km, time_start, time_end):
    """
    Reads PSD data for a specific channel within a given time window from a
    directory of daily HDF5 files.
    """
    time_stamps, psd_list = [], []
    f_axis = None
    
    dates_to_check = [time_start.date() + timedelta(days=d) for d in range((time_end.date() - time_start.date()).days + 1)]

    print(f"-> Searching for data across {len(dates_to_check)} day(s) for channel {channel_km} km...")

    for date in tqdm(dates_to_check, desc=f"Loading Channel {channel_km}km"):
        daily_h5_path = os.path.join(psd_data_dir, f"{date.strftime('%Y-%m-%d')}.h5")
        if not os.path.exists(daily_h5_path):
            continue
        try:
            with h5py.File(daily_h5_path, 'r') as hf:
                if f_axis is None and 'frequency' in hf:
                    f_axis = hf['frequency'][:]
                
                group_name = f"ch_{channel_km:.1f}km"
                if group_name not in hf:
                    continue
                
                for ts_str in hf[group_name].keys():
                    try:
                        ts = datetime.strptime(ts_str, '%Y%m%dT%H%M%S')
                        if time_start <= ts < time_end:
                            time_stamps.append(ts)
                            psd_list.append(hf[group_name][ts_str][:])
                    except (ValueError, KeyError):
                        continue
        except Exception as e:
            print(f"Warning: Could not read file {daily_h5_path}. Error: {e}")

    if not psd_list:
        return [], None, None
    
    # Data is not guaranteed to be time-sorted, so sort it now.
    sorted_indices = np.argsort(time_stamps)
    sorted_times = [time_stamps[i] for i in sorted_indices]
    sorted_psds = [psd_list[i] for i in sorted_indices]
    
    # Stack into a matrix with shape (n_times, n_frequencies)
    psd_matrix = np.array(sorted_psds)
    
    print(f"-> Loaded {psd_matrix.shape[0]} PSD segments for channel {channel_km} km.")
    return sorted_times, f_axis, psd_matrix


# ==============================================================================
# SECTION 2: PPSD CALCULATION AND PLOTTING FUNCTION
# ==============================================================================

def plot_single_channel_ppsd(channel_km, f_axis, psd_matrix, cfg: Config):
    """
    Calculates and plots the PPSD for a single channel's data.
    """
    if psd_matrix is None or psd_matrix.size == 0 or f_axis is None:
        print(f"  因没有数据，跳过通道 {channel_km} 的PPSD绘图。")
        return

    print(f"-> 正在为通道 {channel_km} km 计算PPSD...")
    
    # --- Font Sizes ---
    # 在这里集中调整所有绘图相关的字体大小
    title_fontsize = 24
    label_fontsize = 24
    tick_fontsize = 24
    cbar_label_fontsize = 24
    
    # --- 1. Prepare Data ---
    # Convert PSD to dB relative to 1 (strain/s)^2/Hz
    # Use a small epsilon to avoid log(0) errors
    psd_db = 10 * np.log10(psd_matrix + 1e-25) 
    
    # Define bins for the 2D histogram
    db_bins = np.linspace(-220, -80, 141) # Y-axis (dB)
    
    # Create logarithmically spaced bins for the frequency axis
    valid_freq_idx = f_axis > 0
    if not np.any(valid_freq_idx):
        print(f"  跳过通道 {channel_km} 的PPSD：未找到正频率。")
        return
    min_freq = f_axis[valid_freq_idx].min()
    max_freq = f_axis.max()
    freq_bins_log = np.logspace(np.log10(min_freq), np.log10(max_freq), 100)

    # --- 2. Calculate 2D Histogram ---
    hist2d = np.zeros((len(db_bins)-1, len(freq_bins_log)-1))

    # Iterate through each frequency bin to build the histogram
    for i in tqdm(range(len(freq_bins_log)-1), desc="  正在构建直方图"):
        f_start, f_end = freq_bins_log[i], freq_bins_log[i+1]
        
        # Find indices of the original frequency axis that fall in this bin
        freq_indices = np.where((f_axis >= f_start) & (f_axis < f_end))[0]
        if freq_indices.size == 0:
            continue
            
        # Get all PSD values (in dB) within this frequency slice
        psd_in_f_bin = psd_db[:, freq_indices].flatten()
        
        # Calculate the histogram for this slice and add it to our 2D array
        hist, _ = np.histogram(psd_in_f_bin, bins=db_bins)
        hist2d[:, i] = hist
        
    # --- 3. Normalize to get Probability Density ---
    total_counts_per_freq = hist2d.sum(axis=0)
    # Avoid division by zero for frequency bins with no data
    total_counts_per_freq[total_counts_per_freq == 0] = 1 
    
    ppsd = (hist2d / total_counts_per_freq) * 100 # Normalize to 100%

    # --- 4. Plotting ---
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Use pcolormesh to display the 2D histogram
    X, Y = np.meshgrid(freq_bins_log, db_bins)
    im = ax.pcolormesh(X, Y, ppsd, cmap='YlGnBu', shading='auto', 
                       vmax=np.percentile(ppsd[ppsd>0], 95) if np.any(ppsd > 0) else 1) # Cap color at 95th percentile

    # --- 5. Formatting and Labels ---
    ax.set_xscale('log')
    ax.set_xlabel('Frequency [Hz]', fontsize=label_fontsize)
    ax.set_ylabel(r'PSD [dB, rel. to 1 $(strain/s)^2/Hz$]', fontsize=label_fontsize)
    
    
    title_str = f"PPSD for Channel {channel_km:.1f} km\n" \
                # f"Time Window: {cfg.time_window[0].strftime('%Y-%m-%d')} to {cfg.time_window[1].strftime('%Y-%m-%d')}"
    ax.set_title(title_str, fontsize=title_fontsize)
    
    ax.set_xlim(min_freq, max_freq)
    ax.set_ylim(db_bins.min(), db_bins.max())
    ax.grid(True, which='both', linestyle=':', linewidth=0.5)
    
    # Set tick font size
    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize)

    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label('PSD Probability Density [%]', size=cbar_label_fontsize)
    cbar.ax.tick_params(labelsize=cbar_label_fontsize)
    
    fig.tight_layout()
    
    # --- 6. Save Figure ---
    channel_str = f"ch{channel_km:.1f}km".replace('.', 'p')
    output_filename = os.path.join(cfg.plot_output_dir, f'ppsd_{channel_str}_{cfg.time_window_str}.png')
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ 通道 {channel_km} km 的 PPSD 图像已保存至：{output_filename}")
# ==============================================================================
# SECTION 3: SUMMARY PPSD PERCENTILE PLOT (MODIFIED FOR FONT SIZE)
# ==============================================================================

def plot_summary_ppsd_percentile(all_psd_data, f_axis, percentile, cfg: Config):
    """
    Plots a 2D heatmap of a specific PSD percentile across all channels.

    The plot shows channel distance on the y-axis, frequency on the x-axis,
    and the color represents the specified percentile of PSD in dB.
    """
    if not all_psd_data or f_axis is None:
        print("Warning: Cannot generate summary plot. No data loaded for any channel.")
        return

    print("\n" + "="*50)
    print(f"Generating summary plot for the {percentile}th percentile...")

    # --- 1. Prepare data for plotting ---
    sorted_channels = sorted(all_psd_data.keys())
    percentile_values_list = []

    for chan_km in tqdm(sorted_channels, desc="  Calculating percentiles"):
        psd_matrix = all_psd_data[chan_km]
        
        if psd_matrix is None or psd_matrix.size == 0:
            percentile_values_list.append(np.full(f_axis.shape, np.nan))
            continue
        
        # Convert to dB (same unit as individual PPSD plots)
        psd_db = 10 * np.log10(psd_matrix + 1e-25)
        
        percentiles_for_channel = np.percentile(psd_db, percentile, axis=0)
        percentile_values_list.append(percentiles_for_channel)
    
    percentile_matrix = np.vstack(percentile_values_list)

    # --- 2. Create the plot ---
    fig, ax = plt.subplots(figsize=(10, 6))

    valid_freq_idx = f_axis > 0
    f_axis_plot = f_axis[valid_freq_idx]
    percentile_matrix_plot = percentile_matrix[:, valid_freq_idx]

    if f_axis_plot.size == 0:
        print("Warning: No valid positive frequencies to plot.")
        plt.close(fig)
        return

    X, Y = np.meshgrid(f_axis_plot, sorted_channels)
    
    vmin = np.nanpercentile(percentile_matrix_plot, 5)
    vmax = np.nanpercentile(percentile_matrix_plot, 95)

    im = ax.pcolormesh(X, Y, percentile_matrix_plot, cmap='viridis', shading='auto', vmin=vmin, vmax=vmax)

    # --- 3. Formatting and Labels (ADJUSTED FOR LARGER FONTS) ---
    # Define font sizes for better readability in documents
    TICK_LABEL_SIZE = 24
    AXIS_LABEL_SIZE = 24
    TITLE_SIZE = 24
    CBAR_LABEL_SIZE = 24

    ax.set_xscale('log')
    ax.set_xlabel('Frequency [Hz]', fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel('Channel Distance [km]', fontsize=AXIS_LABEL_SIZE)
    ylabel = np.arange(1, 15.5, 2)  # Y-axis labels for channel distances
    ax.set_yticks(ylabel)
    ax.set_yticklabels([f"{y:.0f}" for y in ylabel], fontsize=TICK_LABEL_SIZE)
    
    title_str = f"PPSD at {percentile}th Percentile\n" \
                # f"Time Window: {cfg.time_window[0].strftime('%Y-%m-%d')} to {cfg.time_window[1].strftime('%Y-%m-%d')}"
    ax.set_title(title_str, fontsize=TITLE_SIZE, fontweight='bold')
    
    ax.set_xlim(f_axis_plot.min(), f_axis_plot.max())
    ax.grid(True, which='both', linestyle=':', linewidth=0.5)

    # Adjust tick label sizes
    ax.tick_params(axis='both', which='major', labelsize=TICK_LABEL_SIZE)
    

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    # The colorbar label clearly states the unit is dB, consistent with the other plot
    cbar.set_label(f'[dB, rel. to 1 $(strain/s)^2/Hz$]', size=CBAR_LABEL_SIZE)
    cbar.ax.tick_params(labelsize=TICK_LABEL_SIZE) # Set colorbar tick label size
    
    fig.tight_layout(pad=1.5) # Add padding to ensure labels don't get cut off

    # --- 4. Save Figure ---
    output_filename = os.path.join(cfg.plot_output_dir, f'ppsd_summary_{percentile}percentile_{cfg.time_window_str}.png')
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ Summary percentile plot saved to: {output_filename}")


# ==============================================================================
# SECTION 4: MAIN EXECUTION BLOCK (MODIFIED)
# ==============================================================================
if __name__ == '__main__':
    #* 1. Initialize configuration
    cfg = Config()
    
    #* 2. Create the main output directory
    os.makedirs(cfg.plot_output_dir, exist_ok=True)
    
    print("--- Starting PPSD Plotting Script ---")
    print(f"Reading data from: {cfg.psd_data_dir}")
    print(f"Saving plots to:   {cfg.plot_output_dir}")
    
    if not os.path.isdir(cfg.psd_data_dir):
        print(f"Error: PSD data directory not found at '{cfg.psd_data_dir}'")
        sys.exit("Please check the path in the Config class.")

    #* 3. Phase 1: Load data for all specified channels
    all_psd_data = {}
    f_axis = None
    
    
    print("\n--- Phase 1: Loading data for all specified channels ---")
    for chan_km in cfg.channels_to_plot:
        # Load all PSD data for the current channel and time window
        times, current_f_axis, psd_matrix = load_psd_data_from_daily_h5(
            cfg.psd_data_dir, 
            chan_km, 
            cfg.time_window[0], 
            cfg.time_window[1]
        )
        
        # Store the loaded data
        all_psd_data[chan_km] = psd_matrix
        
        # Store the frequency axis (should be the same for all)
        if f_axis is None and current_f_axis is not None:
            f_axis = current_f_axis
            
            
            
    #* 4. Phase 2: Generate individual PPSD plots for each channel
    print("\n--- Phase 2: Generating individual PPSD plots ---")
    for chan_km in cfg.channels_to_plot:
        psd_matrix = all_psd_data.get(chan_km) # Use .get for safety
        
        # Calculate and plot the PPSD (original function)
        plot_single_channel_ppsd(chan_km, f_axis, psd_matrix, cfg)
    
    
    
    #* 5. Phase 3: Generate the summary percentile plot
    print("\n--- Phase 3: Generating summary percentile plot ---")
    # You can change the percentile value here, e.g., 50, 90, 95
    percentile_to_plot = 50
    plot_summary_ppsd_percentile(all_psd_data, f_axis, percentile_to_plot, cfg)
    
    # Example: Plot another percentile if needed
    # plot_summary_ppsd_percentile(all_psd_data, f_axis, 95, cfg)

    #* 6. Save the configuration for reproducibility
    cfg.save_config()
    print("\n--- All channels processed. Script finished. ---")