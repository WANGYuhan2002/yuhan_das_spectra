# -*- coding: utf-8 -*-
"""
analyze_noise_trends.py

VERSION: v1.1
DESCRIPTION:
This script analyzes the temporal evolution of seismic noise for specified DAS channels.
It processes PSD data in discrete time steps, calculates a noise metric within a
defined frequency band (e.g., the 95th percentile of power), and converts it to
an amplitude proxy. The results for all channels over the entire analysis period
are stored in a single HDF5 file and visualized in two ways:
1. A 2D heatmap showing all channels.
2. A multi-panel plot with individual time series for each channel.

This script relies on the data loading function from the original PPSD script.
"""
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from datetime import datetime, timedelta
import os
import sys
import h5py
from tqdm import tqdm
import json
from matplotlib.colors import LogNorm

# --- [去硬编码] 包根推导 + functions 路径 (相对 __file__) ---
# 包结构: check_minmag_pkg/{functions, noise_psd, ...}; 本文件在 noise_psd/ 下。
_PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(_PKG_ROOT, "functions"))
import load_das_h5_PE2024 as load_das_h5_PE2024
from plot import Ticks_xaxis_day_or_hour, Ticks_xaxis_minute_or_second

# ==============================================================================
# SECTION 0: CONFIGURATION CLASS
# Manages all paths and parameters for the noise trend analysis.
# ==============================================================================
class Config:
    def __init__(self):
        # --- 1. Base Paths and Data Location ---
        # [去硬编码] 输出基目录改为包内 noise_psd/ 自身目录;
        #   原: os.path.join(HOME_dir, "Run_4_amp_spec_final/")
        # 重算产出的 noise_trends_*.h5 即落在此目录下 (与 products/ 内冻结成品同结构)。
        self.output_dir_base = os.path.join(_PKG_ROOT, "noise_psd")

        # Parameters from the PSD calculation script (must match)
        self.nperseg_seconds = 200
        self.chunk_duration = timedelta(minutes=40)

        # Path to the directory containing daily H5 files
        params_str = f"chunk{self.chunk_duration}_nperseg{int(self.nperseg_seconds)}s"
        self.psd_data_dir = os.path.join(self.output_dir_base, params_str, 'psd_daily_h5')

        # --- 2. Noise Trend Analysis Parameters ---
        # Define the overall time window for the entire analysis
        self.analysis_time_window = (
            # datetime(2024, 5, 2, 12, 0, 0),
            # datetime(2024, 5, 5, 12, 0, 0)
            
            datetime.strptime('2024-04-24 10:20:00', '%Y-%m-%d %H:%M:%S'),
            datetime.strptime('2024-06-22 17:10:00', '%Y-%m-%d %H:%M:%S')
            
            # datetime(2024, 5, 4, 0, 0, 0),
            # datetime(2024, 5, 5, 0, 0, 0)
        )

        # Define the time resolution for the analysis (e.g., 1 day, 6 hours)
        self.time_step_duration =  self.chunk_duration# timedelta(hours=0, minutes=1)  # 1 minute

        # x-axis tick formatting
        self.Day_tick, self.Hour_tick = 5, 0
        self.Minute_tick, self.Second_tick = 0, 0

        # Define the frequency band of interest for noise calculation [Hz]
        self.frequency_band =(8, 20)
        # self.frequency_band =(0.02, 1)  # For lower frequency noise analysis

        # Define the percentile to represent the upper limit of noise (e.g., 95 for P95)
        self.percentile_to_calculate = 5

        # Define channels to analyze
        # self.channels_to_analyze = [2.0, 4.0, 6.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
        self.channels_to_analyze = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8, 8.3, 8.6,
                                          9, 9.2, 9.4, 9.6, 9.8, 10, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 12.2, 12.4, 12.6, 12.8,
                                          13, 13.2, 13.4, 13.6, 13.8, 14,
                                          14.2, 14.4, 14.6, 14.8, 15, 15.2])
        
        self.channels_to_average_noise_level = [8.3, 8.6,
                                          9, 9.2, 9.4, 9.6, 9.8, 10, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 12.2, 12.4, 12.6, 12.8,
                                          13, 13.2, 13.4, 13.6, 13.8, 14]
        # --- 3. Output Parameters ---
        self.time_window_str = f"{self.analysis_time_window[0].strftime('%Y%m%d_%H%M%S')}-{self.analysis_time_window[1].strftime('%Y%m%d_%H%M%S')}"
        
        # Create a dedicated output directory for this analysis
        self.analysis_output_dir = os.path.join(self.output_dir_base, params_str, 'noise_trend_plot', self.time_window_str + '_step' + str(self.time_step_duration))

        # Define the name for the output HDF5 file and the final plots
        plot_basename = f'f{self.frequency_band[0]}-{self.frequency_band[1]}Hz_p{self.percentile_to_calculate}'
        self.output_h5_filename = os.path.join(self.analysis_output_dir, f'noise_trends_{plot_basename}.h5')
        self.output_heatmap_plot_filename = os.path.join(self.analysis_output_dir, f'noise_heatmap_{plot_basename}.png')
        self.output_subplots_plot_filename = os.path.join(self.analysis_output_dir, f'noise_subplots_{plot_basename}.png')


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
        """Saves the configuration to a JSON file in the analysis output directory."""
        os.makedirs(self.analysis_output_dir, exist_ok=True)
        path = os.path.join(self.analysis_output_dir, 'noise_analysis_config.json')
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(self.to_serializable_dict(), f, indent=4, ensure_ascii=False)
        print(f"\nConfiguration saved to: {path}")

# ==============================================================================
# SECTION 1: DATA LOADING FUNCTION (Identical to plot_ppsd.py)
# ==============================================================================

def load_psd_data_from_daily_h5(psd_data_dir, channel_km, time_start, time_end):
    """
    Reads PSD data for a specific channel within a given time window from a
    directory of daily HDF5 files.
    """
    time_stamps, psd_list = [], []
    f_axis = None
    dates_to_check = [time_start.date() + timedelta(days=d) for d in range((time_end.date() - time_start.date()).days + 1)]

    for date in dates_to_check:
        daily_h5_path = os.path.join(psd_data_dir, f"{date.strftime('%Y-%m-%d')}.h5")
        if not os.path.exists(daily_h5_path):
            continue
        try:
            with h5py.File(daily_h5_path, 'r') as hf:
                if f_axis is None and 'frequency' in hf:
                    f_axis = hf['frequency'][:]
                group_name = f"ch_{channel_km:.1f}km"
                if group_name not in hf: continue
                for ts_str in hf[group_name].keys():
                    try:
                        ts = datetime.strptime(ts_str, '%Y%m%dT%H%M%S')
                        if time_start <= ts < time_end:
                            time_stamps.append(ts)
                            psd_list.append(hf[group_name][ts_str][:])
                    except (ValueError, KeyError): continue
        except Exception as e:
            print(f"Warning: Could not read file {daily_h5_path}. Error: {e}")

    if not psd_list: return [], None, None
    
    psd_matrix = np.array(psd_list)
    return time_stamps, f_axis, psd_matrix


# ==============================================================================
# SECTION 2: NOISE CALCULATION AND STORAGE
# ==============================================================================

def calculate_noise_level(psd_matrix, f_axis, freq_band, percentile):
    """
    Calculates the noise amplitude proxy from a PSD matrix for a given frequency band.

    Returns:
        float: The noise amplitude proxy. Returns np.nan if data is insufficient.
    """
    if psd_matrix is None or psd_matrix.size == 0 or f_axis is None:
        return np.nan

    # 1. Find frequency indices within the specified band
    freq_indices = np.where((f_axis >= freq_band[0]) & (f_axis <= freq_band[1]))[0]
    if freq_indices.size == 0:
        return np.nan

    # 2. Extract all PSD values within this frequency band
    psd_in_band = psd_matrix[:, freq_indices]

    # 3. Calculate the specified percentile of the power
    # This gives us the Pxx value that represents the high-noise level
    power_at_percentile = np.percentile(psd_in_band, percentile)

    # 4. Convert power to amplitude proxy by taking the square root.
    # The unit is (strain/s)/sqrt(Hz). It's an amplitude spectral density.
    noise_amplitude = np.sqrt(power_at_percentile)
    
    return noise_amplitude

def save_results_to_h5(filepath, results_matrix, time_steps, channels_km, cfg):
    """
    Saves the noise trend results to a well-structured HDF5 file.
    """
    print(f"\n-> Saving results to HDF5 file: {filepath}")
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with h5py.File(filepath, 'w') as f:
        # --- Store main data ---
        f.create_dataset('noise_amplitude', data=results_matrix, 
                         chunks=(10, len(channels_km)), compression="gzip")
        
        # --- Store coordinate data ---
        # Store time steps (start of each window) as UTC strings
        time_data = np.array([dt.strftime('%Y-%m-%dT%H:%M:%S') for dt in time_steps], dtype='S')
        f.create_dataset('timestamps_utc_str', data=time_data)
        f.create_dataset('channels_km', data=np.array(channels_km))

        # --- Store metadata as attributes for reproducibility ---
        f.attrs['description'] = "Noise amplitude trend analysis results."
        f.attrs['amplitude_unit'] = "(strain/s)/sqrt(Hz)"
        f.attrs['amplitude_description'] = f"sqrt(P{cfg.percentile_to_calculate} of PSD)"
        f.attrs['frequency_band_hz'] = str(cfg.frequency_band)
        f.attrs['time_step_duration'] = str(cfg.time_step_duration)
        f.attrs['creation_date_utc'] = datetime.utcnow().isoformat()
    print("✅ HDF5 file saved successfully.")

# ==============================================================================
# SECTION 3A: HEATMAP PLOTTING FUNCTION (CORRECTED VERSION)
# ==============================================================================
from matplotlib.colors import LogNorm
def plot_noise_heatmap(h5_filepath, output_plot_path, cfg: Config):
    """
    Reads the HDF5 results file and creates a 2D heatmap of noise trends.
    (This version is corrected to handle pcolormesh dimension requirements.)
    """
    print(f"\n-> 正在从以下文件生成热力图: {h5_filepath}")
    try:
        with h5py.File(h5_filepath, 'r') as f:
            noise_matrix = f['noise_amplitude'][:]
            channels = f['channels_km'][:]
            time_strs = [t.decode('utf-8') for t in f['timestamps_utc_str'][:]]
            time_steps = [datetime.strptime(ts, '%Y-%m-%dT%H:%M:%S') for ts in time_strs]

    except Exception as e:
        print(f"错误: 无法读取HDF5文件 '{h5_filepath}'。 {e}")
        return

    # --- Font Sizes ---
    # 在这里集中调整所有绘图相关的字体大小
    title_fontsize = 24
    label_fontsize = 24
    cbar_label_fontsize = 24
    day_hour_tick_fontsize = 24
    min_sec_tick_fontsize = 24
    
    # fig, ax = plt.subplots(figsize=(15, 8))
    fig = plt.figure(figsize=(18, 10))
    # fig = plt.figure(figsize=(90, 32))
    
    gs = plt.GridSpec(2,2,figure=fig,height_ratios=[8,3],width_ratios=[11,1]) # adapt the gs according to the number of templates
    axes = [fig.add_subplot(gs[0,:]), fig.add_subplot(gs[1,0])] 
    # axes = [fig.add_subplot(211), fig.add_subplot(212)]  # Two subplots: heatmap and average noise level
    
    # <<< pcolormesh 坐标修正部分 (保持不变) >>>
    # pcolormesh requires coordinate arrays to be one larger than the data array.
    # We need to define the corners of the pixels.

    # 1. For the X-axis (time), create N+1 coordinates for N data points.
    x_coords = time_steps + [time_steps[-1] + cfg.time_step_duration]

    # 2. For the Y-axis (channels), create M+1 coordinates for M channels.
    channels_arr = np.array(channels)
    if len(channels_arr) > 1:
        # Calculate spacing from the difference between the last two channels
        d_chan = channels_arr[-1] - channels_arr[-2]
    else:
        # If only one channel, assume a default spacing (e.g., 1.0)
        d_chan = 1.0
    # Append the top edge of the last channel's pixel
    y_coords = np.append(channels_arr, channels_arr[-1] + d_chan)

    # The data matrix needs to be transposed to (channels, time) for pcolormesh
    C = noise_matrix.T
    # <<< 修正部分结束 >>>
    
    # Find a robust min/max for the color scale, ignoring NaNs and extreme outliers
    valid_data = C[~np.isnan(C)]
    if valid_data.size == 0:
        print("警告: 没有有效数据可以绘制。")
        return
        
    vmin = np.percentile(valid_data, 5)
    vmax = np.percentile(valid_data, 95)
    
    ax = axes[0]
    # Use the corrected 1D coordinate arrays directly. matplotlib handles the grid.
    im = ax.pcolormesh(x_coords, y_coords, C, shading='auto', norm=LogNorm(vmin=vmin, vmax=vmax), cmap='magma')

    # --- Formatting and Labels ---
    ax.set_ylabel('Channel Distance [km]', fontsize=label_fontsize)
    
    title_str = (f"Data Amplitude Trend ({cfg.frequency_band[0]}-{cfg.frequency_band[1]} Hz, P{cfg.percentile_to_calculate})\n"
                 f"Time Resolution: {str(cfg.time_step_duration)}")
    ax.set_title(title_str, fontsize=title_fontsize)

    # Format the x-axis to display dates nicely
    if cfg.Day_tick > 0 or cfg.Hour_tick > 0:
        Ticks_xaxis_day_or_hour(ax, cfg.Day_tick, cfg.Hour_tick, rotation=30, font=day_hour_tick_fontsize)
    if cfg.Minute_tick > 0 or cfg.Second_tick > 0:
        Ticks_xaxis_minute_or_second(ax, cfg.Minute_tick, cfg.Second_tick, rotation=0, font=min_sec_tick_fontsize)

    cbar = fig.colorbar(im, ax=ax, pad=0.01)
    cbar.set_label(f'Data Amplitude Proxy [$((m/m)/s)/\\sqrt{{Hz}}$]', size=cbar_label_fontsize)
    cbar.ax.tick_params(labelsize=cbar_label_fontsize)
    
    ax.set_ylim(y_coords.min(), y_coords.max()) # Use the new y_coords for limits
    ax.grid(True, which='major', axis='x', linestyle=':', linewidth=0.5, color='white', alpha=0.5)
    ax.set_xticks([])
    ax.tick_params(axis='y', labelsize=day_hour_tick_fontsize)
    
    ax = axes[1]
    all_channels_list = list(channels)
    indices_to_avg = [all_channels_list.index(c) for c in cfg.channels_to_average_noise_level if c in all_channels_list]

    if not indices_to_avg:
        print("警告: 未找到用于平均的通道，将跳过平均值图的绘制。")
        ax.text(0.5, 0.5, 'No channels found for averaging.', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    else:
        # 2. 选择这些通道的数据列
        data_to_avg = noise_matrix[:, indices_to_avg]
        
        # 3. 计算每个时间步长上这些通道的平均值，忽略NaN
        average_noise = np.nanmean(data_to_avg, axis=1)
        
    # 4. 绘制平均噪声水平
    ax.plot(time_steps, average_noise, color='k', marker='.', markersize=2, linestyle='-', linewidth=1.5)

    ax.set_title(f'Data amplitude averaged over {cfg.channels_to_average_noise_level[0]} - {cfg.channels_to_average_noise_level[-1]} km [$((m/m)/s)/\\sqrt{{Hz}}$]', fontsize=label_fontsize)
    
    ax.grid(True, which='both', linestyle=':', linewidth=0.8)
    ax.set_yscale('log')
    ax.tick_params(axis='y', labelsize=day_hour_tick_fontsize)

    # --- 通用坐标轴格式化 ---
    # ax.set_xlabel('Date', fontsize=label_fontsize)

    # 格式化共享的x轴刻度
    if cfg.Day_tick > 0 or cfg.Hour_tick > 0:
        Ticks_xaxis_day_or_hour(ax, cfg.Day_tick, cfg.Hour_tick, rotation=15, font=day_hour_tick_fontsize)
    if cfg.Minute_tick > 0 or cfg.Second_tick > 0:
        Ticks_xaxis_minute_or_second(ax, cfg.Minute_tick, cfg.Second_tick, rotation=0, font=day_hour_tick_fontsize)
    
    # 为两个图设置x轴限制
    ax.set_xlim(time_steps[0], time_steps[-1])
    ax.tick_params(axis='x', labelsize=day_hour_tick_fontsize)
    

        
    
    # fig.tight_layout()

    # --- Save Figure ---
    fig.savefig(output_plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ 热力图已保存至: {output_plot_path}")


# ==============================================================================
# SECTION 3B: SUBPLOTS PLOTTING FUNCTION
# ==============================================================================

def plot_noise_subplots(h5_filepath, output_plot_path, cfg):
    """
    Reads the HDF5 results file and creates a multi-panel plot where each 
    channel's noise trend is in its own subplot.
    """
    print(f"\n-> Generating subplots from: {h5_filepath}")
    try:
        with h5py.File(h5_filepath, 'r') as f:
            noise_matrix = f['noise_amplitude'][:]
            channels = f['channels_km'][:]
            time_strs = [t.decode('utf-8') for t in f['timestamps_utc_str'][:]]
            
            # --- Convert time strings back to datetime objects for plotting ---
            time_steps = [datetime.strptime(ts, '%Y-%m-%dT%H:%M:%S') for ts in time_strs]

    except Exception as e:
        print(f"Error: Could not read HDF5 file '{h5_filepath}'. {e}")
        return

    num_channels = len(channels)
    if num_channels == 0:
        print("Warning: No channels to plot.")
        return

    # --- Dynamically adjust figure height based on the number of channels ---
    # We allocate ~2.5 inches of height per subplot for clarity
    fig_height = num_channels * 2.0
    fig, axes = plt.subplots(
        nrows=num_channels, 
        ncols=1, 
        figsize=(18, fig_height), 
        sharex=True  # Share the x-axis (time) across all subplots
    )
    fsz = 20

    # If there's only one channel, plt.subplots returns a single Axes object, not an array
    if num_channels == 1:
        axes = [axes]

    # --- Main plotting loop for each channel ---
    for i, ax in enumerate(axes):
        channel_km = channels[i]
        noise_data = noise_matrix[:, i] # Get the data column for this channel

        ax.plot(time_steps, noise_data, marker='o', markersize=3, linestyle='-', label=f'Channel {channel_km:.1f} km')
        
        # Set a specific y-label for each subplot indicating the channel distance
        ax.set_ylabel(f"{channel_km:.1f} km", fontsize=fsz)
        ax.grid(True, which='major', linestyle=':', linewidth=0.8)
        ax.tick_params(axis='y', labelsize=fsz)
        # Only show data points that are not NaN
        ax.set_xlim(time_steps[0], time_steps[-1])
        
        if cfg.Day_tick > 0 or cfg.Hour_tick > 0:
            Ticks_xaxis_day_or_hour(ax, cfg.Day_tick, cfg.Hour_tick, rotation=30, font=fsz)
        if cfg.Minute_tick > 0 or cfg.Second_tick > 0:
            Ticks_xaxis_minute_or_second(ax, cfg.Minute_tick, cfg.Second_tick, rotation=0,font=fsz)
            
        ax.set_yscale('log')  # Use logarithmic scale for better visibility of noise levels
        ax.set_ylim(bottom=1e-10)  # Set a lower limit for the y-axis

    # --- Formatting the entire figure ---
    
    # Add a shared, main y-axis label for the whole figure
    fig.text(0.06, 0.5, f'Noise Amplitude Proxy [$(strain/s)/\\sqrt{{Hz}}$]', 
             va='center', rotation='vertical', fontsize=fsz)

    # Add a main title to the figure
    title_str = (f"Noise Amplitude Trend per Channel ({cfg.frequency_band[0]}-{cfg.frequency_band[1]} Hz, P{cfg.percentile_to_calculate})\n"
                 f"Time Resolution: {str(cfg.time_step_duration)}")
    fig.suptitle(title_str, fontsize=fsz, y=0.99)


    # Adjust layout to prevent labels/titles from overlapping
    fig.tight_layout(rect=[0.08, 0.05, 1, 0.97]) # Adjust rect to make space for fig.text and suptitle

    # --- Save Figure ---
    fig.savefig(output_plot_path, dpi=300)
    plt.close(fig)
    print(f"✅ Subplots plot saved to: {output_plot_path}")


# ==============================================================================
# SECTION 4: MAIN EXECUTION BLOCK
# ==============================================================================
if __name__ == '__main__':
    # 1. Initialize configuration and create output directory
    cfg = Config()
    os.makedirs(cfg.analysis_output_dir, exist_ok=True)
    
    print("--- Starting Noise Trend Analysis Script ---")
    print(f"Reading data from: {cfg.psd_data_dir}")
    print(f"Saving results to: {cfg.analysis_output_dir}")
    
    if not os.path.isdir(cfg.psd_data_dir):
        print(f"Error: PSD data directory not found at '{cfg.psd_data_dir}'")
        sys.exit("Please check the path in the Config class.")

    # 2. Generate time steps for the analysis
    time_steps = []
    current_time = cfg.analysis_time_window[0]
    while current_time < cfg.analysis_time_window[1]:
        time_steps.append(current_time)
        current_time += cfg.time_step_duration
    
    num_steps = len(time_steps)
    num_channels = len(cfg.channels_to_analyze)
    print(f"\nAnalysis will be run for {num_channels} channels over {num_steps} time steps.")

    # 3. Initialize results matrix with NaNs
    results_matrix = np.full((num_steps, num_channels), np.nan)

    # 4. Main processing loop
    # Outer loop over channels, inner loop over time (often more efficient for data loading)
    for j, chan_km in enumerate(tqdm(cfg.channels_to_analyze, desc="Processing Channels")):
        for i, t_start in enumerate(tqdm(time_steps, desc=f"  Channel {chan_km}km", leave=False)):
            t_end = t_start + cfg.time_step_duration
            
            # Load data for this specific channel and time window
            _, f_axis, psd_matrix = load_psd_data_from_daily_h5(
                cfg.psd_data_dir, 
                chan_km, 
                t_start, 
                t_end
            )
            
            # Calculate the noise metric
            noise_val = calculate_noise_level(
                psd_matrix, 
                f_axis, 
                cfg.frequency_band, 
                cfg.percentile_to_calculate
            )
            
            results_matrix[i, j] = noise_val

    # 5. Save the computed results to an HDF5 file
    save_results_to_h5(
        cfg.output_h5_filename, 
        results_matrix, 
        time_steps, 
        cfg.channels_to_analyze, 
        cfg
    )

    # 6. Plot the results from the HDF5 file as a heatmap
    plot_noise_heatmap(cfg.output_h5_filename, cfg.output_heatmap_plot_filename, cfg)
    
    # 7. Plot the results in subplots 
    plot_noise_subplots(cfg.output_h5_filename, cfg.output_subplots_plot_filename, cfg)

    # 8. Save the configuration for reproducibility
    cfg.save_config()
    print("\n--- Analysis and plotting complete. Script finished. ---")