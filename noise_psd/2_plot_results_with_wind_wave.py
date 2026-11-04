# -*- coding: utf-8 -*-
"""
plot_results.py

VERSION: v3.0
DESCRIPTION:
This script visualizes pre-computed PSD data stored in a daily HDF5 file structure.
It loads data from an arbitrary time window, which can span multiple days/files,
and generates spectrograms and binned time-series plots.

Key Features:
- A rewritten `load_psd_data_from_daily_h5` function to read from the new
  daily HDF5 file structure.
- A `Config` class to manage all plotting parameters.
- Saves plots and the run configuration for full reproducibility.
"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
import numpy as np
from datetime import datetime, timedelta
import os
import sys
import h5py # Added for HDF5 reading
from tqdm import tqdm
import json

# --- [去硬编码] 包根推导 + functions / meteo 路径 (相对 __file__) ---
# 包结构: check_minmag_pkg/{functions, meteo, noise_psd, ...}; 本文件在 noise_psd/ 下。
_PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(_PKG_ROOT, "functions"))
import load_das_h5_PE2024 as load_das_h5_PE2024
from plot import Ticks_xaxis_day_or_hour, Ticks_xaxis_minute_or_second

# [去硬编码] 气象加载模块改指本包 meteo/ (原指 check_minmag_validity/.../wind_speed/)。
sys.path.append(os.path.join(_PKG_ROOT, "meteo"))
from load_kefa_meteo import load_kefa_meteo_data

# ####################################################################################
# # --- 1. CORE DATA LOADING FUNCTION (Rewritten for Daily HDF5) ---
# ####################################################################################

def load_psd_data_from_daily_h5(psd_data_dir, channel_km, time_start, time_end):
    """
    Reads PSD data for a specific channel within a given time window from a
    directory of daily HDF5 files.

    Args:
        psd_data_dir (str): The directory containing the daily HDF5 files.
        channel_km (float): The channel distance to load.
        time_start (datetime): The start of the time window.
        time_end (datetime): The end of the time window.

    Returns:
        tuple: A tuple containing (time_stamps, f_axis, psd_matrix).
               Returns ([], None, None) if no data is found.
    """
    time_stamps, psd_list = [], []
    f_axis = None
    
    # 1. Generate the list of dates to check based on the time window
    dates_to_check = []
    current_date = time_start.date()
    while current_date <= time_end.date():
        dates_to_check.append(current_date)
        current_date += timedelta(days=1)

    print(f"  -> Searching for data across {len(dates_to_check)} day(s) for channel {channel_km} km...")

    # 2. Iterate through each potential daily file
    for date in dates_to_check:
        daily_h5_path = os.path.join(psd_data_dir, f"{date.strftime('%Y-%m-%d')}.h5")
        
        if not os.path.exists(daily_h5_path):
            continue

        try:
            with h5py.File(daily_h5_path, 'r') as hf:
                # Read frequency vector if we haven't already
                if f_axis is None and 'frequency' in hf:
                    f_axis = hf['frequency'][:]
                
                group_name = f"ch_{channel_km:.1f}km"
                if group_name not in hf:
                    continue
                
                channel_group = hf[group_name]
                
                # 3. Iterate through datasets (timestamps) in the channel group
                for ts_str in channel_group.keys():
                    try:
                        ts = datetime.strptime(ts_str, '%Y%m%dT%H%M%S')
                        # 4. Filter by the precise time window
                        if time_start <= ts < time_end:
                            time_stamps.append(ts)
                            psd_list.append(channel_group[ts_str][:])
                    except ValueError:
                        continue # Ignore keys that are not valid timestamps
        except Exception as e:
            print(f"    Warning: Could not read file {daily_h5_path}. Error: {e}")

    # 5. Final processing and sorting
    if not psd_list:
        return [], None, None
    
    # Sort the data by timestamp, as it may be unordered from reading multiple files
    sorted_indices = np.argsort(time_stamps)
    sorted_times = [time_stamps[i] for i in sorted_indices]
    sorted_psds = [psd_list[i] for i in sorted_indices]
    
    # Transpose to get shape (n_frequencies, n_times)
    psd_matrix = np.array(sorted_psds).T
    
    print(f"  -> Loaded {psd_matrix.shape[1]} time steps for channel {channel_km} km.")
    return sorted_times, f_axis, psd_matrix

# ####################################################################################
# # --- 2. VISUALIZATION FUNCTIONS (Unchanged) ---
# ####################################################################################

def plot_spectrogram(cfg):
    """
    为指定通道和时间窗口绘制光谱图。
    This function now calls the new HDF5 data loader.
    """
    print("\n--- Plotting Spectrograms ---")
    t_start, t_end = cfg.zoom_window
    n_channels = len(cfg.channels_to_plot)

    if n_channels == 0:
        print("No channels specified for plotting.")
        return
    
    fig, axes = plt.subplots(n_channels, 1, figsize=(8, 4 * n_channels))
    fig.suptitle('Spectrogram', fontsize=20, y=0.99)
    
    if n_channels == 1:
        axes = [axes]
        
    for i, chan_km in enumerate(cfg.channels_to_plot):
        ax = axes[i]
        
        print(f"Processing Spectrogram for channel: {chan_km} km")
        
        # MODIFICATION: Calling the new data loading function
        times, f, Pxx_matrix = load_psd_data_from_daily_h5(cfg.psd_data_dir, chan_km, t_start, t_end)
        
        if not times:
            print(f"  No data found for channel {chan_km} km in the specified time zoom. Skipping.")
            ax.text(0.5, 0.5, 'No data in time zoom', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Channel {chan_km} km')
            continue
            
        # The rest of the plotting logic is identical
        asd_matrix = np.sqrt(Pxx_matrix)
        
        positive_asd = asd_matrix[asd_matrix > 0]
        if positive_asd.size > 0:
            vmin = np.percentile(positive_asd, 5)
            vmax = np.percentile(positive_asd, 95)
        else:
            vmin, vmax = 0.001, 1.0
        log_norm = colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)

        im = ax.pcolormesh(mdates.date2num(times), f, asd_matrix, 
                           shading='gouraud', cmap='viridis', norm=log_norm)

        ax.axhline(y=8, color='r', linestyle='--', linewidth=0.5)
        ax.axhline(y=20, color='r', linestyle='--', linewidth=0.5)

        ax.set_yscale('log')
        ax.set_ylim(cfg.freq_bottom, cfg.freq_top)
        ax.set_title(f'Channel {chan_km} km', fontsize=16)
        ax.set_ylabel('Frequency [Hz]', fontsize=12)
        
        if cfg.Day_tick > 0 or cfg.Hour_tick > 0:
            Ticks_xaxis_day_or_hour(ax, cfg.Day_tick, cfg.Hour_tick, rotation=30, font=10)
        if cfg.Minute_tick > 0 or cfg.Second_tick > 0:
            Ticks_xaxis_minute_or_second(ax, cfg.Minute_tick, cfg.Second_tick, rotation=0,font=14)
        
        ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.5)

        cbar = fig.colorbar(im, ax=ax, pad=0.01)
        cbar.set_label('ASD [Strain Rate/$\\sqrt{Hz}$]', fontsize=12)
        
        ax.set_xlim(cfg.zoom_window)

    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    
    time_str = f"_{t_start.strftime('%Y%m%dT%H%M')}_to_{t_end.strftime('%Y%m%dT%H%M')}"
    plot_filename = f'spectrogram_subplots{time_str}.png'
    fig.savefig(os.path.join(cfg.plot_output_dir, plot_filename))
    plt.close(fig)
        
    print("--- Spectrogram plots complete. ---")


def plot_binned_timeseries(cfg):
    """
    在一张图中以子图形式绘制指定通道的分频段时序图。
    This function now calls the new HDF5 data loader.
    """
    print("\n--- Plotting Binned Time-Series ---")
    t_start, t_end = cfg.zoom_window
    n_channels = len(cfg.channels_to_plot)
    if n_channels == 0:
        print("No channels specified for plotting.")
        return

    fig, axes = plt.subplots(n_channels, 1, figsize=(8, 3 * n_channels))
    fig.suptitle('Binned Mean ASD Variation', fontsize=20, y=0.99)
    
    if n_channels == 1:
        axes = [axes]
        
    for i, chan_km in enumerate(cfg.channels_to_plot):
        ax = axes[i]
        print(f"Processing Binned Time-Series for channel: {chan_km} km")
        
        # MODIFICATION: Calling the new data loading function
        times, f, Pxx_matrix = load_psd_data_from_daily_h5(cfg.psd_data_dir, chan_km, t_start, t_end)
        
        if not times:
            ax.text(0.5, 0.5, 'No data in time zoom', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Channel {chan_km} km')
            continue

        # The rest of the plotting logic is identical
        asd_matrix = np.sqrt(Pxx_matrix)
        
        for f_bin in cfg.frequency_bins:
            f_low, f_high = f_bin
            indices = np.where((f >= f_low) & (f < f_high))[0]
            if indices.size > 0:
                mean_asd_over_bin = np.mean(asd_matrix[indices, :], axis=0)
                label=f'{f_low}-{f_high} Hz'
                ax.plot(times, mean_asd_over_bin, marker=None, linestyle='-', markersize=3, linewidth= 1,label=label, alpha=0.5)
                # ax.scatter(times, mean_asd_over_bin, marker='.', markersize=4, label=label, alpha=0.5)
        
        ax.set_title(f'Channel {chan_km} km')
        ax.set_ylabel('Mean ASD [Strain Rate/$\\sqrt{Hz}$]')
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend()
        ax.set_yscale('log')
        
        ax.set_xlim(cfg.zoom_window)
        
        if cfg.Day_tick > 0 or cfg.Hour_tick > 0:
            Ticks_xaxis_day_or_hour(ax, cfg.Day_tick, cfg.Hour_tick, rotation=30, font=10)

        if cfg.Minute_tick > 0 or cfg.Second_tick > 0:
            Ticks_xaxis_minute_or_second(ax, cfg.Minute_tick, cfg.Second_tick, rotation=0,font=14)

    fig.autofmt_xdate(rotation=30, ha='right')
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])

    time_str = f"_{t_start.strftime('%Y%m%dT%H%M')}_to_{t_end.strftime('%Y%m%dT%H%M')}"
    plot_filename = f'binned_timeseries_subplots{time_str}.png'
    fig.savefig(os.path.join(cfg.plot_output_dir, plot_filename))
    plt.close(fig)
    print("--- Binned time-series plot complete. ---")

def plot_spectrogram_and_binned_timeseries(cfg): # ax参数不再需要从外部传入
    """
    Plot both spectrograms and binned time-series for a list of channels. Each save as a independent figure.
    This function is a wrapper to call the individual plotting functions.
    """
    print("\n--- Plotting Spectrograms and Binned Time-Series ---")
            
    # --- Font Sizes ---
    # 在这里集中调整所有绘图相关的字体大小
    title_fontsize = 29
    label_fontsize = 24
    cbar_label_fontsize = 24
    legend_fontsize = 20
    day_hour_tick_fontsize = 22
    min_sec_tick_fontsize = 22
    
    t_start, t_end = cfg.zoom_window
    
    for ii, chan_km in enumerate(cfg.channels_to_plot):
        
        fig = plt.figure(figsize=(22, 18))
        
        # [修正] 修正GridSpec以创建3个对齐的子图
        gs = plt.GridSpec(3, 2, figure=fig, height_ratios=[6, 3, 3], width_ratios=[20, 0.3], hspace=0.1)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
        ax2 = fig.add_subplot(gs[2, 0], sharex=ax0)
        cax = fig.add_subplot(gs[0, 1]) # Colorbar的独立坐标轴
        axes = [ax0, ax1, ax2]

        # --------------------------------------------------------------------
        # --- 子图1: Spectrogram ---
        # --------------------------------------------------------------------
        ax = axes[0]
        times, f, Pxx_matrix = load_psd_data_from_daily_h5(cfg.psd_data_dir, chan_km, t_start, t_end)
        
        if not times:
            print(f"  No data found for channel {chan_km} km in the specified time zoom. Skipping.")
            ax.text(0.5, 0.5, 'No data in time zoom', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'Channel {chan_km} km', fontsize=title_fontsize)
            plt.close(fig) # 关闭未完成的图
            continue
            
        asd_matrix = np.sqrt(Pxx_matrix)
        
        positive_asd = asd_matrix[asd_matrix > 0]
        vmin, vmax = (np.percentile(positive_asd, 5), np.percentile(positive_asd, 95)) if positive_asd.size > 0 else (0.001, 1.0)
        log_norm = colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)

        im = ax.pcolormesh(mdates.date2num(times), f, asd_matrix, shading='gouraud', cmap='viridis', norm=log_norm)

        # ax.axhline(y=8, color='r', linestyle='--', linewidth=1)
        # ax.axhline(y=20, color='r', linestyle='--', linewidth=1)

        ax.set_yscale('log')
        ax.set_ylim(cfg.freq_bottom, cfg.freq_top)
        ax.tick_params(axis='y', labelsize=day_hour_tick_fontsize)
        ax.set_title(f'Spectrogram at Channel {chan_km} km', fontsize=title_fontsize)
        ax.set_ylabel('Frequency [Hz]', fontsize=label_fontsize)
        ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.5)
        
        # 隐藏X轴标签，因为它与下面的图共享
        plt.setp(ax.get_xticklabels(), visible=False)

        cbar = fig.colorbar(im, cax=cax) # 将colorbar画在指定的cax上
        cbar.set_label('ASD [Strain Rate/$\\sqrt{Hz}$]', fontsize=cbar_label_fontsize)
        cbar.ax.tick_params(labelsize=day_hour_tick_fontsize)

        # --------------------------------------------------------------------
        # --- 子图2: Binned time-series & Significant Wave Height ---
        # --------------------------------------------------------------------
        ax = axes[1]
        
        # 绘制分频段时序数据
        asd_matrix = np.sqrt(Pxx_matrix)
        for ii, f_bin in enumerate(cfg.frequency_bins):
            f_low, f_high = f_bin
            indices = np.where((f >= f_low) & (f < f_high))[0]
            if indices.size > 0:
                mean_asd_over_bin = np.mean(asd_matrix[indices, :], axis=0)
                label=f'ASD average in {f_low}-{f_high} Hz'
                if ii ==0:
                    cc = 'b'
                elif ii ==1:
                    cc = 'g'
                ax.plot(times, mean_asd_over_bin, linestyle='-', linewidth=2, label=label, alpha=0.8, color=cc)
        
        ax.set_ylabel('Mean ASD', fontsize=label_fontsize)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend(loc='upper left', fontsize=legend_fontsize)
        ax.set_yscale('log')
        ax.tick_params(axis='y', labelsize=day_hour_tick_fontsize)
        plt.setp(ax.get_xticklabels(), visible=False) # 隐藏X轴标签

        # [新增] 加载气象数据
        atmos_df, wave_df = load_kefa_meteo_data(t_start, t_end)
        
        # [新增] 在右侧Y轴绘制有效波高 (SWH)
        ax_twin = ax.twinx()
        color = 'r'
        ax_twin.set_ylabel('Significant wave height (m)', color=color, fontsize=label_fontsize)
        ax_twin.plot(wave_df.index, wave_df['swh'], color=color, linestyle='--', linewidth=3, label='SWH')
        ax_twin.tick_params(axis='y', labelcolor=color, labelsize=day_hour_tick_fontsize)
        ax_twin.legend(loc='upper right', fontsize=legend_fontsize)
        
        # --------------------------------------------------------------------
        # --- 子图3: Wind Speed & Mean Wave Period ---
        # --------------------------------------------------------------------
        ax = axes[2]
        
        # [新增] 绘制风速
        color = 'C5'
        ax.set_ylabel('Wind Speed (m/s)', color=color, fontsize=label_fontsize)
        ax.plot(atmos_df.index, atmos_df['wind_speed'], color=color, linestyle='-', linewidth=2, label='Wind Speed')
        ax.tick_params(axis='y', labelcolor=color, labelsize=day_hour_tick_fontsize)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend(loc='upper left', fontsize=legend_fontsize)

        # [新增] 在右侧Y轴绘制平均波周期 (MWP)
        ax_twin = ax.twinx()
        color = 'tab:orange'
        ax_twin.set_ylabel('Mean Wave Period (s)', color=color, fontsize=label_fontsize)
        ax_twin.plot(wave_df.index, wave_df['mwp'], color=color, linestyle='--', linewidth=2, label='MWP')
        ax_twin.tick_params(axis='y', labelcolor=color, labelsize=day_hour_tick_fontsize)
        ax_twin.legend(loc='upper right', fontsize=legend_fontsize)

        # 为最下面的图设置X轴的时间刻度
        ax.set_xlim(cfg.zoom_window)
        if cfg.Day_tick > 0 or cfg.Hour_tick > 0:
            Ticks_xaxis_day_or_hour(ax, cfg.Day_tick, cfg.Hour_tick, rotation=30, font=day_hour_tick_fontsize)
        if cfg.Minute_tick > 0 or cfg.Second_tick > 0:
            Ticks_xaxis_minute_or_second(ax, cfg.Minute_tick, cfg.Second_tick, rotation=0, font=min_sec_tick_fontsize)
        ax.set_xlabel('Time', fontsize=label_fontsize)
        
        # --------------------------------------------------------------------
        # --- Final Adjustments ---
        # --------------------------------------------------------------------
        # plt.tight_layout(rect=[0, 0, 1, 0.98]) # 调整布局以适应主标题
        
        # 保存图像
        save_path = f"{cfg.plot_output_dir}/sg_meteo_channel_{chan_km}.png"
        plt.savefig(save_path)
        # print(f"  Figure saved to {save_path}")
        plt.show() # 用于在Jupyter或交互式环境中显示
        plt.close(fig) # 绘制完并保存后关闭图形，释放内存
        

    print("--- All plots complete. ---")
    
# ####################################################################################
# # --- 0. CONFIGURATION CLASS (Updated) ---
# ####################################################################################

class Config:
    """
    A class to hold all script parameters and save them to a JSON file.
    """
    def __init__(self):
        # --- 1. Locate the pre-computed PSD data ---
        # [去硬编码] 改为包内 noise_psd/ 自身目录 (相对包根);
        #   原: f"{HOME_dir}Run_4_amp_spec_final/"。PSD daily h5 默认在此目录下。
        self.output_dir_base = os.path.join(_PKG_ROOT, "noise_psd")
        
        # !! Must match the parameters from the calculation script !!
        self.nperseg_seconds = 200
        self.chunk_duration = timedelta(minutes=40)
        
        # Reconstruct the path to the PSD data directory
        params_str = f"chunk{self.chunk_duration}_nperseg{int(self.nperseg_seconds)}s"
        
        # !! CRITICAL: This path now points to the directory with daily H5 files !!
        self.psd_data_dir = os.path.join(self.output_dir_base, params_str, 'psd_daily_h5')


        # --- 3. Define Plotting Parameters ---
        # Select channels to plot
        self.channels_to_plot = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8, 8.3, 8.6,
                                          9, 9.2, 9.4, 9.6, 9.8, 10, 10.2, 10.4, 10.6, 10.8, 11, 11.2, 11.4, 11.6, 11.8, 12, 12.2, 12.4, 12.6, 12.8,
                                          13, 13.2, 13.4, 13.6, 13.8, 14,
                                          14.2, 14.4, 14.6, 14.8, 15, 15.2])

        # Define the time window for plotting
        # self.zoom_window = (
        #     datetime(2024, 5, 2, 12, 0, 0),
        #     datetime(2024, 5, 5, 12, 0, 0)
        # )
        self.zoom_window = (
            datetime.strptime('2024-04-24 10:20:00', '%Y-%m-%d %H:%M:%S'),
            datetime.strptime('2024-06-22 17:10:00', '%Y-%m-%d %H:%M:%S')
        )
        # self.zoom_window = (
        #     datetime.strptime('2024-05-03 04:40:00', '%Y-%m-%d %H:%M:%S'),
        #     datetime.strptime('2024-05-03 04:50:00', '%Y-%m-%d %H:%M:%S')
        # )
        
        
        self.zoom_tstr = f"{self.zoom_window[0].strftime('%Y%m%dT%H%M')}_to_{self.zoom_window[1].strftime('%Y%m%dT%H%M')}"
        
        # Define frequency bins for the time-series plot
        self.frequency_bins = [(0.02, 1), (4, 35)]
        # self.frequency_bins = [(1, 5), (5, 10), (10,20), (20, 35)]
        
        # Define the y-axis (frequency) range for the spectrogram
        self.freq_bottom, self.freq_top = 2/self.nperseg_seconds, 40.0
        
        # x-axis tick formatting
        self.Day_tick, self.Hour_tick = 5, 0
        self.Minute_tick, self.Second_tick = 0, 0

        # --- 2. Create the plot output directory ---
        # The 'plots' directory will be a sibling to the data directory
        self.plot_output_dir = os.path.join(os.path.dirname(self.psd_data_dir), 'plots',f'{self.zoom_tstr}')

    def _to_serializable_dict(self):
        """Converts object attributes to a JSON-serializable dictionary."""
        # This function remains unchanged
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
        """Saves the configuration to a JSON file in the plot output directory."""
        # This function remains unchanged
        os.makedirs(self.plot_output_dir, exist_ok=True)
        path = os.path.join(self.plot_output_dir, 'plot_config.json')
        serializable_dict = self._to_serializable_dict()
        
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(serializable_dict, f, indent=4, ensure_ascii=False)
        print(f"\nConfiguration saved to: {path}")


# ####################################################################################
# # --- 3. EXECUTION SCRIPT (Unchanged) ---
# ####################################################################################
if __name__ == '__main__':
    
    # 1. Instantiate Configuration
    cfg = Config()
    
    # 2. Create output directory
    os.makedirs(cfg.plot_output_dir, exist_ok=True)
    
    print(f"--- Starting Plotting Stage ---")
    print(f"Reading data from: {cfg.psd_data_dir}")
    print(f"Saving plots to:   {cfg.plot_output_dir}")

    # 3. Check for data and Call plotting functions
    if not os.path.isdir(cfg.psd_data_dir):
        print(f"Error: PSD data directory not found at '{cfg.psd_data_dir}'")
        print("Please run the 'calculate_psd_daily_hdf5' script first or check the path in the Config class.")
    else:
        # Plot spectrograms
        # plot_spectrogram(cfg)
        
        # # Plot binned time-series
        # plot_binned_timeseries(cfg)
        
        plot_spectrogram_and_binned_timeseries(cfg)
        # 4. Save configuration to JSON
        cfg.save_to_json()