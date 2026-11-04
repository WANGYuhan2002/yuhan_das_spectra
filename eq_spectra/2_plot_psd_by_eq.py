# -*- coding: utf-8 -*-
"""
advanced_psd_analysis_plot.py

VERSION: v2.5
DESCRIPTION:
This script creates a comprehensive PSD analysis plot on a per-channel basis.
For each specified channel, it combines two distinct datasets:
1. Long-term, continuous PSD data for that specific channel to generate a
   Probabilistic Power Spectral Density (PPSD) background, representing the
   ambient noise field.
2. Event-specific PSD data from the same channel, which is filtered by time and
   magnitude, and then overlaid as median spectral lines, grouped by magnitude.

This approach provides a clear visualization of specific event signals relative
to the statistical distribution of background noise for individual channels.

MODIFICATIONS (v2.5):
- ADDED: Inline text labels are now plotted next to each median spectral line to show the magnitude range directly on the chart.
- ADDED: A semi-transparent background is added to the inline labels for better readability.
- MODIFIED: The legend title is now more descriptive, explicitly stating that lines represent the "Median PSD".
- (from v2.4) ADDED: Amplitude Extraction Band parameters (frequency and PSD range) are now defined in the Config class.
- (from v2.4) ADDED: The minimum number of unique events required to plot a median spectrum is now a parameter in the Config class.
- (from v2.4) MODIFIED: The plotting function now uses the new Config parameters to draw the extraction band and check the event count.
- (from v2.3) MODIFIED: The main execution block now iterates through each channel defined in `CHANNELS_FOR_PPSD`.
- (from v2.3) MODIFIED: For each iteration, both event PSDs and continuous PPSD data are filtered/loaded for that single channel.
- (from v2.3) MODIFIED: A separate plot is generated for each channel.
- (from v2.3) MODIFIED: The plot title now includes the channel number and the time window used for the PPSD background calculation.
- (from v2.3) MODIFIED: The output filename is now unique for each channel.
"""
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.dates as mdates
import numpy as np
from datetime import datetime, timedelta
import os
import sys
import h5py
import json
import pandas as pd
from tqdm import tqdm

# ==============================================================================
# SECTION 0: CONFIGURATION CLASS
# 管理所有分析的路径和参数
# ==============================================================================
class Config:
    def __init__(self):
        # --- [去硬编码] 包根推导 (相对 __file__); 本脚本不 import functions ---
        # 包结构: check_minmag_pkg/{eq_spectra, ...}; 本文件在 eq_spectra/ 下。
        _pkg_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

        # --- 1. 连续PSD数据路径 (用于PPSD背景) ---
        # [注意] 这份按天连续 PSD daily h5 (数十 GB) 不随包附带; 需先用 noise_psd/
        #   或 Run_4 系列脚本对外部 DAS 重算才能生成。这里给一个包内默认占位路径,
        #   换成你重算所得目录即可。原指 check_minmag_validity/Run_4_amp_spec_more_channels/
        #   chunk0:10:00_nperseg40s/psd_daily_h5/。
        self.CONTINUOUS_PSD_DIR = os.path.join(
            _pkg_root, "noise_psd", "chunk0:10:00_nperseg40s", "psd_daily_h5")

        # --- 2. 事件PSD数据路径 (用于前景谱线) ---
        self.Template = 'template_0423_80Hz'
        # self.MagList_subdir = 'eq_list_MF_v3_40_ampMag_10-35Hz'
        self.MagList_subdir = 'eq_list_MF_v3_40_RAWCCsMag_5-20Hz'  # [MODIFIED] 使用新的震级列表目录

        # [去硬编码] 事件 PSD 目录改指本包 eq_spectra/eq_list/<MagList>/<Template>/。
        #   随包只带 eq_list/<MagList>/Events_classified_Real_Mag.csv;
        #   psd_data_channel_*.h5 需先跑 1_calculate_psd_by_eq*.py 重算生成。
        self.EVENT_PSD_DIR = os.path.join(
            _pkg_root, "eq_spectra", "eq_list", self.MagList_subdir, self.Template)


        # --- 3. 图像输出目录 ---
        self.PLOT_OUTPUT_DIR = os.path.join(self.EVENT_PSD_DIR, 'psd_plots_per_channel/') # 修改了输出目录以避免混淆

        # --- 背景噪声数据加载参数 ---
        # self.CONTINUOUS_TIME_WINDOW = (
        #     datetime(2024, 5, 2, 12, 0, 0),
        #     datetime(2024, 5, 5, 12, 0, 0)
        # )
        
        self.CONTINUOUS_TIME_WINDOW = (
            datetime.strptime('2024-04-24 10:20:00', '%Y-%m-%d %H:%M:%S'),
            datetime.strptime('2024-06-22 17:10:00', '%Y-%m-%d %H:%M:%S')
        )
        
        # --- [USER REQUEST] 定义要单独分析和绘图的通道列表 ---
        # self.CHANNELS_FOR_PPSD = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
        self.CHANNELS_FOR_PPSD = [8.3, 8.6, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0]

        # ======================================================================
        # --- 事件筛选参数 ---
        # ======================================================================
        self.EVENT_TIME_WINDOW = (
            datetime.strptime('2024-04-24 10:20:00', '%Y-%m-%d %H:%M:%S'),
            datetime.strptime('2024-06-22 17:10:00', '%Y-%m-%d %H:%M:%S')
        )
        
        self.EVENT_MAG_WINDOW = (-1.0, 5.0)

        # --- 绘图参数 ---
        self.MAGNITUDE_BINS = [
            (0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.0),
            (1.0, 1.2), (1.2, 1.4), (1.4, 1.6), (1.6, 2.0),
            (2.0, 2.5), (2.5, 3.0),
            (3.0, 3.5), (3.5, 4.0), (4.0, 4.5),
            (4.5, 5.0)
        ]

        # [NEW] 定义触发震级分箱绘图所需的最小独立事件数
        self.MIN_EVENTS_FOR_MEDIAN_PLOT = 3

        # [NEW] 振幅提取带参数 (用于绘图)
        self.AMP_EXTRACT_FREQ_RANGE = (8, 20) # [Hz]
        self.AMP_EXTRACT_PSD_RANGE = (-200, -110) # [dB]

    def to_serializable_dict(self):
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
        os.makedirs(self.PLOT_OUTPUT_DIR, exist_ok=True)
        path = os.path.join(self.PLOT_OUTPUT_DIR, 'config.json')
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(self.to_serializable_dict(), f, indent=4, ensure_ascii=False)
        print(f"配置已保存至: {path}")

# ==============================================================================
# SECTION 1: 数据加载与处理函数
# ==============================================================================

# --- Function 1.1: 加载并筛选事件PSD数据 ---
def load_event_psds(event_psd_dir: str, time_window: tuple, mag_window: tuple):
    print("\n--- 正在加载并筛选所有通道的事件PSD数据 ---")
    print(f"时间筛选范围: {time_window[0].strftime('%Y-%m-%d %H:%M:%S')} 至 {time_window[1].strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"震级筛选范围: {mag_window[0]} 至 {mag_window[1]}")
    psd_data_list = []
    freq_vector = None
    if not os.path.isdir(event_psd_dir):
        print(f"错误: 事件PSD目录 '{event_psd_dir}' 未找到")
        return None, None
    psd_h5_files = [f for f in os.listdir(event_psd_dir) if f.startswith('psd_data_channel_') and f.endswith('.h5')]
    if not psd_h5_files:
        print(f"错误: 在 '{event_psd_dir}' 中未找到 'psd_data_channel_*.h5' 文件")
        return None, None
    unique_event_ids = set()
    for h5_filename in tqdm(psd_h5_files, desc="扫描所有事件PSD文件"):
        h5_path = os.path.join(event_psd_dir, h5_filename)
        with h5py.File(h5_path, 'r') as hf:
            if freq_vector is None:
                freq_vector = hf['frequency'][:]
            channel_km_val = hf.attrs['channel_km']
            for item_name in hf:
                if item_name.startswith('event_'):
                    dataset = hf[item_name]
                    magnitude = dataset.attrs['magnitude']
                    trigger_time_str = dataset.attrs['trigger_time']
                    trigger_time = datetime.fromisoformat(trigger_time_str)
                    if not (time_window[0] <= trigger_time <= time_window[1]):
                        continue
                    if not (mag_window[0] <= magnitude <= mag_window[1]):
                        continue
                    psd_data_list.append({
                        'magnitude': magnitude,
                        'psd': dataset[:],
                        'channel': channel_km_val,
                        'trigger_time': trigger_time
                    })
                    unique_event_ids.add(item_name)
    total_events_after_filter = len(unique_event_ids)
    if not psd_data_list:
        print("警告: 没有任何事件PSD数据通过筛选条件。")
        return pd.DataFrame(), None
    print(f"筛选后，共找到并加载了 {total_events_after_filter} 个独立事件记录 (跨所有通道)。")
    return pd.DataFrame(psd_data_list), freq_vector

# --- Function 1.2: 加载连续PSD数据 (用于PPSD背景) ---
def load_continuous_psds_for_ppsd(psd_data_dir, channels_to_load, time_start, time_end):
    print(f"\n--- 正在为通道 {channels_to_load} 加载连续PSD数据用于PPSD背景计算 ---")
    all_psds_list = []
    f_axis = None
    if not os.path.isdir(psd_data_dir):
        print(f"错误: 连续PSD目录 '{psd_data_dir}' 未找到")
        return None, None
    dates_to_check = [time_start.date() + timedelta(days=d) for d in range((time_end.date() - time_start.date()).days + 1)]
    # 此函数现在预期 channels_to_load 是一个包含单个通道的列表
    for channel_km in tqdm(channels_to_load, desc=f"加载通道 {channels_to_load[0]}km 数据"):
        channel_psd_list = []
        for date in dates_to_check:
            daily_h5_path = os.path.join(psd_data_dir, f"{date.strftime('%Y-%m-%d')}.h5")
            if not os.path.exists(daily_h5_path): continue
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
                                channel_psd_list.append(hf[group_name][ts_str][:])
                        except (ValueError, KeyError): continue
            except Exception as e:
                print(f"警告: 无法读取 {daily_h5_path} (通道 {channel_km})。错误: {e}")
        if channel_psd_list:
            all_psds_list.extend(channel_psd_list)
    if not all_psds_list:
        print(f"警告: 未能为通道 {channels_to_load[0]} 加载任何用于PPSD背景的连续PSD数据。")
        return None, None
    psd_matrix = np.array(all_psds_list)
    print(f"为通道 {channels_to_load[0]} 的PPSD计算加载了 {psd_matrix.shape[0]} 条连续PSD记录。")
    return psd_matrix, f_axis

# --- Function 1.3: 插值PSD矩阵 ---
def interpolate_psd_matrix(psd_matrix, original_freq, target_freq):
    print(f"\n--- 正在插值连续PSD数据，从 {len(original_freq)} 个频率点到 {len(target_freq)} 个频率点 ---")
    # 避免 log(0)
    valid_orig_idx = original_freq > 0
    valid_targ_idx = target_freq > 0
    log_original_freq = np.log(original_freq[valid_orig_idx])
    log_target_freq = np.log(target_freq[valid_targ_idx])

    interpolated_matrix = np.zeros((psd_matrix.shape[0], len(target_freq)))

    for i in tqdm(range(psd_matrix.shape[0]), desc="插值PSD"):
        log_psd = np.log(psd_matrix[i, valid_orig_idx])
        # 使用左右边界值进行填充，以处理目标频率超出原始频率范围的情况
        interpolated_log_psd = np.interp(log_target_freq, log_original_freq, log_psd, left=log_psd[0], right=log_psd[-1])
        interpolated_matrix[i, valid_targ_idx] = np.exp(interpolated_log_psd)
        # 如果目标频率从0开始，而原始频率不是，则需要处理第一个点
        if target_freq[0] == 0:
             # 可以选择将第一个点设置为插值后的第一个有效点的值，或者保持为0
             if original_freq[0] != 0:
                 interpolated_matrix[i, 0] = interpolated_matrix[i, 1]
             else:
                 interpolated_matrix[i, 0] = psd_matrix[i, 0]

    print("插值完成。")
    return interpolated_matrix

# ==============================================================================
# SECTION 2: 最终绘图函数
# ==============================================================================
def plot_combined_psd_analysis(event_df, common_freq, continuous_psd_matrix, channel_km: float, cfg: Config):
    """
    生成最终的分析图。事件和连续数据现在共享一个公共频率轴。
    所有图形文本均为英文。
    """
    # --- 新增：字体大小集中控制 ---
    font_size_base = 20  # 基础字号，可调整
    font_size_title = font_size_base + 2       # 图表标题
    font_size_label = font_size_base           # 坐标轴标签和颜色条标签
    font_size_legend_title = font_size_base    # 图例标题
    font_size_legend = font_size_base      # 图例内容
    font_size_annotation = font_size_base-2  # 注释文本 (如图中谱线旁的标签)
    font_size_tick = font_size_base - 2        # 刻度标签

    print(f"\n--- 正在为通道 {channel_km}km 生成组合PSD图 ---")
    fig, ax = plt.subplots(figsize=(14, 10))
    colors_ = plt.cm.plasma(np.linspace(0, 1, len(cfg.MAGNITUDE_BINS)))
    if common_freq is None or len(common_freq) < 2:
        print("错误: 公共频率向量无效，无法继续绘图。")
        return
    valid_freq_idx = common_freq > 0

    # --- PPSD 背景 (来自插值后的连续数据) ---
    if continuous_psd_matrix is not None and continuous_psd_matrix.size > 0:
        print("正在计算连续数据的PPSD...")
        all_psds_db = 10 * np.log10(continuous_psd_matrix[:, valid_freq_idx])
        db_bins = np.linspace(-220, -80, 141)
        # 确保频率分箱从有效频率的最小值开始
        freq_bins_log = np.logspace(np.log10(common_freq[valid_freq_idx].min()), np.log10(common_freq.max()), 100)
        hist2d = np.zeros((len(db_bins)-1, len(freq_bins_log)-1))

        for i in tqdm(range(len(freq_bins_log)-1), desc="计算PPSD"):
            f_start, f_end = freq_bins_log[i], freq_bins_log[i+1]
            # 在有效频率内查找索引
            freq_indices_in_valid = np.where((common_freq[valid_freq_idx] >= f_start) & (common_freq[valid_freq_idx] < f_end))[0]
            if not freq_indices_in_valid.any(): continue
            psd_in_f_bin = all_psds_db[:, freq_indices_in_valid].flatten()
            hist, _ = np.histogram(psd_in_f_bin, bins=db_bins)
            hist2d[:, i] = hist

        total_counts = hist2d.sum(axis=0) + 1e-9 # 避免除以零
        ppsd = (hist2d / total_counts) * 100

        X, Y = np.meshgrid(freq_bins_log, db_bins)
        im = ax.pcolormesh(X, Y, ppsd, cmap='YlGnBu', shading='auto', vmax=np.percentile(ppsd[ppsd>0], 95))
        cbar = fig.colorbar(im, ax=ax, pad=0.01)
        cbar.set_label('PSD Probability Density [%]', size=font_size_label)
        cbar.ax.tick_params(labelsize=font_size_tick)
    else:
        print("由于没有加载连续数据，跳过PPSD背景绘制。")

    # --- 中位谱线 (来自筛选后的事件数据) ---
    print("正在计算并绘制筛选后事件的中位谱线...")
    plotted_lines = 0
    for i, (m_min, m_max) in enumerate(cfg.MAGNITUDE_BINS):
        bin_df = event_df[(event_df['magnitude'] >= m_min) & (event_df['magnitude'] < m_max)]
        if bin_df.empty: continue

        num_unique_events = bin_df['trigger_time'].nunique()
        # [MODIFIED] 使用Config中定义的变量
        if num_unique_events < cfg.MIN_EVENTS_FOR_MEDIAN_PLOT:
            print(f"跳过震级分箱 M {m_min:.1f}-{m_max:.1f}，独立事件数量不足 ({num_unique_events} < {cfg.MIN_EVENTS_FOR_MEDIAN_PLOT})。")
            continue

        median_psd_db = 10 * np.log10(np.median(np.vstack(bin_df['psd'].values), axis=0)) #! choose median
        # median_psd_db = 10 * np.log10(np.mean(np.vstack(bin_df['psd'].values), axis=0)) #! choose average
        
        legend_label = f'M {m_min:.1f} to {m_max:.1f} ({num_unique_events})'
        ax.plot(common_freq, median_psd_db, color=colors_[i], lw=2.5, label=legend_label, zorder=10)
        

        # --- [NEW] 在谱线旁边添加文字说明震级范围 ---
        # 通过交错x轴位置来避免标签重叠
        # text_x_pos = 22 + (plotted_lines % 4) * 4  # Stagger between 22, 26, 30, 34 Hz
        text_x_pos = 4 + plotted_lines * 0.1  # 每条线的x位置间隔2Hz
        if text_x_pos > common_freq.max(): # 确保文本不会超出绘图范围
            text_x_pos = common_freq.max() * 0.8
            
        # 插入以找到曲线上我们所选x点的y值
        text_y_pos = np.interp(text_x_pos, common_freq[valid_freq_idx], median_psd_db[valid_freq_idx])

        # 添加文本标签，颜色与线条匹配
        inline_label = f'M {m_min:.1f}-{m_max:.1f}'
        if m_min > 1.4:
            ax.text(text_x_pos, text_y_pos, inline_label,
                    color=colors_[i],
                    fontsize=font_size_annotation,
                    weight='bold',
                    ha='left',
                    va='center', # 在线上垂直居中
                    bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', boxstyle='round,pad=0.1'), # 添加微弱的白色背景以提高可读性
                    zorder=20)
            
        plotted_lines += 1


    # --- [MODIFIED] 样式和注释: 绘制振幅提取带 (全部使用英文) ---
    # 从配置中获取参数
    freq_min, freq_max = cfg.AMP_EXTRACT_FREQ_RANGE
    psd_min, psd_max = cfg.AMP_EXTRACT_PSD_RANGE
    rect_width = freq_max - freq_min
    rect_height = psd_max - psd_min

    # 创建并添加矩形
    amp_rect = patches.Rectangle(
        (freq_min, psd_min), rect_width, rect_height,
        lw=1.5, edgecolor='r', facecolor='r', alpha=0.07, zorder=5
    )
    ax.add_patch(amp_rect)

    # 计算文本位置并添加 (确保文本在矩形上方且清晰可见)
    text_x = (freq_min + freq_max) / 2
    text_y = psd_max + 5  # 将文本垂直位置设置在矩形框上方5个单位处
    text_label = f'Amplitude Extraction Band\n({freq_min}-{freq_max} Hz)'
    ax.text(text_x, text_y, text_label, color='r', ha='center', va='bottom', fontsize=font_size_annotation, weight='bold', zorder=15)


    # --- [USER REQUEST] 动态生成包含通道和PPSD时间信息的英文标题 ---
    ppsd_time_str = f"PPSD Background from: {cfg.CONTINUOUS_TIME_WINDOW[0].strftime('%Y-%m-%d')} to {cfg.CONTINUOUS_TIME_WINDOW[1].strftime('%Y-%m-%d')}"
    plot_title = f'PSD Analysis of earthquakes for Channel {channel_km:.1f} km\n{ppsd_time_str}'

    ax.set_xscale('log')
    ax.set_xlabel('Frequency [Hz]', fontsize=font_size_label)
    ax.set_ylabel(r'PSD [dB, rel. to 1 $(strain/s)^2/Hz$]', fontsize=font_size_label)
    ax.set_title(plot_title, fontsize=font_size_title)
    ax.tick_params(axis='both', which='major', labelsize=font_size_tick)


    # 确保xlim从有效频率的最小值开始
    ax.set_xlim(common_freq[valid_freq_idx].min(), common_freq.max())
    ax.set_ylim(-210, -70) # 调整Y轴范围以更好地显示
    ax.grid(True, which='both', linestyle=':', linewidth=0.5)
    
    # --- [MODIFIED] 更新图例标题以明确说明是中位数 ---
    legend_title = f'Median PSD by Magnitude Bin\n(>{cfg.MIN_EVENTS_FOR_MEDIAN_PLOT-1} Events per bin)'
    legend = ax.legend(title=legend_title, fontsize=font_size_legend)
    plt.setp(legend.get_title(), fontsize=font_size_legend_title)

    fig.tight_layout()

    # --- [USER REQUEST] 生成包含通道信息的唯一文件名 ---
    time_str = f"{cfg.EVENT_TIME_WINDOW[0].strftime('%Y%m%d')}-{cfg.EVENT_TIME_WINDOW[1].strftime('%Y%m%d')}"
    mag_str = f"M{cfg.EVENT_MAG_WINDOW[0]:.1f}_{cfg.EVENT_MAG_WINDOW[1]:.1f}".replace('.', 'p')
    channel_str = f"ch{channel_km:.1f}km".replace('.', 'p')
    output_filename = os.path.join(cfg.PLOT_OUTPUT_DIR, f'psd_analysis_{channel_str}_{time_str}_{mag_str}.png')

    fig.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"\n✅ 通道 {channel_km}km 的分析图已保存至: {output_filename}")
    plt.close(fig)

# ==============================================================================
# SECTION 3: 主执行块
# ==============================================================================
if __name__ == '__main__':
    # 1. 初始化配置
    cfg = Config()
    cfg.save_config()

    # 2. 加载一次所有通道的事件数据以提高效率
    master_event_df, event_freq = load_event_psds(
        cfg.EVENT_PSD_DIR,
        cfg.EVENT_TIME_WINDOW,
        cfg.EVENT_MAG_WINDOW
    )

    if master_event_df is None or master_event_df.empty:
        print("\n严重错误: 无法加载任何事件数据。正在中止脚本。")
        sys.exit()
    if event_freq is None:
        print("\n严重错误: 未能从事件数据中加载频率向量。正在中止脚本。")
        sys.exit()

    # --- [USER REQUEST] 主循环：为每个指定的通道生成一个图 ---
    for channel in cfg.CHANNELS_FOR_PPSD:
        print(f"\n{'='*25} 开始处理通道 {channel} km {'='*25}")

        # 3a. 从主数据框中筛选当前通道的事件数据
        channel_event_df = master_event_df[master_event_df['channel'] == channel]
        if channel_event_df.empty:
            print(f"警告: 通道 {channel}km 没有找到符合条件的事件数据，跳过此通道。")
            continue
        print(f"通道 {channel}km: 找到 {len(channel_event_df)} 条事件记录。")

        # 3b. 仅为当前通道加载背景噪声数据
        continuous_psd_matrix, continuous_freq = load_continuous_psds_for_ppsd(
            cfg.CONTINUOUS_PSD_DIR,
            [channel], # 传递单个通道的列表
            cfg.CONTINUOUS_TIME_WINDOW[0],
            cfg.CONTINUOUS_TIME_WINDOW[1]
        )

        # 4. 验证数据，处理频率不匹配问题
        if continuous_psd_matrix is None or continuous_freq is None:
            print(f"警告: 未能加载通道 {channel}km 的连续背景数据，跳过此通道。")
            continue

        # --- 核心逻辑：处理频率向量不匹配 ---
        common_freq = event_freq
        interpolated_continuous_psd_matrix = continuous_psd_matrix

        if not np.array_equal(event_freq, continuous_freq):
            print(f"\n警告: 通道 {channel}km 的事件({len(event_freq)}pts)和连续数据({len(continuous_freq)}pts)频率向量不匹配。正在对连续数据进行插值...")
            interpolated_continuous_psd_matrix = interpolate_psd_matrix(
                continuous_psd_matrix,
                original_freq=continuous_freq,
                target_freq=common_freq
            )
        else:
            print(f"\n信息: 通道 {channel}km 的事件和连续数据频率向量已匹配，无需插值。")

        # 5. 所有数据准备就绪，为当前通道绘图
        plot_combined_psd_analysis(
            channel_event_df,
            common_freq,
            interpolated_continuous_psd_matrix,
            channel, # 传递当前通道号
            cfg
        )

    print("\n所有指定通道处理完毕。")