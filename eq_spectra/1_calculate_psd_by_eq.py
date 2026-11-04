import pandas as pd
import numpy as np
import scipy.signal as signal
from datetime import datetime, timedelta
import os
import h5py
import json
from tqdm import tqdm
import gc

# ##################################################################
# # 重要：请确保您有一个名为 load_das_h5_PE2024.py 的文件在您的项目中，
# # 且该文件包含一个名为 load_das_custom 的函数。
# ##################################################################
import sys
# --- [去硬编码] 包根推导 + functions 路径 (相对 __file__) ---
# 包结构: check_minmag_pkg/{functions, eq_spectra, ...}; 本文件在 eq_spectra/ 下。
_PKG_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(_PKG_ROOT, "functions"))
import load_das_h5_PE2024 as load_das_h5_PE2024

# --- [显式可见] 外部大 DAS 数据目录 (事件 PSD 重算才需要; 换机改这里) ---
DAS_DATA_DIR = os.path.join(_PKG_ROOT, 'IN_DAS_data') + os.sep  # [GitHub 版] 默认 = 包内 IN_DAS_data/（放入或软链接 80 Hz 归档 DATA_DEC_80Hz_16m_SORTED/）；换机改这一行


class Config:
    """管理所有输入、输出和计算参数的配置类。"""
    def __init__(self):
        # --- 1. 基本路径设置 ---
        # [去硬编码] DAS 数据目录取顶层显式变量; 事件列表目录改指本包 eq_spectra/eq_list/。
        self.DAS_DATA_DIR = DAS_DATA_DIR

        # [去硬编码] 指向本包内随包附带的事件列表 (只含 Events_classified_Real_Mag.csv;
        # 重算产出的 psd_data_channel_*.h5 也会落在此目录下的 template 子目录)。
        # 默认事件列表 = RAWCCsMag_5-20Hz（与 step-2 出图脚本 2_plot_psd_by_eq.py 统一，端到端一致）。
        #   另一份 eq_list_MF_v3_40_ampMag_10-35Hz/ 也随包附带，切换改下面这行即可。
        self.EVENT_LIST_DIR = os.path.join(
            _PKG_ROOT, "eq_spectra", "eq_list", "eq_list_MF_v3_40_RAWCCsMag_5-20Hz")
        self.OUTPUT_BASE_DIR = self.EVENT_LIST_DIR

        # --- 2. 输入文件和模板设置 ---
        self.EVENT_CSV_FILE = os.path.join(self.EVENT_LIST_DIR, 'Events_classified_Real_Mag.csv')
        self.TEMPLATE_NAME = 'template_0423_80Hz'

        # --- 3. 通道和数据参数 ---
        self.SAMPLING_RATE_HZ = 80.0
        # self.CHANNELS_KM = [8.3, 9, 10, 11, 12, 13, 14.0] 
        self.CHANNELS_KM = [8.3, 8.6, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13,13.5, 14.0] 
        self.DX_METERS = 16.0 * 1.020951986312866
        self.CHANNELS_IDX = [int(km * 1000 / self.DX_METERS) for km in self.CHANNELS_KM]

        # --- 4. PSD计算参数 ---
        self.TIME_WINDOW_SECONDS = (1.5, 2.5) # 事件前后时间窗口, pre & post
        self.NPERSEG_SECONDS = 4

        # --- 5. 自动生成的路径 ---
        run_folder_name = f"{self.TEMPLATE_NAME}"
        self.run_output_dir = os.path.join(self.OUTPUT_BASE_DIR, run_folder_name)

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

    def save_to_json(self):
        os.makedirs(self.run_output_dir, exist_ok=True)
        path = os.path.join(self.run_output_dir, 'config.json')
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(self.to_serializable_dict(), f, indent=4, ensure_ascii=False)
        print(f"配置已保存至: {path}")

def main_calculation_optimized(config: Config):
    """
    优化后的计算函数。
    循环逻辑：先遍历事件，对每个事件读取一次数据，然后在内存中处理所有相关通道。
    """
    print("--- 开始PSD计算和存储 (优化流程) ---")
    try:
        events_df = pd.read_csv(config.EVENT_CSV_FILE, sep='\t', parse_dates=['Rt_trig'])
        events_df = events_df[events_df['Template'] == config.TEMPLATE_NAME]
        print(f"成功加载 {len(events_df)} 个事件，来源： {config.EVENT_CSV_FILE}")
    except FileNotFoundError:
        print(f"错误：事件文件未找到于 {config.EVENT_CSV_FILE}")
        return
        
    nperseg = int(config.NPERSEG_SECONDS * config.SAMPLING_RATE_HZ)
    freq_vector = np.fft.rfftfreq(nperseg, d=1.0/config.SAMPLING_RATE_HZ)

    # --- 高效文件管理：预先为每个通道创建并打开HDF5文件 ---
    h5_files = {}
    print("--- 正在初始化输出文件... ---")
    for i, chan_idx in enumerate(config.CHANNELS_IDX):
        chan_km = config.CHANNELS_KM[i]
        h5_path = os.path.join(config.run_output_dir, f'psd_data_channel_{chan_km:.1f}km.h5')
        
        # 以写入模式打开文件
        hf = h5py.File(h5_path, 'w')
        hf.create_dataset('frequency', data=freq_vector)
        hf.attrs['channel_km'] = chan_km
        hf.attrs['channel_idx'] = chan_idx
        hf.attrs['source_event_file'] = config.EVENT_CSV_FILE
        # 将打开的文件句柄存入字典，以便后续使用
        h5_files[chan_idx] = hf
    print("--- 输出文件初始化完成 ---")

    # --- 主循环：首先遍历事件 (Event-centric loop) ---
    print("--- 开始处理事件... ---")
    # 使用 iterrows() 配合 tqdm 来显示进度
    for _, event in tqdm(list(events_df.iterrows()), desc="处理事件"):
        try:
            t_trig = event['Rt_trig']
            magnitude = event['Real_Mag']
            
            t_start = t_trig - timedelta(seconds=config.TIME_WINDOW_SECONDS[0])
            t_end = t_trig + timedelta(seconds=config.TIME_WINDOW_SECONDS[1])
            
            # --- 关键优化：对每个事件，只读取一次数据 ---
            # load_das_custom 应该返回一个 (时间样本, 所有通道) 的数据块
            das_data, headers, _ = load_das_h5_PE2024.load_das_custom(
                t_start, t_end, input_dir=config.DAS_DATA_DIR
            )
            
            # 检查数据有效性
            if das_data is None or das_data.shape[0] < nperseg:
                print(f"警告：在 {t_trig} 跳过事件（数据不足或缺失）")
                continue

            fs_actual = headers.get('fs', config.SAMPLING_RATE_HZ)

            # --- 内层循环：处理这个事件对应的所有通道 ---
            # das_data 已经在内存中，这里不再有磁盘读取
            for i, chan_idx in enumerate(config.CHANNELS_IDX):
                chan_km = config.CHANNELS_KM[i]

                # 检查所需通道索引是否在加载的数据范围内
                if chan_idx >= das_data.shape[1]:
                    print(f"警告：通道索引 {chan_idx} (对应 {chan_km}km) 超出数据范围 (0-{das_data.shape[1]-1})。跳过此通道。")
                    continue
                
                # 从内存中提取单个通道的数据
                trace = das_data[:, chan_idx].astype(np.float64)

                # 计算PSD
                f, pxx = signal.welch(trace, fs=fs_actual, nperseg=nperseg, scaling='density')
                
                # 从字典中获取正确的文件句柄并写入数据
                hf = h5_files[chan_idx]
                dataset_name = f"event_{t_trig.strftime('%Y%m%dT%H%M%S')}"
                dset = hf.create_dataset(dataset_name, data=pxx, compression="gzip")
                dset.attrs['trigger_time'] = t_trig.isoformat()
                dset.attrs['magnitude'] = float(magnitude)

        except Exception as e:
            print(f"\n错误：处理事件 {event.get('Rt_trig')} 时失败: {e}")
        finally:
            # 释放大的数据块内存，并进行垃圾回收
            if 'das_data' in locals():
                del das_data
            gc.collect()

    # --- 收尾工作：关闭所有HDF5文件 ---
    print("\n--- 计算完成，正在关闭所有文件... ---")
    for chan_idx, hf in h5_files.items():
        try:
            hf.close()
        except Exception as e:
            print(f"关闭文件 (通道索引 {chan_idx}) 时出错: {e}")
    
    print("\n--- 所有操作成功完成 ---")


if __name__ == '__main__':
    config = Config()
    config.save_to_json()
    # 调用优化后的函数
    main_calculation_optimized(config)