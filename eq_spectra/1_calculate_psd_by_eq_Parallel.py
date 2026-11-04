import pandas as pd
import numpy as np
import scipy.signal as signal
from datetime import datetime, timedelta
import os
import h5py
import json
from tqdm import tqdm
import gc
import multiprocessing
import traceback

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
        # --- 0. 并行计算设置 ---
        self.NUM_WORKERS = 10  # 设置工作的CPU核心数量

        # --- 1. 基本路径设置 ---
        # [去硬编码] DAS 数据目录取顶层显式变量; 事件列表目录改指本包 eq_spectra/eq_list/。
        self.DAS_DATA_DIR = DAS_DATA_DIR

        # 默认事件列表 = RAWCCsMag_5-20Hz（与 step-2 出图脚本统一，端到端一致）。
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
        # self.CHANNELS_KM = [8.3, 8.6, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14.0] 
        self.CHANNELS_KM = np.array([1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,   8.3, 8.6, 9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14.0])

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

def process_single_event(args):
    """
    由单个工作进程执行的函数，处理单个事件。
    包含了数据加载、PSD计算和加锁文件写入。
    """
    event_info, config, locks = args
    t_trig, magnitude = event_info

    try:
        t_start = t_trig - timedelta(seconds=config.TIME_WINDOW_SECONDS[0])
        t_end = t_trig + timedelta(seconds=config.TIME_WINDOW_SECONDS[1])
        nperseg = int(config.NPERSEG_SECONDS * config.SAMPLING_RATE_HZ)

        # 1. 加载数据
        das_data, headers, _ = load_das_h5_PE2024.load_das_custom(
            t_start, t_end, input_dir=config.DAS_DATA_DIR
        )

        if das_data is None or das_data.shape[0] < nperseg:
            return f"Skipped (data insufficient or missing) for event at {t_trig}"

        fs_actual = headers.get('fs', config.SAMPLING_RATE_HZ)

        # 2. 遍历通道，计算PSD并写入
        for i, chan_idx in enumerate(config.CHANNELS_IDX):
            chan_km = config.CHANNELS_KM[i]

            if chan_idx >= das_data.shape[1]:
                # 此警告仅在控制台打印，不中断流程
                print(f"警告：通道索引 {chan_idx} (对应 {chan_km}km) 超出数据范围 (0-{das_data.shape[1]-1})。跳过此通道。")
                continue
            
            trace = das_data[:, chan_idx].astype(np.float64)
            f, pxx = signal.welch(trace, fs=fs_actual, nperseg=nperseg, scaling='density')
            
            dataset_name = f"event_{t_trig.strftime('%Y%m%dT%H%M%S')}"
            
            # 3. 加锁写入HDF5文件，确保线程安全
            h5_path = os.path.join(config.run_output_dir, f'psd_data_channel_{chan_km:.1f}km.h5')
            
            # 从共享字典中获取此通道对应的锁
            lock = locks[chan_idx]
            with lock:
                # 以追加模式打开文件，写入数据集
                with h5py.File(h5_path, 'a') as hf:
                    dset = hf.create_dataset(dataset_name, data=pxx, compression="gzip")
                    dset.attrs['trigger_time'] = t_trig.isoformat()
                    dset.attrs['magnitude'] = float(magnitude)

        del das_data # 释放内存
        gc.collect()
        return True # 表示成功

    except Exception as e:
        # 捕获详细的异常信息并返回
        error_info = traceback.format_exc()
        return f"Failed to process event at {t_trig}: {e}\n{error_info}"


def main_parallel_computation(config: Config):
    """
    使用多进程并行处理事件的主函数。
    """
    print("--- 开始并行PSD计算 ---")
    
    # 1. 加载事件列表
    try:
        events_df = pd.read_csv(config.EVENT_CSV_FILE, sep='\t', parse_dates=['Rt_trig'])
        events_df = events_df[events_df['Template'] == config.TEMPLATE_NAME]
        print(f"成功加载 {len(events_df)} 个事件进行处理，来源： {config.EVENT_CSV_FILE}")
    except FileNotFoundError:
        print(f"错误：事件文件未找到于 {config.EVENT_CSV_FILE}")
        return

    # 2. 在主进程中预先创建HDF5文件和元数据，避免并行创建时的竞争条件
    print("--- 正在初始化输出文件... ---")
    nperseg = int(config.NPERSEG_SECONDS * config.SAMPLING_RATE_HZ)
    freq_vector = np.fft.rfftfreq(nperseg, d=1.0/config.SAMPLING_RATE_HZ)
    os.makedirs(config.run_output_dir, exist_ok=True)

    for i, chan_idx in enumerate(config.CHANNELS_IDX):
        chan_km = config.CHANNELS_KM[i]
        h5_path = os.path.join(config.run_output_dir, f'psd_data_channel_{chan_km:.1f}km.h5')
        with h5py.File(h5_path, 'w') as hf: # 'w'模式会清空旧文件
            hf.create_dataset('frequency', data=freq_vector)
            hf.attrs['channel_km'] = chan_km
            hf.attrs['channel_idx'] = chan_idx
            hf.attrs['source_event_file'] = config.EVENT_CSV_FILE
    print("--- 输出文件初始化完成 ---")

    # 3. 设置多进程管理器和锁
    # Manager用于在进程之间共享对象，这里我们用它来共享每个文件的锁
    manager = multiprocessing.Manager()
    locks = manager.dict()
    for chan_idx in config.CHANNELS_IDX:
        locks[chan_idx] = manager.Lock()

    # 4. 准备要传递给每个工作进程的参数列表
    tasks = [
        ((row['Rt_trig'], row['Real_Mag']), config, locks) 
        for _, row in events_df.iterrows()
    ]

    # 5. 创建并运行进程池
    print(f"--- 使用 {config.NUM_WORKERS} 个CPU核心开始处理事件, 来自模版 {config.TEMPLATE_NAME} ---")
    results = []
    
    # 使用 with 语句确保进程池被正确关闭
    # 使用 imap_unordered 以便在任务完成时立即更新进度条
    with multiprocessing.Pool(processes=config.NUM_WORKERS) as pool:
        with tqdm(total=len(tasks), desc="处理事件进度") as pbar:
            for result in pool.imap_unordered(process_single_event, tasks):
                results.append(result)
                pbar.update(1)

    # 6. 收尾工作和结果报告
    print("\n--- 所有任务已完成，正在分析结果... ---")
    
    success_count = sum(1 for r in results if r is True)
    failed_tasks = [r for r in results if r is not True]

    print(f"--- 操作总结 ---")
    print(f"总任务数: {len(tasks)}")
    print(f"成功处理: {success_count}")
    print(f"失败或跳过: {len(failed_tasks)}")

    if failed_tasks:
        print("\n--- 失败/跳过任务详情 ---")
        for i, failure_msg in enumerate(failed_tasks):
            print(f"{i+1}: {failure_msg}")
    
    print("\n--- 并行计算流程全部结束 ---")


if __name__ == '__main__':
    # 设置多进程启动方法为 'fork'，在Linux/macOS上通常是默认且最高效的
    # 在某些环境下（如macOS的新版本或Windows）可能需要显式设置
    try:
        multiprocessing.set_start_method('fork')
    except RuntimeError:
        print("多进程启动方法已设置，继续执行。")

    config = Config()
    config.save_to_json()
    
    # 调用并行计算函数
    main_parallel_computation(config)