# -----------------------------------------------------------------------------
# 脚本名称: 数据探索、提取与可视化工具 (最终版)
# 脚本功能:
# 1. 加载并探索业务化NetCDF数据。
# 2. 合并数据，并定义函数提取指定时空点的时间序列。
# 3. 新增强大的可视化函数，将提取的数据绘制成图表。
#
# --- 改动说明 ---
# 1. 增加了 matplotlib 库的导入。
# 2. 新增了第四部分和 plot_data 函数，用于数据可视化。
# 3. 图表严格按照要求，所有标签和标题使用英文。
# -----------------------------------------------------------------------------

# --- 导入所需库 ---
import xarray as xr
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# =============================================================================
# 第一部分: 数据加载与初步探索
# (此部分与之前相同，保持不变)
# =============================================================================
print("===== 正在加载并探索数据... =====")

data_dir = "kefa_wind"
accum_file = os.path.join(data_dir, "data_stream-oper_stepType-accum.nc")
instant_file = os.path.join(data_dir, "data_stream-oper_stepType-instant.nc")
wave_file = os.path.join(data_dir, "data_stream-wave_stepType-instant.nc")

try:
    ds_accum = xr.open_dataset(accum_file)
    ds_instant = xr.open_dataset(instant_file)
    ds_wave = xr.open_dataset(wave_file)
    print("--- 所有数据文件加载成功！ ---")
except Exception as e:
    print(f"加载数据时发生错误: {e}")
    # 在出错时退出，防止后续代码执行
    exit()

# =============================================================================
# 第二部分: 数据预处理与函数定义
# (此部分与之前相同，保持不变)
# =============================================================================

def merge_and_prepare_atmos_data(ds_inst, ds_acc):
    print("\n===== 正在合并大气数据并计算风速风向... =====")
    atmos_ds = xr.merge([ds_inst, ds_acc])
    print("--- 大气数据合并完成 ---")

    if 'u10' in atmos_ds and 'v10' in atmos_ds:
        print("--- 正在计算风速和风向... ---")
        wind_speed = np.sqrt(atmos_ds['u10']**2 + atmos_ds['v10']**2)
        wind_direction = (270 - np.rad2deg(np.arctan2(atmos_ds['v10'], atmos_ds['u10']))) % 360
        atmos_ds = atmos_ds.assign(wind_speed=wind_speed, wind_direction=wind_direction)
        atmos_ds.wind_speed.attrs['long_name'] = '10m wind speed'
        atmos_ds.wind_speed.attrs['units'] = 'm/s'
        atmos_ds.wind_direction.attrs['long_name'] = '10m wind direction (from)'
        atmos_ds.wind_direction.attrs['units'] = 'degree'
        print("--- 风速和风向计算并添加完成 ---")
    else:
        print("[警告] 数据集中未找到 'u10' 或 'v10'，无法计算风速风向。")
    
    return atmos_ds

def extract_point_data(ds, lat, lon, start_time, end_time, variables=None):
    if variables is None:
        variables_to_extract = list(ds.data_vars)
    else:
        if isinstance(variables, str):
            variables = [variables]
        missing_vars = [var for var in variables if var not in ds.data_vars]
        if missing_vars:
            print(f"[错误] 请求的变量 {missing_vars} 在数据集中不存在！可用变量: {list(ds.data_vars)}")
            return pd.DataFrame()
        variables_to_extract = variables

    print(f"\n===== 准备为地点 (Lat≈{lat}, Lon≈{lon}) 提取变量: {variables_to_extract} =====")
    print(f"--- 时间范围: {start_time} 到 {end_time} ---")

    try:
        point_selection = ds.sel(latitude=lat, longitude=lon, method='nearest')
        actual_lat = point_selection.latitude.item()
        actual_lon = point_selection.longitude.item()
        print(f"--- 诊断: 请求的经纬度 ({lat}, {lon}) -> 找到最接近的格点 ({actual_lat:.4f}, {actual_lon:.4f}) ---")
        
        time_slice = slice(start_time, end_time)
        data_selection = point_selection[variables_to_extract].sel(valid_time=time_slice)

        if data_selection.valid_time.size == 0:
            print("--- [警告] 在指定的时间范围内没有找到任何数据点！---")
            return pd.DataFrame()

        df = data_selection.to_dataframe()
        df = df.reset_index()
        df = df.drop(columns=[col for col in ['latitude', 'longitude', 'number', 'expver', 'step', 'surface'] if col in df.columns], errors='ignore')
        df = df.set_index('valid_time')
        print(f"--- [成功] 成功提取到 {len(df)} 条数据记录 ---")
        return df
    except Exception as e:
        print(f"[错误] 提取数据时发生未知错误: {e}")
        return pd.DataFrame()

# =============================================================================
# 第三部分: 示例：提取特定时空点的数据
# (此部分与之前相同，仅打印DataFrame的头部信息)
# =============================================================================

atmos_data = merge_and_prepare_atmos_data(ds_instant, ds_accum)

target_lat = 38.5
target_lon = 20.5
start_date = '2024-04-24 10:20:00' 
end_date = '2024-06-22 17:10:00'


atmos_vars_to_get = ['sst', 'msl', 'wind_speed', 'wind_direction']
atmos_df = extract_point_data(atmos_data, target_lat, target_lon, start_date, end_date, variables=atmos_vars_to_get)

if not atmos_df.empty:
    print("\n--- [输出] 提取到的大气数据 (DataFrame) ---")
    print(atmos_df.head())
else:
    print("\n--- [输出] 未能提取到大气数据 ---")

wave_vars_to_get = ['swh', 'mwp']
wave_df = extract_point_data(ds_wave, target_lat, target_lon, start_date, end_date, variables=wave_vars_to_get)

if not wave_df.empty:
    print("\n--- [输出] 提取到的海浪数据 (DataFrame) ---")
    print(wave_df.head())
else:
    print("\n--- [输出] 未能提取到海浪数据 ---")
    
    
#! 保存 pd。DataFrame 到 CSV 文件， 方便外部程序读取
atmos_df.to_csv(os.path.join(data_dir, 'kefa_atmos_data.csv'))
wave_df.to_csv(os.path.join(data_dir, 'kefa_wave_data.csv'))

# =============================================================================
# 第四部分: 新增 - 数据可视化
# =============================================================================

def plot_data(atmos_df, wave_df, lat, lon):
    """
    将提取出的大气和海浪数据进行可视化。
    """
    print("\n===== 正在生成数据图表... =====")

    # 检查是否有数据可供绘图
    if atmos_df.empty and wave_df.empty:
        print("--- [警告] 没有可供绘图的数据。 ---")
        return

    # 创建一个3行1列的图表布局，共享X轴
    fig, axes = plt.subplots(3, 1, figsize=(15, 12), sharex=True)
    fig.suptitle(f'Time Series Analysis for Location (Lat={lat:.2f}, Lon={lon:.2f})', fontsize=16)

    # --- 子图1: 大气条件 (SST & MSL) ---
    ax1 = axes[0]
    ax1_twin = ax1.twinx() # 创建双Y轴

    # 绘制海面温度
    color1 = 'tab:red'
    ax1.set_ylabel('Sea Surface Temperature (K)', color=color1, fontsize=12)
    ax1.plot(atmos_df.index, atmos_df['sst'], color=color1, label='SST')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.grid(True, linestyle='--', alpha=0.6)

    # 绘制海平面气压
    color2 = 'tab:blue'
    ax1_twin.set_ylabel('Mean Sea Level Pressure (Pa)', color=color2, fontsize=12)
    ax1_twin.plot(atmos_df.index, atmos_df['msl'], color=color2, label='MSL')
    ax1_twin.tick_params(axis='y', labelcolor=color2)
    
    ax1.set_title('Atmospheric Conditions', fontsize=14)
    # 合并两个Y轴的图例
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1_twin.legend(lines + lines2, labels + labels2, loc='upper right')

    # --- 子图2: 风速与风向 ---
    ax2 = axes[1]
    ax2_twin = ax2.twinx()

    # 绘制风速
    color1 = 'tab:green'
    ax2.set_ylabel('Wind Speed (m/s)', color=color1, fontsize=12)
    ax2.plot(atmos_df.index, atmos_df['wind_speed'], color=color1, label='Wind Speed')
    ax2.tick_params(axis='y', labelcolor=color1)
    ax2.grid(True, linestyle='--', alpha=0.6)

    # 绘制风向 (使用散点图更合适)
    color2 = 'tab:purple'
    ax2_twin.set_ylabel('Wind Direction (degrees)', color=color2, fontsize=12)
    ax2_twin.scatter(atmos_df.index, atmos_df['wind_direction'], color=color2, label='Wind Direction', s=5, alpha=0.7)
    ax2_twin.tick_params(axis='y', labelcolor=color2)
    ax2_twin.set_ylim(0, 360) # 风向范围是0-360度

    ax2.set_title('Wind Analysis', fontsize=14)
    lines, labels = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2_twin.legend(lines + lines2, labels + labels2, loc='upper right')


    # --- 子图3: 海浪参数 ---
    ax3 = axes[2]
    ax3_twin = ax3.twinx()

    # 绘制有效波高
    color1 = 'tab:cyan'
    ax3.set_ylabel('Significant Wave Height (m)', color=color1, fontsize=12)
    ax3.plot(wave_df.index, wave_df['swh'], color=color1, label='SWH')
    ax3.tick_params(axis='y', labelcolor=color1)
    ax3.grid(True, linestyle='--', alpha=0.6)

    # 绘制平均波周期
    color2 = 'tab:orange'
    ax3_twin.set_ylabel('Mean Wave Period (s)', color=color2, fontsize=12)
    ax3_twin.plot(wave_df.index, wave_df['mwp'], color=color2, label='MWP')
    ax3_twin.tick_params(axis='y', labelcolor=color2)

    ax3.set_title('Wave Properties', fontsize=14)
    lines, labels = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3_twin.legend(lines + lines2, labels + labels2, loc='upper right')
    
    # --- 格式化X轴 (日期) ---
    ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate() # 自动旋转日期标签防止重叠
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97]) # 调整布局以适应主标题
    print("--- 图表生成完毕，即将显示... ---")
    plt.show()
    fig.savefig(os.path.join(data_dir, 'kefalonia_data_analysis_plot.png'))

# 仅在DataFrame不为空时才调用绘图函数
if not atmos_df.empty or not wave_df.empty:
    plot_data(atmos_df, wave_df, target_lat, target_lon)
else:
    print("\n[信息] 由于未能提取到任何数据，跳过绘图步骤。")