# -----------------------------------------------------------------------------
# 脚本名称: CDS API 整合版下载与验证脚本
# 脚本功能: 
# 1. 严格按照用户从CDS网站生成的请求，拆分为两次独立的下载任务。
# 2. 任务一：下载大气数据。
# 3. 任务二：下载海洋波浪数据。
# 4. 下载完成后，自动分析并可视化大气数据中的风速，以作验证。
# -----------------------------------------------------------------------------

# --- 导入所需库 ---
import cdsapi
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os

# --- 全局设置 ---
# 定义输出文件名
atmospheric_filename = "kefalonia_atmospheric_data.nc"
wave_filename = "kefalonia_wave_data.nc"
# 建立API客户端
client = cdsapi.Client()

# --- 准备工作：删除可能存在的旧文件 ---
for f in [atmospheric_filename, wave_filename]:
    if os.path.exists(f):
        os.remove(f)
        print(f"已删除旧的数据文件: {f}")


# =============================================================================
# ===== 开始任务 1/2: 下载大气数据... =====
# =============================================================================
print("\n===== 任务 1/2: 开始请求大气数据... =====")
try:
    # 这是您请求中的大气部分
    client.retrieve(
        "reanalysis-era5-single-levels-monthly-means", # 大气数据集
        {
            "product_type": "monthly_averaged_reanalysis_by_hour_of_day",
            "variable": [
                "10m_u_component_of_wind", "10m_v_component_of_wind", "2m_dewpoint_temperature",
                "2m_temperature", "mean_sea_level_pressure", "surface_pressure",
                "total_precipitation",
            ],
            "year": "2024",
            "month": ["04", "05", "06"],
            "time": [
                "00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
                "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
                "12:00", "13:00", "14:00", "15:00", "16:00", "17:00",
                "18:00", "19:00", "20:00", "21:00", "22:00", "23:00",
            ],
            "area": [38.6, 20.4, 38.1, 21],
            "format": "netcdf",  # 使用NetCDF格式，更易于处理
        },
    ).download(atmospheric_filename) # 下载并指定文件名
    print(f"--- 任务1成功！大气数据已保存为: {atmospheric_filename} ---")
except Exception as e:
    print(f"!!!!!! 任务1失败: 下载大气数据时出错 !!!!!!")
    print(f"错误详情: {e}")

# =============================================================================
# ===== 开始任务 2/2: 下载海洋波浪数据... =====
# =============================================================================
print("\n===== 任务 2/2: 开始请求海洋波浪数据... =====")
try:
    # 这是您请求中的海洋波浪部分
    client.retrieve(
        "reanalysis-era5-wave-monthly-means",  # 注意：这里必须使用专门的海洋波浪数据集
        {
            "product_type": "monthly_averaged_reanalysis",
            "variable": [
                "mean_wave_direction", "mean_wave_period", 
                "significant_height_of_combined_wind_waves_and_swell",
                "sea_surface_temperature"
            ],
            "year": "2024",
            "month": ["04", "05", "06"],
            # 注意: 月平均海洋数据通常不区分小时
            "time": "00:00",
            "area": [38.6, 20.4, 38.1, 21],
            "format": "netcdf",
        },
    ).download(wave_filename) # 下载并指定文件名
    print(f"--- 任务2成功！海洋波浪数据已保存为: {wave_filename} ---")
except Exception as e:
    print(f"!!!!!! 任务2失败: 下载海洋波浪数据时出错 !!!!!!")
    print(f"错误详情: {e}")

# =============================================================================
# ===== 下载任务完成，开始验证和可视化... =====
# =============================================================================
if os.path.exists(atmospheric_filename):
    print(f"\n===== 开始处理和可视化大气数据: {atmospheric_filename} =====")
    try:
        ds = xr.open_dataset(atmospheric_filename)
        wind_speed = np.sqrt(ds.u10**2 + ds.v10**2)
        data_to_plot = wind_speed.isel(time=0)

        fig = plt.figure(figsize=(12, 10))
        ax = plt.axes(projection=ccrs.PlateCarree())
        data_to_plot.plot(ax=ax, transform=ccrs.PlateCarree(), robust=True, cmap='viridis', 
                          cbar_kwargs={'label': 'Wind Speed (m/s)'})
        ax.coastlines()
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
        plot_time = data_to_plot.time.dt.strftime('%Y-%m-%d %H:%M UTC').item()
        ax.set_title(f"10m Wind Speed over Kefalonia\n{plot_time}", fontsize=16)
        
        plt.show()
    except Exception as e:
        print(f"\n处理或绘图过程中发生错误: {e}")
else:
    print("\n未找到大气数据文件，无法进行可视化。")