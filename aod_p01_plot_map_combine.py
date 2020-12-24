#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-12-23 11:49
# @Author  : NingAnMe <ninganme@qq.com>
import os
from datetime import datetime
import argparse

import numpy as np
import matplotlib.pyplot as plt
from pylab import subplots_adjust

from DV.dv_map import dv_map

from lib.aod import AodFy3d5km, AodCombine

from config import AOD_COMBINE_DIR, AOD_FY3D_5KM_DIR, AOD_PICTURE_DIR
from aod_h01_combine import get_day_str, get_month_str, get_season_str, get_year_str
from config import get_area_range, get_areas


def plot_map_picture(data, lons, lats, title='', vmin=-np.inf, vmax=np.inf, areas=None, box=None, ticks=None,
                     file_out=None, ptype='pcolormesh', mksize=5, nanhai=False):
    if file_out is None:
        print('没有指定输出文件：file_out is None')
        return

    latitude = lats
    longitude = lons
    value = data
    out_file = file_out
    if vmin is not None and vmax is not None:
        valid = np.logical_and(value > vmin, value < vmax)
        value[~valid] = np.nan
    else:
        value[value == 0] = np.nan
        vmin = np.nanmin(value)
        vmax = np.nanmax(value)

    print(np.nanmin(data), np.nanmax(data), np.nanmean(data))

    # 开始画图-----------------------

    fig = plt.figure(figsize=(9, 8))  # 图像大小

    p = dv_map(fig=fig)

    subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.15)
    p.show_colorbar = False
    p.show_countries = False
    p.show_coastlines = False
    # p.show_line_of_latlon = False
    # p.show_china = True
    if nanhai:
        p.nanhai_loc = [0.83, 0.25, 0.15, 0.17]
        p.nanhai_minimap()
    p.show_china_province = True
    p.show_inside_china = True
    p.show_inside_china_mini = True

    if box:
        if abs(box[1] - box[0]) < 10:
            p.delat = 2
        else:
            p.delat = 5
        if abs(box[2] - box[3]) < 10:
            p.delon = 2
        elif abs(box[2] - box[3]) < 40:
            p.delon = 5
        else:
            p.delon = 10
    # else:
    #     p.delat = 2  # 纬度刻度线分辨率
    #     p.delon = 2  # 经度刻度线分辨率
    p.color_coast = "#3a3a3a"  # 海岸线颜色
    p.color_contry = "#3a3a3a"  # 国家颜色
    p.fontsize_tick = 15

    # set color map
    p.valmin = vmin
    p.valmax = vmax
    p.colormap = plt.get_cmap('jet')  # mpl.cm.rainbow, summer, jet, bwr
    # p.colorbar_extend = "max"

    # plot
    p.easyplot(latitude, longitude, value, vmin=vmin, vmax=vmax, box=box, markersize=mksize, ptype=ptype)

    if areas is not None and len(areas) > 0:
        print('设置地区 ：{}'.format(areas))
        # aeres = ["江苏省", "安徽省", "浙江省", "上海市"]
        for aere in areas:
            p.city_boundary(aere, linewidth=1.2, shape_name='中国省级行政区')
        p.set_area(areas)

    # 色标 ---------------------------
    cb_loc = [0.12, 0.07, 0.76, 0.03]
    # unit = r"$\mathregular{(10^{15}\/\/molec/cm^2)}$"
    fontsize = 16
    # p.add_custom_colorbar(cb_loc, p.valmin, p.valmax,
    #                       fmt="%d",
    #                       unit="(1E16 molec/cm^2)",
    #                       fontsize=fontsize)
    c_ax = fig.add_axes(cb_loc)
    # cbar = fig.colorbar(p.cs, cax=c_ax, ticks=np.arange(0, 1.6, 0.3), orientation='horizontal')
    fig.colorbar(p.cs, cax=c_ax, ticks=ticks, orientation='horizontal')
    for l in c_ax.xaxis.get_ticklabels():
        l.set_fontproperties(p.font_mid)
        l.set_fontsize(fontsize)
        l.set_color(p.color_ticker)
    # cbar写单位
    # cbar.ax.set_title(unit, x=1.0382, y=0, color=p.color_ticker,
    #                   ha='left', va='center',
    #                   fontproperties=p.font_mid, fontsize=fontsize)

    # 标题 ---------------------------
    p.w_title = p.suptitle(title, fontsize=14, y=0.97)

    # save
    p.savefig(out_file, dpi=300)
    print(">>> {}".format(out_file))
    p.clean()


def plot_map(datetime_start, datetime_end, data_dir=None, out_dir=None, data_loader=AodCombine,
             data_type=None, date_type=None):
    print("plot_map")
    print("datetime_start === {}".format(datetime_start))
    print("datetime_end === {}".format(datetime_end))
    print("data_dir === {}".format(data_dir))
    print("out_dir === {}".format(out_dir))
    print("data_loader === {}".format(data_loader))
    print("data_type === {}".format(data_type))
    print("date_type === {}".format(date_type))
    filelist = list()
    for root, dirs, files in os.walk(data_dir):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            date_ = data_loader(name).dt
            if not (datetime_start <= date_ <= datetime_end):
                continue
            filelist.append(os.path.join(root, name))
    filelist.sort()
    for in_file in filelist:
        print("<<< : {}".format(in_file))
        for area_type in ["China", "YRD"]:
            loader = data_loader(in_file)
            title = get_title(data_type, date_type, loader.dt, area_type)
            filename = "_".join(title.split()) + '.png'
            out_file = os.path.join(out_dir, filename)

            data = loader.get_aod()
            lons, lats = loader.get_lon_lat()

            vmin = 0
            vmax = 1.5

            ticks = np.arange(0, 1.51, 0.3)

            if area_type == 'China':
                nanhai = True
            else:
                nanhai = False

            mksize = 5

            areas = get_areas(area_type)

            longitude_range, latitude_range = get_area_range(area_type)
            box = [latitude_range[1], latitude_range[0], longitude_range[0], longitude_range[1]]

            plot_map_picture(data, lons, lats, title=title, vmin=vmin, vmax=vmax,
                             areas=areas, box=box, ticks=ticks, file_out=out_file,
                             mksize=mksize, nanhai=nanhai)


get_dt_str = {
    'Daily': get_day_str,
    'Monthly': get_month_str,
    'Seasonly': get_season_str,
    'Yearly': get_year_str,
}


def get_title(data_type, date_type, dt, area_type):
    satellite, sensor = data_type.split('_')[:2]
    date_str = get_dt_str[date_type](dt)

    title = "{} {} {} AOD(550nm) over {}".format(date_str, satellite, sensor, area_type)
    return title


def main(data_type=None, date_start=None, date_end=None, date_type=None):
    datetime_start = datetime.strptime(date_start, "%Y%m%d")
    datetime_end = datetime.strptime(date_end, "%Y%m%d")

    combine_dir = os.path.join(AOD_COMBINE_DIR, 'AOD_COMBINE_{}'.format(data_type))
    picture_dir = os.path.join(AOD_PICTURE_DIR, 'AOD_MAP_{}'.format(data_type))

    if data_type in {'FY3D_MERSI_1KM', 'AQUA_MODIS_3KM', 'AQUA_MODIS_10KM'}:

        if date_type in {'Daily', 'Monthly', 'Seasonly', 'Yearly'}:
            data_dir = os.path.join(combine_dir, date_type)
            out_dir = os.path.join(picture_dir, date_type)
            plot_map(datetime_start, datetime_end, data_dir=data_dir, out_dir=out_dir,
                     data_type=data_type, date_type=date_type)
        else:
            parser.print_help()
            raise ValueError(date_type)
    elif data_type == 'FY3D_MERSI_5KM':
        if date_type in {'Daily'}:
            data_dir = AOD_FY3D_5KM_DIR
            out_dir = os.path.join(picture_dir, date_type)
            plot_map(datetime_start, datetime_end, data_dir=data_dir, out_dir=out_dir, data_loader=AodFy3d5km,
                     data_type=data_type, date_type=date_type)
        elif date_type in {'Monthly', 'Seasonly', 'Yearly'}:
            data_dir = os.path.join(combine_dir, date_type)
            out_dir = os.path.join(picture_dir, date_type)
            plot_map(datetime_start, datetime_end, data_dir=data_dir, out_dir=out_dir,
                     data_type=data_type, date_type=date_type)
        else:
            parser.print_help()
            raise ValueError(date_type)
    else:
        parser.print_help()
        raise ValueError(data_type)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Schedule')
    parser.add_argument('--dataType', help='数据类型：FY3D_MERSI_1KM、FY3D_MERSI_5KM、AQUA_MODIS_3KM、AQUA_MODIS_5KM',
                        required=True)
    parser.add_argument('--dateStart', help='开始时间（8位时间）：YYYYMMDD(20190101)', required=True)
    parser.add_argument('--dateEnd', help='结束时间（8位时间）：YYYYMMDD(20190102)', required=True)
    parser.add_argument('--dateType', help='合成的类型日、月、季、年：Daily、Monthly、Seasonly、Yearly', required=True)
    args = parser.parse_args()

    dataType = args.dataType
    dateStart = args.dateStart
    dateEnd = args.dateEnd
    dateType = args.dateType

    main(data_type=dataType, date_start=dateStart, date_end=dateEnd, date_type=dateType)

    """
    绘制 FY3D_MERSI_1KM 日数据的分布图
    python3 aod_p01_plot_map.py --dataType FY3D_MERSI_1KM --dateStart 20190101 --dateEnd 20190131 --dateType Daily

    绘制 AQUA_MODIS_3KM 日数据的分布图
    python3 aod_p01_plot_map.py --dataType AQUA_MODIS_3KM --dateStart 20190101 --dateEnd 20190131 --dateType Daily

    绘制 AQUA_MODIS_10KM 日数据的分布图
    python3 aod_p01_plot_map.py --dataType AQUA_MODIS_10KM --dateStart 20190101 --dateEnd 20190131 --dateType Daily

    绘制 FY3D_MERSI_1KM 月数据的分布图
    python3 aod_p01_plot_map.py --dataType FY3D_MERSI_1KM --dateStart 20190101 --dateEnd 20190131 --dateType Monthly

    绘制 AQUA_MODIS_3KM 月数据的分布图
    python3 aod_p01_plot_map.py --dataType FY3D_MERSI_5KM --dateStart 20190101 --dateEnd 20190131 --dateType Monthly

    绘制 AQUA_MODIS_3KM 月数据的分布图
    python3 aod_p01_plot_map.py --dataType AQUA_MODIS_3KM --dateStart 20190101 --dateEnd 20190131 --dateType Monthly
    
    绘制 AQUA_MODIS_10KM 月数据的分布图
    python3 aod_p01_plot_map.py --dataType AQUA_MODIS_10KM --dateStart 20190101 --dateEnd 20190131 --dateType Monthly
    """