#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-10-29 15:23
# @Author  : NingAnMe <ninganme@qq.com>
"""
绘制AOD数据月平均和季度平均的空间分布图
"""
import os

from lib.hdf5 import get_hdf5_data
from config import *
from DV.dv_map import dv_map
from pylab import *


# ############################################## mean map ############################################

def plot_map_mean_5km(frequency='monthly'):
    dir_out = os.path.join(AOD_MEAN_DIR, SATELLITE, frequency.upper())
    print("<<< {}".format(dir_out))
    filenames = os.listdir(dir_out)
    files = [os.path.join(dir_out, x) for x in filenames]
    for aod_mean_file in files:
        print("<<< {}".format(aod_mean_file))

        month = os.path.basename(aod_mean_file).split('_')[5]

        box = [LATITUDE_RANGE[1], LATITUDE_RANGE[0], LONGITUDE_RANGE[0], LONGITUDE_RANGE[1]]  # nlat, slat, wlon, elon:北（小），南（大），东（大），西（小）
        title = '{} {}'.format(month, AREA)
        vmin = 0
        vmax = 1.5
        markersize = 5
        filename_out = '{}_'.format(AREA) + os.path.basename(aod_mean_file) + '.png'
        dir_out = os.path.join(AOD_MAP_DIR, SATELLITE, frequency.upper())
        file_out = os.path.join(dir_out, filename_out)
        if os.path.isfile(file_out):
            print('already exist {}'.format(file_out))
            continue

        aod = get_hdf5_data(aod_mean_file, 'aod_mean', 1, 0, [0, 1.5])
        lons = get_hdf5_data(aod_mean_file, 'lon', 1, 0, [-180, 180])
        lats = get_hdf5_data(aod_mean_file, 'lat', 1, 0, [-90, 90])

        latitude = lats
        longitude = lons
        value = aod
        out_file = file_out

        # 开始画图-----------------------

        fig = plt.figure(figsize=(9., 8.))  # 图像大小

        p = dv_map(fig=fig)

        subplots_adjust(left=0.07, right=0.98, top=0.90, bottom=0.15)
        p.show_colorbar = False
        p.show_countries = False
        p.show_coastlines = False
        # p.show_line_of_latlon = False
        # p.show_china = True
        p.show_china_province = True
        p.show_inside_china = True
        p.show_inside_china_mini = False

        p.delat = 2  # 纬度刻度线分辨率
        p.delon = 2  # 经度刻度线分辨率
        p.color_coast = "#3a3a3a"  # 海岸线颜色
        p.color_contry = "#3a3a3a"  # 国家颜色
        p.fontsize_tick = 15

        # set color map
        p.valmin = vmin
        p.valmax = vmax
        p.colormap = plt.get_cmap('jet')  # mpl.cm.rainbow, summer
        # p.colorbar_extend = "max"

        if AREA == 'ChangSanJiao':
            print('设置地区')
            p.set_area(["江苏省", "安徽省", "浙江省", "上海市"])

        # plot
        p.easyplot(latitude, longitude, value, vmin=vmin, vmax=vmax, box=box, markersize=markersize, ptype="pcolormesh")

        # 色标 ---------------------------
        cb_loc = [0.12, 0.07, 0.76, 0.03]
        # unit = r"$\mathregular{(10^{15}\/\/molec/cm^2)}$"
        fontsize = 16
        # p.add_custom_colorbar(cb_loc, p.valmin, p.valmax,
        #                       fmt="%d",
        #                       unit="(1E16 molec/cm^2)",
        #                       fontsize=fontsize)
        c_ax = fig.add_axes(cb_loc)
        cbar = fig.colorbar(p.cs, cax=c_ax, orientation='horizontal')
        for l in c_ax.xaxis.get_ticklabels():
            l.set_fontproperties(p.font_mid)
            l.set_fontsize(fontsize)
            l.set_color(p.color_ticker)
        # cbar写单位
        # cbar.ax.set_title(unit, x=1.0382, y=0, color=p.color_ticker,
        #                   ha='left', va='center',
        #                   fontproperties=p.font_mid, fontsize=fontsize)

        # 标题 ---------------------------
        p.w_title = p.suptitle(title, fontsize=22, y=0.97)

        # save
        p.savefig(out_file, dpi=300)
        print(">>> {}".format(out_file))
        p.clean()


if __name__ == '__main__':
    for SATELLITE in ['AOD_MEAN_FY3D_5KM', 'AOD_MEAN_MODIS_10KM']:
        AREAs = ['ChangSanJiao', 'ZhuSanJiao', 'FenWei', 'JingJinJi']

        for AREA in AREAs:

            if AREA == 'China':
                LONGITUDE_RANGE = LONGITUDE_RANGE_China
                LATITUDE_RANGE = LATITUDE_RANGE_China
            elif AREA == 'ChangSanJiao':
                LONGITUDE_RANGE = LONGITUDE_RANGE_ChangSanJiao
                LATITUDE_RANGE = LATITUDE_RANGE_ChangSanJiao
            elif AREA == 'ZhuSanJiao':
                LONGITUDE_RANGE = LONGITUDE_RANGE_ZhuSanJiao
                LATITUDE_RANGE = LATITUDE_RANGE_ZhuSanJiao
            elif AREA == 'FenWei':
                LONGITUDE_RANGE = LONGITUDE_RANGE_FenWei
                LATITUDE_RANGE = LATITUDE_RANGE_FenWei
            elif AREA == 'JingJinJi':
                LONGITUDE_RANGE = LONGITUDE_RANGE_JingJinJi
                LATITUDE_RANGE = LATITUDE_RANGE_JingJinJi
            else:
                LONGITUDE_RANGE = None
                LATITUDE_RANGE = None

            plot_map_mean_5km(frequency='monthly')
