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


def get_season(ym):
    season = {
        '201812': '2018 DJF',
        '201903': '2019 MAM',
        '201906': '2019 JJA',
        '201909': '2019 SON',
        '201912': '2019 DJF',
        '202003': '2020 MAM',
    }
    return season[ym]


# ############################################## mean map ############################################

def plot_map_mean(frequency='monthly'):
    dir_out = os.path.join(AOD_MEAN_DIR, SATELLITE, frequency.upper())
    print("dir_out <<< {}".format(dir_out))
    filenames = os.listdir(dir_out)
    files = [os.path.join(dir_out, x) for x in filenames]

    mode = None
    if 'MEAN' in SATELLITE:
        mode = 'MEAN'
    elif 'DIFF' in SATELLITE:
        mode = 'DIFF'
    else:
        assert KeyError('不支持的绘图类型')

    for aod_mean_file in files:
        print("<<< {}".format(aod_mean_file))

        if 'FY3D' in SATELLITE:
            sate_ = "FY3D MERSI"
        elif 'MODIS' in SATELLITE:
            sate_ = "AQUA MODIS"
        else:
            sate_ = ''
        month = os.path.basename(aod_mean_file).split('_')[5][:6]

        if frequency == 'monthly':
            d_ = month
        elif frequency == 'seasonly':
            d_ = get_season(month)
        elif frequency == 'all':
            d_ = '201901-202005'
        else:
            d_ = ''

        box_ = [LATITUDE_RANGE[1], LATITUDE_RANGE[0], LONGITUDE_RANGE[0], LONGITUDE_RANGE[1]]  # nlat, slat, wlon, elon:北（小），南（大），东（大），西（小）
        if 'MERSI' in SATELLITE and 'MODIS' in SATELLITE:
            title_ = '{} Diff.(MERSI-MODIS) AOD over {}'.format(d_, AREA)
            vmin = -0.5
            vmax = 0.5
            ticks = np.arange(-0.5, 0.51, 0.1)
        else:
            title_ = '{} {} AOD (550nm) over {}'.format(d_, sate_, AREA)
            vmin = 0
            vmax = 1.5
            ticks = np.arange(0, 1.6, 0.3)

        markersize = 5
        filename_out = '{}_'.format(AREA) + os.path.basename(aod_mean_file) + '.png'
        dir_out = os.path.join(AOD_PICTURE_DIR, 'MAP', SATELLITE, frequency.upper())
        file_out = os.path.join(dir_out, filename_out)

        # 是否重处理
        # if os.path.isfile(file_out):
        #     print('already exist {}'.format(file_out))
        #     continue

        aod = get_hdf5_data(aod_mean_file, 'aod_mean', 1, 0, [vmin, vmax], np.nan)
        lons = get_hdf5_data(aod_mean_file, 'lon', 1, 0, [-180, 180])
        lats = get_hdf5_data(aod_mean_file, 'lat', 1, 0, [-90, 90])
        print(np.nanmin(aod), np.nanmax(aod), np.nanmean(aod))

        latitude = lats
        longitude = lons
        value = aod
        out_file = file_out
        valid = np.logical_and(aod > vmin, aod < vmax)
        value[~valid] = np.nan

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
        if mode == 'MEAN':
            p.colormap = plt.get_cmap('jet')  # mpl.cm.rainbow, summer, jet, bwr
        elif mode == 'DIFF':
            p.colormap = plt.get_cmap('bwr')
        # p.colorbar_extend = "max"

        # plot
        p.easyplot(latitude, longitude, value, vmin=vmin, vmax=vmax, box=box_, markersize=markersize, ptype="pcolormesh")

        print('设置地区 ：{}'.format(AREA))
        if AREA == 'YRD':
            citys = ["江苏省", "安徽省", "浙江省", "上海市"]
        elif AREA == 'BTH':
            citys = ["北京市", "天津市", "河北省"]
        elif AREA == 'FWP':
            citys = ["陕西省", "山西省", "河南省"]
        elif AREA == 'FWP':
            citys = ["广东省"]
        else:
            citys = []
        for city in citys:
            p.city_boundary(city, linewidth=1.2, shape_name='中国省级行政区')

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
        p.w_title = p.suptitle(title_, fontsize=14, y=0.97)

        # save
        p.savefig(out_file, dpi=300)
        print(">>> {}".format(out_file))
        p.clean()


if __name__ == '__main__':
    # SATELLITEs = ['AOD_MEAN_FY3D_5KM', 'AOD_MEAN_MODIS_10KM', 'AOD_DIFF_MERSI_MODIS']
    SATELLITEs = ['AOD_DIFF_MERSI_MODIS']
    for SATELLITE in SATELLITEs:
        AREAs = ['YRD', 'PRD', 'FWP', 'BTH']

        for AREA in AREAs:

            if AREA == 'China':
                LONGITUDE_RANGE = LONGITUDE_RANGE_China
                LATITUDE_RANGE = LATITUDE_RANGE_China
            elif AREA == 'YRD':
                LONGITUDE_RANGE = LONGITUDE_RANGE_ChangSanJiao
                LATITUDE_RANGE = LATITUDE_RANGE_ChangSanJiao
            elif AREA == 'PRD':
                LONGITUDE_RANGE = LONGITUDE_RANGE_ZhuSanJiao
                LATITUDE_RANGE = LATITUDE_RANGE_ZhuSanJiao
            elif AREA == 'FWP':
                LONGITUDE_RANGE = LONGITUDE_RANGE_FenWei
                LATITUDE_RANGE = LATITUDE_RANGE_FenWei
            elif AREA == 'BTH':
                LONGITUDE_RANGE = LONGITUDE_RANGE_JingJinJi
                LATITUDE_RANGE = LATITUDE_RANGE_JingJinJi
            else:
                LONGITUDE_RANGE = None
                LATITUDE_RANGE = None

            plot_map_mean(frequency='monthly')
            plot_map_mean(frequency='seasonly')
            plot_map_mean(frequency='all')
