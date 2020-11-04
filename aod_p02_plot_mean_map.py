#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-10-29 15:23
# @Author  : NingAnMe <ninganme@qq.com>
"""
绘制AOD数据月平均和季度平均的空间分布图
"""
import os
from lib.hdf5 import get_hdf5_data
from lib.plot_map import plot_map_project
from config import *


# ############################################## mean map ############################################

def plot_map_mean_5km(frequency='monthly'):
    aod_mean_file = ''
    aod = get_hdf5_data(aod_mean_file, 'aod_mean', 1, 0, [0, 1.5])
    lons = get_hdf5_data(aod_mean_file, 'lon', 1, 0, [-180, 180])
    lats = get_hdf5_data(aod_mean_file, 'lat', 1, 0, [-90, 90])

    box = [LATITUDE_RANGE[0], LATITUDE_RANGE[1], LONGITUDE_RANGE[1], LONGITUDE_RANGE[0]],  # nlat, slat, wlon, elon:北（小），南（大），东（大），西（小）
    title = '',
    vmin = 0,
    vmax = 1.5,
    markersize = 5
    filename_out = os.path.basename(aod_mean_file) + '.png'
    dir_out = os.path.join(AOD_MAP_DIR, 'AOD_MEAN_FY3D_5KM', frequency.upper())
    file_out = os.path.join(dir_out, filename_out)
    plot_map_project(lats, lons, aod, file_out,
                     box=box, title=title, vmin=vmin, vmax=vmax, markersize=markersize)


if __name__ == '__main__':
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

