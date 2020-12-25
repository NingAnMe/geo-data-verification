#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-10-29 15:23
# @Author  : NingAnMe <ninganme@qq.com>
"""
绘制AOD数据月平均和季度平均的空间分布图
"""

import os
from datetime import datetime
import argparse
from collections import defaultdict

import numpy as np

from lib.aod import AodCombine
from lib.proj_aod import proj_china
from lib.province_mask import get_province_mask

from config import AOD_COMBINE_DIR, AOD_PICTURE_DIR
from aod_h01_combine import get_day_str, get_month, get_season, get_year, get_month_str, get_season_str, get_year_str
from config import get_area_range, get_areas
from aod_p02_plot_map_origin import plot_map_picture

import warnings

warnings.filterwarnings('ignore')


def plot_map(datetime_start, datetime_end, data_dir_x=None, data_dir_y=None, out_dir=None,
             data_type_x=None, data_type_y=None, date_type=None):
    print("plot_map")
    print("datetime_start === {}".format(datetime_start))
    print("datetime_end === {}".format(datetime_end))
    print("data_dir_x === {}".format(data_dir_x))
    print("data_dir_y === {}".format(data_dir_y))
    print("out_dir === {}".format(out_dir))
    print("data_type_x === {}".format(data_type_x))
    print("data_type_y === {}".format(data_type_y))
    print("date_type === {}".format(date_type))

    file_dict_x = defaultdict(list)
    file_dict_y = defaultdict(list)
    for root, dirs, files in os.walk(data_dir_x):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            date_ = AodCombine(name).dt
            if not (datetime_start <= date_ <= datetime_end):
                continue
            if date_type == 'Monthly':
                date_str = get_month(date_)
            elif date_type == 'Seasonly':
                date_str = get_season(date_)
            elif date_type == 'Yearly':
                date_str = get_year(date_)
            else:
                raise ValueError(date_type)
            file_dict_x[date_str].append(file_dict_x)
    for root, dirs, files in os.walk(data_dir_y):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            date_ = AodCombine(name).dt
            if not (datetime_start <= date_ <= datetime_end):
                continue
            if date_type == 'Monthly':
                date_str = get_month(date_)
            elif date_type == 'Seasonly':
                date_str = get_season(date_)
            elif date_type == 'Yearly':
                date_str = get_year(date_)
            else:
                raise ValueError(date_type)
            file_dict_y[date_str].append(file_dict_y)

    for date_str in file_dict_x.keys():
        in_file_x = file_dict_x[date_str]
        if date_str not in file_dict_y:
            continue
        in_file_y = file_dict_y[date_str]
        print("<<< : {}".format(in_file_x))
        print("<<< : {}".format(in_file_y))
        for area_type in ["China", "YRD"]:
            loader_x = AodCombine(in_file_x)
            loader_y = AodCombine(in_file_y)
            title = get_title(data_type_x, data_type_y, date_type, loader_x.dt, area_type)
            filename = "_".join(title.split()) + '.png'
            out_file = os.path.join(out_dir, filename)

            datax = loader_x.get_aod()
            lonsx, latsx = loader_x.get_lon_lat()
            datax, lonsx, latsx = proj_china(datax, lonsx, latsx)

            datay = loader_y.get_aod()
            lonsy, latsy = loader_y.get_lon_lat()
            datay, lonsy, latsy = proj_china(datay, lonsy, latsy)

            data = datax - datay  # 真实值 - 参考值
            lons = lonsx
            lats = latsx

            vmin = 0
            vmax = 1.5

            ticks = np.arange(0, 1.51, 0.3)

            if area_type == 'China':
                nanhai = True
            else:
                nanhai = False

            mksize = 5

            areas = get_areas(area_type)
            mask = get_province_mask(areas)

            valid = np.logical_and.reduce((data > vmin, data < vmax, mask))

            data_mask = data[valid]
            lons_mask = lons[valid]
            lats_mask = lats[valid]

            longitude_range, latitude_range = get_area_range(area_type)
            box = [latitude_range[1], latitude_range[0], longitude_range[0], longitude_range[1]]

            plot_map_picture(data_mask, lons_mask, lats_mask, title=title, vmin=vmin, vmax=vmax,
                             areas=areas, box=box, ticks=ticks, file_out=out_file,
                             mksize=mksize, nanhai=nanhai)


get_dt_str = {
    'Daily': get_day_str,
    'Monthly': get_month_str,
    'Seasonly': get_season_str,
    'Yearly': get_year_str,
}


def get_title(data_type_x, data_type_y, date_type, dt, area_type):
    satellitex, sensorx = data_type_x.split('_')[:2]
    satellitey, sensory = data_type_y.split('_')[:2]
    date_str = get_dt_str[date_type](dt)

    title = "{} {}_{} {}_{} AOD(550nm) over {}".format(date_str, satellitex, sensorx, satellitey, sensory, area_type)
    return title


def main(data_type_x=None, data_type_y=None, date_start=None, date_end=None, date_type=None):
    datetime_start = datetime.strptime(date_start, "%Y%m%d")
    datetime_end = datetime.strptime(date_end, "%Y%m%d")

    combine_dir_x = os.path.join(AOD_COMBINE_DIR, 'AOD_COMBINE_{}'.format(data_type_x), date_type)
    combine_dir_y = os.path.join(AOD_COMBINE_DIR, 'AOD_COMBINE_{}'.format(data_type_y), date_type)
    picture_dir = os.path.join(AOD_PICTURE_DIR, 'AOD_MAP_{}_{}'.format(data_type_x, data_type_y), date_type)
    dataTypes = {'FY3D_MERSI_1KM', 'FY3D_MERSI_5KM', 'AQUA_MODIS_3KM', 'AQUA_MODIS_10KM'}

    if data_type_x in dataTypes and data_type_y in dataTypes:
        if date_type in {'Daily', 'Monthly', 'Seasonly', 'Yearly'}:
            combine_dir_x = combine_dir_x
            combine_dir_y = combine_dir_y
            out_dir = picture_dir
            plot_map(datetime_start, datetime_end, data_dir_x=combine_dir_x, data_dir_y=combine_dir_y, out_dir=out_dir,
                     data_type_x=data_type_x, data_type_y=data_type_y, date_type=date_type)
        else:
            parser.print_help()
            raise ValueError(date_type)
    else:
        parser.print_help()
        raise ValueError(data_type_x, data_type_y)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Schedule')
    parser.add_argument('--dataTypeX', help='数据类型：FY3D_MERSI_1KM、FY3D_MERSI_5KM、AQUA_MODIS_3KM、AQUA_MODIS_10KM',
                        required=True)
    parser.add_argument('--dataTypeY', help='数据类型：FY3D_MERSI_1KM、FY3D_MERSI_5KM、AQUA_MODIS_3KM、AQUA_MODIS_10KM',
                        required=True)
    parser.add_argument('--dateStart', help='开始时间（8位时间）：YYYYMMDD(20190101)', required=True)
    parser.add_argument('--dateEnd', help='结束时间（8位时间）：YYYYMMDD(20190102)', required=True)
    parser.add_argument('--dateType', help='合成的类型日、月、季、年：Daily、Monthly、Seasonly、Yearly', required=True)
    args = parser.parse_args()

    dataTypeX = args.dataTypeX
    dataTypeY = args.dataTypeY
    dateStart = args.dateStart
    dateEnd = args.dateEnd
    dateType = args.dateType

    main(data_type_x=dataTypeX, data_type_y=dataTypeY, date_start=dateStart, date_end=dateEnd, date_type=dateType)
