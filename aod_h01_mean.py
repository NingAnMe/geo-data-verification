#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-10-29 10:17
# @Author  : NingAnMe <ninganme@qq.com>
"""
计算AOD数据月平均和季度平均
"""
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta
from collections import defaultdict
from scipy.interpolate import griddata

from lib.path import *
from lib.aod import *
from lib.hdf5 import write_hdf5_and_compress

from config import *


def get_season(ym):
    season = {
        '201901': '201812',
        '201902': '201812',
        '201903': '201903',
        '201904': '201903',
        '201905': '201903',
        '201906': '201906',
        '201907': '201906',
        '201908': '201906',
        '201909': '201909',
        '201910': '201909',
        '201911': '201909',
        '201912': '201912',
        '202001': '201912',
        '202002': '201912',
        '202003': '202003',
        '202004': '202003',
        '202005': '202003',
    }
    if ym in season:
        return season[ym]
    else:
        return


def get_sum_count(aod_sum, aod_count, aod):
    if aod_sum is None:
        aod_sum = np.zeros_like(aod, dtype=np.float)
    if aod_count is None:
        aod_count = np.zeros_like(aod, dtype=np.int)
    index_valid = aod > 0
    aod_sum[index_valid] += aod[index_valid]
    aod_count[index_valid] += 1
    return aod_sum, aod_count


def get_mean(aod_sum, aod_count):
    index_valid = aod_count > 0
    aod_mean = np.full_like(aod_sum, -999, dtype=np.float)
    aod_mean[index_valid] = aod_sum[index_valid] / aod_count[index_valid]
    return aod_mean


def mean_modis_100km(frequency='monthly'):
    file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_MODIS_100KM_DIR):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            date_j_str = name.split('.')[1][1:]
            date_str = datetime.strptime(date_j_str, "%Y%j").strftime("%Y%m")
            file_dict[date_str].append(os.path.join(root, name))

    for season, files in file_dict.items():
        aod_sum = None
        aod_count = None
        filename_out = 'AOD_MODIS_10KM_MEAN_{}_{}.HDF'.format(frequency, season)
        dir_out = os.path.join(AOD_MEAN_DIR, 'AOD_MEAN_MODIS_10KM', frequency.upper())
        file_out = os.path.join(dir_out, filename_out)
        # if os.path.isfile(file_out):
        #     print('already exist {}'.format(file_out))
        #     continue

        for file_ in files:
            print('<<< {}'.format(file_))
            loader = AodModis100km(file_)
            aod = loader.get_aod()  # 无效值赋值为-999
            if aod is None:
                print('aod 为 None: {}'.format(file_))
                continue
            aod = aod[1]
            aod_sum, aod_count = get_sum_count(aod_sum, aod_count, aod)
        if aod_sum is not None and aod_count is not None:
            aod_mean = get_mean(aod_sum, aod_count)
        else:
            continue

        make_sure_path_exists(dir_out)
        lon, lat = AodModis100km.get_lon_lat()
        lon = lon.reshape(-1, 1)
        lat = lat.reshape(-1, 1)
        aod_mean = aod_mean.reshape(-1, 1)
        lat_lon = np.concatenate((lat, lon), axis=1)
        lon_5km, lat_5km = AodModis100km.get_lon_lat_5km()
        aod_mean_5km = griddata(lat_lon, aod_mean, (lat_5km, lon_5km), method='cubic')
        data_write = {
            'aod_mean': aod_mean_5km[:, :, 0],
            'lon': lon_5km,
            'lat': lat_5km
        }

        write_hdf5_and_compress(data_write, file_out)


def mean_fy3d_5km(frequency='monthly'):
    file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_FY3D_5KM_DIR):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            date_str = name.split('_')[7]
            if frequency == 'monthly':
                file_dict[date_str[:6]].append(os.path.join(root, name))
            elif frequency == 'seasonly':
                season = get_season(date_str[:6])
                if season is not None:
                    file_dict[season].append(os.path.join(root, name))
            else:
                continue

    for season, files in file_dict.items():
        aod_sum = None
        aod_count = None
        filename_out = 'AOD_FY3D_5KM_MEAN_{}_{}.HDF'.format(frequency, season)
        dir_out = os.path.join(AOD_MEAN_DIR, 'AOD_MEAN_FY3D_5KM', frequency.upper())
        file_out = os.path.join(dir_out, filename_out)
        if os.path.isfile(file_out):
            print('already exist {}'.format(file_out))
            continue

        for file_ in files:
            print('<<< {}'.format(file_))
            loader = AodFy3d5km(file_)
            aod = loader.get_aod()  # 无效值赋值为-999
            if aod is None:
                print('aod 为 None: {}'.format(file_))
                continue
            aod = aod[1]
            aod_sum, aod_count = get_sum_count(aod_sum, aod_count, aod)
        if aod_sum is not None and aod_count is not None:
            aod_mean = get_mean(aod_sum, aod_count)
        else:
            continue

        make_sure_path_exists(dir_out)
        lon, lat = AodFy3d5km.get_lon_lat()
        data_write = {
            'aod_mean': aod_mean,
            'lon': lon,
            'lat': lat
        }
        write_hdf5_and_compress(data_write, file_out)


if __name__ == '__main__':
    # mean_fy3d_5km('monthly')
    mean_modis_100km('monthly')

