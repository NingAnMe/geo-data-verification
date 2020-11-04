#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-09-04 9:40
# @Author  : NingAnMe <ninganme@qq.com>
"""
数据解析
数据匹配
回归图
偏差图
分布图
时间序列图
"""
from datetime import datetime
from collections import defaultdict
import os
import numpy as np
import pandas as pd

from lib.verification import Verification
from lib.path import make_sure_path_exists
from lib.aod import AodFy3d5km, AodModis
from config import *

LONGITUDE_RANGE = LONGITUDE_RANGE_China
LATITUDE_RANGE = LATITUDE_RANGE_China


def get_lon_lat_range_index(longitude_range, latitude_range, lons, lats):
    return np.logical_and.reduce((lons >= longitude_range[0],
                                 lons <= longitude_range[1],
                                 lats >= latitude_range[0],
                                 lats <= latitude_range[1]))


def match_fy3d_5km_modis_10km(aod_fy3d_5km_file, aod_modis_10km_file, aod_fy3d_5km_modis_10km_file,
                              longitude_range=None, latitude_range=None):
    print('<<<<< FY3D  : {}'.format(aod_fy3d_5km_file))
    print('<<<<< MODIS : {}'.format(aod_modis_10km_file))
    # 获取数据1
    try:
        fy3d_5km = AodFy3d5km(aod_fy3d_5km_file)
        data1 = fy3d_5km.get_aod()[1]
        lons1, lats1 = fy3d_5km.get_lon_lat()
    except Exception as why:
        print(why)
        return
    # 经纬度范围过滤
    if longitude_range is not None and latitude_range is not None:
        range_index1 = get_lon_lat_range_index(longitude_range, latitude_range, lons1, lats1)
        if range_index1.sum() <= 0:
            print('data在经纬度范围内的数据不足')
            return
        else:
            data1 = data1[range_index1]
            lons1 = lons1[range_index1]
            lats1 = lats1[range_index1]

    # 获取数据1
    modis_10km = AodModis(aod_modis_10km_file)
    data2 = modis_10km.get_aod()
    lons2, lats2 = modis_10km.get_lon_lat()
    # 经纬度范围过滤
    if longitude_range is not None and latitude_range is not None:
        range_index2 = get_lon_lat_range_index(longitude_range, latitude_range, lons2, lats2)
        if range_index2.sum() <= 0:
            print('data2在经纬度范围内的数据不足')
            return
        else:
            data2 = data2[range_index2]
            lons2 = lons2[range_index2]
            lats2 = lats2[range_index2]

    print(data1.shape)
    print(lons1.shape)
    print(lats1.shape)
    print(data2.shape)
    print(lons2.shape)
    print(lats2.shape)

    data_kdtree, lons_kdtree, lats_kdtree = data2, lons2, lats2
    data_query, lons_query, lats_query = data1, lons1, lats1

    # 数据匹配
    verif = Verification(lons_kdtree, lats_kdtree, lons_query, lats_query)
    if not verif.get_kdtree():
        return

    verif.get_dist_and_index_kdtree()

    # 获取符合距离阈值的点
    pre_dist = 0.03
    index_dist = verif.get_index_dist(pre_dist=pre_dist)

    # 匹配数据
    if index_dist.sum() > 0:
        print('匹配的数据量: {}'.format(index_dist.sum()))
        result = {
            'lons_y': verif.get_kdtree_data(lons_kdtree)[index_dist],
            'lats_y': verif.get_kdtree_data(lats_kdtree)[index_dist],
            'aod_y': verif.get_kdtree_data(data_kdtree)[index_dist],
            'lons_x': verif.get_query_data(lons_query)[index_dist],
            'lats_x': verif.get_query_data(lats_query)[index_dist],
            'aod_x': verif.get_query_data(data_query)[index_dist],
            'dist': verif.get_query_data(verif.dist)[index_dist],
        }
        pd.DataFrame(result).to_csv(aod_fy3d_5km_modis_10km_file, index=False)
        print('>>>>> {}'.format(aod_fy3d_5km_modis_10km_file))
        return result
    else:
        print('匹配的数据量 < 0')
        return


def get_out_file(fy3d_file, modis_file, out_dir):
    fy3d_file_name = os.path.basename(fy3d_file)
    modis_file_name = os.path.basename(modis_file)
    fy3d_date = fy3d_file_name.split('_')[7]
    modis_date = ''.join(modis_file_name.split('.')[1:3])
    out_file = os.path.join(out_dir, 'FY3D_{}_MODIS_{}.csv'.format(fy3d_date, modis_date))
    return out_file


def match_file():

    fy3d_file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_FY3D_5KM_DIR):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            date_str = name.split('_')[7]
            fy3d_file_dict[date_str].append(os.path.join(root, name))

    modis_file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_MODIS_10KM_DIR):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            date_j_str = name.split('.')[1][1:]
            date_str = datetime.strptime(date_j_str, "%Y%j").strftime("%Y%m%d")
            modis_file_dict[date_str].append(os.path.join(root, name))

    for fy3d_date_str in fy3d_file_dict.keys():
        if fy3d_date_str in modis_file_dict:
            fy3d_files = fy3d_file_dict[fy3d_date_str]
            modis_files = modis_file_dict[fy3d_date_str]
            for fy3d_file in fy3d_files:
                for modis_file in modis_files:
                    match_dir = os.path.join(AOD_FY3D_5KM_MODIS_10KM_DIR, 'MATCH')
                    make_sure_path_exists(match_dir)
                    out_file = get_out_file(fy3d_file, modis_file, match_dir)
                    if os.path.isfile(out_file):
                        print('already exist {}'.format(out_file))
                        continue
                    match_fy3d_5km_modis_10km(fy3d_file, modis_file, out_file,
                                              longitude_range=LONGITUDE_RANGE, latitude_range=LATITUDE_RANGE)


def t():
    aod_fy3d_5km_file = r"C:\F\shanghai\AOD\1FY3D-MERSI-L2-AOD-DAILY-5000M\2019\FY3D_MERSI_GBAL_L2_AOD_MLT_GLL_20190102_POAD_5000M_MS.HDF"
    aod_modis_10km_file = r"C:\F\shanghai\AOD\3MYD04_L2_MODIS_AOD_10KM\2019\MYD04_L2.A2019001.0420.061.2019001165304.hdf"
    aod_fy3d_5km_modis_10km_file = 'test/test.csv'
    longitude_range = [70, 140]
    latitude_range = [15, 56]
    match_fy3d_5km_modis_10km(aod_fy3d_5km_file, aod_modis_10km_file, aod_fy3d_5km_modis_10km_file,
                              longitude_range=longitude_range, latitude_range=latitude_range)


if __name__ == '__main__':
    match_file()
