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
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta
from collections import defaultdict
import numpy as np
import pandas as pd

from lib.verification import Verification
from lib.path import make_sure_path_exists
from lib.aod import AodFy3d1km, AodFy4a4km
from config import LONGITUDE_RANGE_China, LATITUDE_RANGE_China
from config import AOD_FY3D_1KM_DIR, AOD_FY4A_4KM_DIR, GEO_FY3D_1KM_DIR, AOD_MATCH_DIR

import warnings
warnings.filterwarnings('ignore')

LONGITUDE_RANGE = LONGITUDE_RANGE_China
LATITUDE_RANGE = LATITUDE_RANGE_China


def get_lon_lat_range_index(longitude_range, latitude_range, lons, lats):
    return np.logical_and.reduce((lons >= longitude_range[0],
                                  lons <= longitude_range[1],
                                  lats >= latitude_range[0],
                                  lats <= latitude_range[1]))


def match_fy3d_1km_fy4a_4km(aod_fy3d_1km_file, geo_fy3d_1km_file, aod_fy4a_4km_file, aod_fy3d_1km_fy4a_4km_file,
                            longitude_range=None, latitude_range=None):
    print('<<<<< FY3D L1  : {}'.format(aod_fy3d_1km_file))
    print('<<<<< FY3D GEO : {}'.format(geo_fy3d_1km_file))
    print('<<<<< FY4A L1  : {}'.format(aod_fy4a_4km_file))
    # 获取数据1
    try:
        fy3d_1km = AodFy3d1km(aod_fy3d_1km_file, geo_file=geo_fy3d_1km_file)
        data1 = fy3d_1km.get_aod()
        lons1, lats1 = fy3d_1km.get_lon_lat()
    except Exception as why:
        print(why)
        return
    # 经纬度范围过滤
    if longitude_range is not None and latitude_range is not None:
        range_index1 = get_lon_lat_range_index(longitude_range, latitude_range, lons1, lats1)
        if range_index1.sum() <= 0:
            print('data在经纬度范围内的数据不足:{}'.format(range_index1.sum()))
            return
        else:
            data1 = data1[range_index1]
            lons1 = lons1[range_index1]
            lats1 = lats1[range_index1]

    # 有效值过滤
    range_index1 = np.logical_and(data1 > 0, data1 < 1.5)
    if range_index1.sum() <= 0:
        print('data1的有效数据不足')
        return
    else:
        data1 = data1[range_index1]
        lons1 = lons1[range_index1]
        lats1 = lats1[range_index1]

    # 获取数据2
    try:
        fy4a_4km = AodFy4a4km(aod_fy4a_4km_file)
        data2 = fy4a_4km.get_aod()
        lons2, lats2 = fy4a_4km.get_lon_lat()
    except Exception as why:
        print(why)
        return
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

    # 有效值过滤
    range_index2 = np.logical_and(data2 > 0, data2 < 1.5)
    if range_index2.sum() <= 0:
        print('data2的有效数据不足')
        return
    else:
        data2 = data2[range_index2]
        lons2 = lons2[range_index2]
        lats2 = lats2[range_index2]

    data_kdtree, lons_kdtree, lats_kdtree = data1, lons1, lats1
    data_query, lons_query, lats_query = data2, lons2, lats2

    # 数据匹配
    verif = Verification(lons_kdtree, lats_kdtree, lons_query, lats_query)
    if not verif.get_kdtree():
        return
    # kdtree_model_pickle = 'kdtree_model.pickle'
    # if os.path.isfile(kdtree_model_pickle):
    #     with open(kdtree_model_pickle, 'rb') as fp:
    #         verif.kdtree_model = pickle.load(fp)
    # else:
    #     if not verif.get_kdtree():
    #         return
    #     else:
    #         with open(kdtree_model_pickle, 'wb') as fp:
    #             pickle.dump(verif.kdtree_model, fp)

    verif.get_dist_and_index_kdtree()

    # 获取符合距离阈值的点
    pre_dist = 0.009
    index_dist = verif.get_index_dist(pre_dist=pre_dist)

    # 匹配数据
    if index_dist.sum() > 0:
        print('匹配的数据量: {}'.format(index_dist.sum()))
        result = {
            'lons_x': verif.get_kdtree_data(lons_kdtree)[index_dist],
            'lats_x': verif.get_kdtree_data(lats_kdtree)[index_dist],
            'aod_x': verif.get_kdtree_data(data_kdtree)[index_dist],
            'lons_y': verif.get_query_data(lons_query)[index_dist],
            'lats_y': verif.get_query_data(lats_query)[index_dist],
            'aod_y': verif.get_query_data(data_query)[index_dist],
            'dist': verif.get_query_data(verif.dist)[index_dist],
        }
        pd.DataFrame(result).to_csv(aod_fy3d_1km_fy4a_4km_file, index=False)
        print('>>>>> {}'.format(aod_fy3d_1km_fy4a_4km_file))
        return result
    else:
        print('匹配的数据量 < 0')
        return


def get_out_file(fy3d_file, fy4a_file, out_dir):
    fy3d_file_name = os.path.basename(fy3d_file)
    fy4a_file_name = os.path.basename(fy4a_file)
    ymdhm_fy3d = fy3d_file_name.split('_')[7:9]
    ymdhms_fy4a = fy4a_file_name.split('_')[9][0:12]
    out_file = os.path.join(out_dir, 'FY3D_{}_{}_FY4A_{}_{}.csv'.format(ymdhm_fy3d[0], ymdhm_fy3d[1],
                                                                        ymdhms_fy4a[0:8], ymdhms_fy4a[8:12]))
    return out_file


def match_file():
    # FY3D_MERSI_ORBT_L2_AOD_MLT_NUL_20190102_0500_1000M_MS.HDF
    fy3d_file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_FY3D_1KM_DIR):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue
            ymd, hm = name.split('_')[7:9]
            date_str = ymd
            fy3d_file_dict[date_str].append(os.path.join(root, name))

    # SATE_FY4A_AGRI_MULT_4000M_DISK_LDA_GLL_HYT_20190101000000_VPRJ5_L2.NC
    fy4a_file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_FY4A_4KM_DIR):
        for name in files:
            if name[-2:].lower() != 'nc':
                continue
            date_str = name.split('_')[9][0:8]
            fy4a_file_dict[date_str].append(os.path.join(root, name))

    for fy3d_date_str in fy3d_file_dict.keys():
        if fy3d_date_str in fy4a_file_dict:
            fy3d_files = fy3d_file_dict[fy3d_date_str]
            fy4a_files = fy4a_file_dict[fy3d_date_str]
            for fy3d_file in fy3d_files:
                for fy4a_file in fy4a_files:
                    ymdhm_fy3d = os.path.basename(fy3d_file).split('_')[7:9]
                    datetime_fy3d = datetime.strptime(''.join(ymdhm_fy3d), '%Y%m%d%H%M')

                    ymdhms_fy4a = os.path.basename(fy4a_file).split('_')[9][0:12]
                    datetime_fy4a = datetime.strptime(ymdhms_fy4a, '%Y%m%d%H%M')

                    datetime_start = datetime_fy3d - relativedelta(minutes=5)
                    datetime_end = datetime_fy3d + relativedelta(minutes=5)
                    if not (datetime_start <= datetime_fy4a <= datetime_end):
                        continue

                    match_dir = os.path.join(AOD_MATCH_DIR, "FY3D_1KM_FY4A_4KM")
                    make_sure_path_exists(match_dir)
                    out_file = get_out_file(fy3d_file, fy4a_file, match_dir)
                    if os.path.isfile(out_file):
                        print('already exist {}'.format(out_file))
                        continue

                    geo_fy3d_name = 'FY3D_MERSI_GBAL_L1_{}_{}_GEO1K_MS.HDF'.format(ymdhm_fy3d[0], ymdhm_fy3d[1])
                    fy3d_geo_file = os.path.join(GEO_FY3D_1KM_DIR, geo_fy3d_name)

                    match_fy3d_1km_fy4a_4km(fy3d_file, fy3d_geo_file, fy4a_file, out_file,
                                            longitude_range=LONGITUDE_RANGE, latitude_range=LATITUDE_RANGE)


def t():
    aod_fy3d_1km_file = r"/DISK/DATA02/PROJECT/SourceData/ShangHai/AOD/2FY3D_MERSI_L2_AOD-ORBT-1000M/FY3D_MERSI_ORBT_L2_AOD_MLT_NUL_20190520_0810_1000M_MS.HDF"
    geo_fy3d_1km_file = r"/DISK/DATA02/PROJECT/SourceData/ShangHai/AOD/2FY3D-MERSI-L1-GEO1K/FY3D_MERSI_GBAL_L1_20190520_0810_GEO1K_MS.HDF"
    aod_modis_10km_file = r"/DISK/DATA02/PROJECT/SourceData/ShangHai/AOD/4FY4A_AGRI_L2_AOD_4000M/05/20/SATE_FY4A_AGRI_MULT_4000M_DISK_LDA_GLL_HYT_20190520080000_VPRJ5_L2.NC"
    aod_fy3d_5km_modis_10km_file = 'test/test.csv'
    longitude_range = [70, 140]
    latitude_range = [15, 56]
    match_fy3d_1km_fy4a_4km(aod_fy3d_1km_file, geo_fy3d_1km_file, aod_modis_10km_file, aod_fy3d_5km_modis_10km_file,
                            longitude_range=longitude_range, latitude_range=latitude_range)


if __name__ == '__main__':
    match_file()
    """
    
    """