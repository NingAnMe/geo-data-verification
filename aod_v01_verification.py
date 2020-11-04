#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-07-30 12:23
# @Author  : NingAnMe <ninganme@qq.com>
import os

import numpy as np
import pandas as pd
from scipy import stats

from lib.aod import AodFy3d
from lib.aeronet import Aeronet
from lib.plot_stats import plot_regression
from lib.verification import Verification

aeronet_map_file = r'site_info.csv'
fy3d_aod_dir = r'/home/kts_project_v1/qiuh/mod_aod/fy3d_aod/Granule'  # gongsi
aeronet_aod_dir = r'/DATA/PROJECT/SourceData/Aeronet/AOD/AOD20/ALL_POINTS'  # gongsi
fy3d_aeronet_dir = r'/home/kts_project_v1/qiuh/mod_aod/fy3d_aeronet'
fy3d_aeronet_image_dir = r'/home/kts_project_v1/qiuh/mod_aod/fy3d_aeronet_image'


# fy3d_aod_dir = r'/nas02/cma/AEROSOL_1.0/SupportData/FY3D_MERSI/Granule'
# fy3d_aod_dir = r'/RED1BDATA/cma/AEROSOL_1.0/SupportData/FY3D_MERSI/Granule/Granule'
# fy3d_aod_file = r'/home/kts_project_v1/qiuh/mod_cpp/fy3d_aod/FY3D_MERSI_AOD_GRANULE_20190228_0000.HDF5'
# aeronet_aod_dir = r'/RED1BDATA/cma/SourceData/Aeronet/AOD/AOD20/ALL_POINTS'  # liujian
# fy3d_aeronet_dir = r'/home/kts_project_v1/qiuh/mod_cpp/fy3d_aeronet'
# fy3d_aeronet_dir = r'/RED1BDATA/cma/AEROSOL_1.0/SupportData/FY3D_MERSI_AERONET'
# fy3d_aeronet_image_dir = r'/RED1BDATA/cma/AEROSOL_1.0/SupportData/FY3D_MERSI_AERONET_IMAGE'


def print_info(data):
    print(np.nanmin(data), np.nanmax(data), np.nanmean(data))


def verification(fy3d_aod_file, aeronet_map, aeronet_file_dir):
    print(f"<<< {fy3d_aod_file}")
    fy3d_aod_name = os.path.splitext(os.path.basename(fy3d_aod_file))[0]
    out_file = os.path.join(fy3d_aeronet_dir, fy3d_aod_name + '.csv')
    if os.path.exists(out_file):
        print('输出文件已经存在：{}'.format(out_file))
        return
    # 获取数据1
    aod1 = AodFy3d(in_file=fy3d_aod_file, geo_file=fy3d_aod_file)
    c_aod1 = aod1.get_aod()
    valid_index = np.logical_and(c_aod1 > 0, c_aod1 < 10)
    lons1, lats1 = aod1.get_lon_lat()
    dt1 = aod1.dt
    c_aod1 = c_aod1[valid_index]
    lons1 = lons1[valid_index]
    lats1 = lats1[valid_index]

    # 获取数据2
    aod2 = pd.read_csv(aeronet_map, index_col=False)
    aod2.sort_values('dt_e', inplace=True, ascending=False)
    lons2 = aod2['lon'].to_numpy()[:333]
    lats2 = aod2['lat'].to_numpy()[:333]
    name = aod2['name'].to_numpy()[:333]

    # 数据匹配
    verif = Verification(lons1, lats1, lons2, lats2)
    if not verif.get_kdtree():
        return

    verif.get_dist_and_index_kdtree()

    # 剔除距离差距过大的点
    pre_dist = 0.1
    index_dist = verif.get_index_dist(pre_dist=pre_dist)

    # 匹配数据
    if index_dist.sum() > 0:
        result = {
            'lons_s1': verif.get_kdtree_data(lons1)[index_dist],
            'lats_s1': verif.get_kdtree_data(lats1)[index_dist],
            'aod_s1': verif.get_kdtree_data(c_aod1)[index_dist],
            'dt_s1': [dt1] * index_dist.sum(),
            'lons_s2': verif.get_query_data(lons2)[index_dist],
            'lats_s2': verif.get_query_data(lats2)[index_dist],
            'name': verif.get_query_data(name)[index_dist],
            'dist': verif.dist[index_dist],
        }
    else:
        print('匹配的数据量 < 0')
        return

    # 循环匹配到的站点，找到时间最接近的点
    station_name = result['name']
    #   时间阈值
    pre_dts = pd.Timedelta('00:30:00')
    dt2 = list()
    aod2 = list()
    for name in station_name:
        for aeronet_filename in os.listdir(aeronet_file_dir):
            name2 = '_'.join(os.path.splitext(aeronet_filename)[0].split('_')[2:])
            if name == name2:
                aeronet = Aeronet(os.path.join(aeronet_file_dir, aeronet_filename))
                dts = aeronet.get_datetime()
                dts = pd.to_datetime(dts, unit='s')
                dts_delta = np.abs(dts - dt1)
                index_dt = np.argmin(dts_delta)
                if dts_delta[index_dt] < pre_dts:
                    aod2.append(aeronet.get_aod550()[index_dt])
                else:
                    aod2.append(np.nan)
                dt2.append(dts[index_dt])
    result['aod_s2'] = aod2
    result['dt_s2'] = dt2
    out_data_df = pd.DataFrame(result)

    out_data_df.to_csv(out_file)
    print(out_data_df)


def main_month():

    aod_dir_list = os.listdir(fy3d_aod_dir)
    aod_dir_list.sort()
    for ymd in aod_dir_list:
        one_day_dir = os.path.join(fy3d_aod_dir, ymd)
        fy3d_aod_filenames = os.listdir(one_day_dir)
        fy3d_aod_filenames.sort()
        for filename in fy3d_aod_filenames:
            fy3d_aod_file = os.path.join(one_day_dir, filename)
            if not os.path.isfile(fy3d_aod_file):
                continue
            if not os.path.splitext(fy3d_aod_file)[1] == '.HDF5':
                continue
            verification(fy3d_aod_file=fy3d_aod_file, aeronet_map=aeronet_map_file, aeronet_file_dir=aeronet_aod_dir)


def main_day():
    ymd = '20190228'
    one_day_dir = os.path.join(fy3d_aod_dir, ymd)
    fy3d_aod_filenames = os.listdir(one_day_dir)
    fy3d_aod_filenames.sort()
    for filename in fy3d_aod_filenames:
        fy3d_aod_file = os.path.join(one_day_dir, filename)
        if not os.path.isfile(fy3d_aod_file):
            continue
        if not os.path.splitext(fy3d_aod_file)[1] == '.HDF5':
            continue
        verification(fy3d_aod_file=fy3d_aod_file, aeronet_map=aeronet_map_file, aeronet_file_dir=aeronet_aod_dir)


if __name__ == '__main__':
    # main_day()
    main_month()
