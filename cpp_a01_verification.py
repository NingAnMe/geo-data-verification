#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-07-23 10:20
# @Author  : NingAnMe <ninganme@qq.com>
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta
import pandas as pd
import numpy as np

from lib.cpp import CppFy3c, CppModis
from lib.plot_stats import plot_regression
from lib.verification import Verification

# 单数据测试
fy3c_cpp_file = os.path.join('test', 'fy3c_cpp', 'FY3C_VIRRD_ORBT_L2_CPP_MLT_NUL_20200330_0820_1000M_MS.HDF')  # 本地测试路径
fy3c_geo_file = os.path.join('test', 'fy3c_geo', 'FY3C_VIRRX_GBAL_L1_20200330_0820_GEOXX_MS.HDF')  # 本地测试路径
modis_cpp_file = os.path.join('test', 'modis_cpp', 'MOD06_L2.A2020090.0830.061.2020091134740.hdf')  # 本地测试路径


# 批量测试
fy3c_cpp_dir = r'test/fy3c_cpp'  # 本地测试路径
fy3c_geo_dir = 'test/fy3c_geo'  # 本地测试路径
modis_cpp_dir = 'test/modis_cpp'  # 本地测试路径
result_dir = 'test/result'  # 本地测试路径

# 服务器路径
# fy3c_cpp_dir = '/home/kts_project_v1/qiuh/mod_cpp/20200104/20200104'
# fy3c_geo_dir = '/DISK/DATA02/PROJECT/SourceData/FENGYUN-3C/VIRR/L1/ORBIT/20200104'
# modis_cpp_dir = '/home/kts_project_v1/qiuh/mod_cpp/modis_cpp_20200104/MOD06_L2/2020/004'
# result_dir = '/home/kts_project_v1/qiuh/mod_cpp/result/20200104'


def print_info(data):
    print(np.nanmin(data), np.nanmax(data), np.nanmean(data))


def save_result(result, out_file):
    r = pd.DataFrame(result)
    r = r.dropna(axis=0)
    if len(r) > 0:
        r.to_csv(out_file)


def verification(fy3c_cpp_file, fy3c_geo_file, modis_cpp_file):
    # 获取数据1
    cpp1 = CppFy3c(in_file=fy3c_cpp_file, geo_file=fy3c_geo_file)
    lons1, lats1 = cpp1.get_lon_lat()
    c_tmp1 = cpp1.get_ctop_temperature()
    print_info(lons1)
    print_info(lats1)
    print_info(c_tmp1)

    # 获取数据2
    cpp2 = CppModis(in_file=modis_cpp_file, geo_file=modis_cpp_file)
    lons2, lats2 = cpp2.get_lon_lat()
    c_tmp2 = cpp2.get_ctop_temperature()
    print_info(lons2)
    print_info(lats2)
    print_info(c_tmp2)

    # data2 KDtree建模
    verif = Verification(lons1, lats1, lons2, lats2)
    if not verif.get_kdtree():
        return

    verif.get_dist_and_index_kdtree()

    # 剔除距离差距过大的点
    pre_dist = 0.01
    index_dist = verif.get_index_dist(pre_dist=pre_dist)

    if index_dist.sum() > 0:
        result = {
            'lon_s1': verif.get_kdtree_data(lons1)[index_dist],
            'lat_s1': verif.get_kdtree_data(lats1)[index_dist],
            'tmp_s1': verif.get_kdtree_data(c_tmp1)[index_dist],
            'lon_s2': verif.get_query_data(lons2)[index_dist],
            'lat_s2': verif.get_query_data(lats2)[index_dist],
            'tmp_s2': verif.get_query_data(c_tmp2)[index_dist],
        }
        return result
    else:
        return


def multi_verification():
    fy3c_cpp_filenames = os.listdir(fy3c_cpp_dir)
    fy3c_cpp_filenames = [i for i in fy3c_cpp_filenames if i[-3:] == 'HDF']

    modis_cpp_filenames = os.listdir(modis_cpp_dir)

    pairs = list()

    for fy3c_cpp_filename in fy3c_cpp_filenames:
        ymd1 = '20200330'
        hm1 = fy3c_cpp_filename.split('_')[8]
        ymdhm1 = ymd1 + hm1
        dt1 = datetime.strptime(ymdhm1, '%Y%m%d%H%M')
        dt_s = dt1 - relativedelta(minutes=10)
        dt_e = dt1 + relativedelta(minutes=10)
        for modis_cpp_filename in modis_cpp_filenames:
            hm2 = modis_cpp_filename.split('.')[2]
            ymdhm2 = ymd1 + hm2
            dt2 = datetime.strptime(ymdhm2, '%Y%m%d%H%M')
            if dt_s <= dt2 <= dt_e:
                cpp_file1 = os.path.join(fy3c_cpp_dir, fy3c_cpp_filename)
                geo_file1 = os.path.join(fy3c_geo_dir, f'FY3C_VIRRX_GBAL_L1_{ymd1}_{hm1}_GEOXX_MS.HDF')
                cpp_file2 = os.path.join(modis_cpp_dir, modis_cpp_filename)
                result_file = os.path.join(result_dir, f"FY3C+MERSI_{ymdhm1}_TREEA+MODIS_{ymdhm2}.HDF")
                # print(fy3c_cpp_filename, modis_cpp_filename)
                pairs.append((cpp_file1, geo_file1, cpp_file2, result_file))
    # print(pairs)
    # print(len(pairs))
    # pairs.sort()
    for fy3c_cpp_file, fy3c_geo_file, modis_cpp_file, result_file in pairs:
        result = verification(fy3c_cpp_file, fy3c_geo_file, modis_cpp_file)
        if result is not None:
            save_result(result, result_file)
            print(result_file)


def plot_cpp_regression():
    result_files = os.listdir(result_dir)
    for filename in result_files:
        result_file = os.path.join(result_dir, filename)
        result_data = pd.read_csv(result_file)
        result_data = pd.DataFrame(result_data)
        print(type(result_data))
        print(result_data.head())
        result_data = result_data.dropna(axis=0)
        if len(result_data) <= 10:
            print(f'数据数量小于10，无法绘图')
            continue
        print(result_data.head(2))

        title = os.path.splitext(filename)[0]
        x_label = 'FY3C+VIRR (K)'
        y_label = 'TERRA+MODIS (K)'
        x_range = [30, 170]
        y_range = [30, 170]
        out_file = os.path.join('test', 'pictrue', filename+'.PNG')
        plot_regression(
            x=result_data['tmp_s1'].to_numpy(),
            y=result_data['tmp_s2'].to_numpy(),
            out_file=out_file,
            title=title,
            x_label=x_label,
            y_label=y_label,
            x_range=x_range,
            y_range=y_range,
        )


def plot_delta_regression():
    result_files = os.listdir(result_dir)
    for filename in result_files:
        result_file = os.path.join(result_dir, filename)
        result_data = pd.read_csv(result_file)
        result_data = pd.DataFrame(result_data)
        print(type(result_data))
        print(result_data.head())
        result_data = result_data.dropna(axis=0)
        if len(result_data) <= 10:
            print(f'数据数量小于10，无法绘图')
            continue
        print(result_data.head(2))

        title = os.path.splitext(filename)[0]
        x_label = 'FY3C+VIRR (K)'
        y_label = 'MODIS-VIRR (K)'
        x_range = [30, 170]
        y_range = [-30, 30]
        out_file = os.path.join('test', 'pictrue', filename+'_Delta.PNG')
        x = result_data['tmp_s1'].to_numpy()
        y = x - result_data['tmp_s2'].to_numpy()
        plot_regression(
            x=x,
            y=y,
            out_file=out_file,
            title=title,
            x_label=x_label,
            y_label=y_label,
            x_range=x_range,
            y_range=y_range,
        )


def cal_mean_rate():
    result_files = os.listdir(result_dir)
    for filename in result_files:
        result_file = os.path.join(result_dir, filename)
        result_data = pd.read_csv(result_file)
        result_data = pd.DataFrame(result_data)
        print(type(result_data))
        print(result_data.head())
        result_data = result_data.dropna(axis=0)
        if len(result_data) <= 10:
            print(f'数据数量小于10，无法绘图')
            continue
        print(result_data.head(2))
        d1_mean = result_data['tmp_s1'].to_numpy().mean()
        d2_mean = result_data['tmp_s2'].to_numpy().mean()
        r = (d1_mean - d2_mean) / d2_mean
        out_file = os.path.join('test', 'pictrue', filename+'_Mean.txt')
        with open(out_file, 'w') as fp:
            fp.write(str(r))
            print('>>> {}'.format(out_file))


if __name__ == '__main__':
    verification(fy3c_cpp_file, fy3c_geo_file, modis_cpp_file)  # 单个数据的匹配接口
    multi_verification()  # 批量匹配的接口
    plot_cpp_regression()  # 绘制回归图
    plot_delta_regression()  # 绘制偏差图
    cal_mean_rate()  # 计算平均偏差率
