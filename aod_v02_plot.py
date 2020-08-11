#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-08-11 15:05
# @Author  : NingAnMe <ninganme@qq.com>

import os

import numpy as np
import pandas as pd
from scipy import stats

from lib.plot import plot_regression

fy3d_aeronet_dir = r'/RED1BDATA/cma/AEROSOL_1.0/SupportData/FY3D_MERSI_AERONET'
fy3d_aeronet_image_dir = r'/RED1BDATA/cma/AEROSOL_1.0/SupportData/FY3D_MERSI_AERONET_IMAGE'


def plot_aod_regression(ymd):
    file_dir = os.path.join(fy3d_aeronet_dir, ymd)
    result_files = os.listdir(file_dir)
    result_data_all = None
    for filename in result_files:
        result_file = os.path.join(file_dir, filename)
        result_data = pd.read_csv(result_file)
        # print(result_data.head())
        result_data = result_data.dropna(axis=0)
        if result_data_all is None:
            result_data_all = result_data
        else:
            result_data_all = pd.concat((result_data_all, result_data), axis=0)

    print(result_data_all)
    if result_data_all is None or len(result_data_all) <= 10:
        print(f'数据数量小于10，无法绘图')
        return
    print(result_data_all.head(2))

    title = ymd
    x_label = 'FY3D+MERSI'
    y_label = 'AERONET'
    x_range = [0, 1.2]
    y_range = [0, 1.2]
    # x_interval = 0.2
    # y_interval = 0.2
    out_file = os.path.join('test', 'pictrue', ymd+'.PNG')

    x = result_data_all['aod_s1'].to_numpy()
    y = result_data_all['aod_s2'].to_numpy()

    index = np.logical_and.reduce((x > 0, y > 0))
    x = x[index]
    y = y[index]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print(slope, intercept, r_value, r_value ** 2, p_value, std_err)
    plot_regression(
        x=x,
        y=y,
        out_file=out_file,
        title=title,
        x_label=x_label,
        y_label=y_label,
        x_range=x_range,
        y_range=y_range,
        # x_interval=x_interval,
        # y_interval=y_interval,
    )


def plot_aod_regression_month(ym):
    filename_list = os.listdir(fy3d_aeronet_dir)
    result_data_all = None
    for filename in filename_list:
        if ym in filename:
            result_file = os.path.join(fy3d_aeronet_dir, filename)
            result_data = pd.read_csv(result_file)
            result_data = result_data.dropna(axis=0)
            if result_data_all is None:
                result_data_all = result_data
            else:
                result_data_all = pd.concat((result_data_all, result_data), axis=0)

    print(result_data_all)
    if result_data_all is None or len(result_data_all) <= 10:
        print(f'数据数量小于10，无法绘图')
        return
    print(result_data_all.head(2))

    title = ym
    x_label = 'FY3D+MERSI'
    y_label = 'AERONET'
    x_range = [0, 1.2]
    y_range = [0, 1.2]
    # x_interval = 0.2
    # y_interval = 0.2
    out_file = os.path.join(fy3d_aeronet_image_dir, ym+'.PNG')

    x = result_data_all['aod_s1'].to_numpy()
    y = result_data_all['aod_s2'].to_numpy()

    index = np.logical_and.reduce((x > 0, y > 0))
    x = x[index]
    y = y[index]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    print(slope, intercept, r_value, r_value ** 2, p_value, std_err)
    plot_regression(
        x=x,
        y=y,
        out_file=out_file,
        title=title,
        x_label=x_label,
        y_label=y_label,
        x_range=x_range,
        y_range=y_range,
        # x_interval=x_interval,
        # y_interval=y_interval,
    )


def main():
    ym = '201902'
    plot_aod_regression_month(ym)


if __name__ == '__main__':
    main()
