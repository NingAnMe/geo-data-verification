#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-09-16 11:09
# @Author  : NingAnMe <ninganme@qq.com>
"""
绘制match结果的图像：regression，hist，timeseries，map_r2
"""
from datetime import datetime
from dateutil.relativedelta import relativedelta
import os
from collections import defaultdict

import pandas as pd
import numpy as np
from scipy import stats

from lib.plot_stats import plot_regression, plot_timeseries, plot_histogram
from lib.plot_map import plot_shanghai
from lib.scical import rmse
from lib.path import make_sure_path_exists
from lib.proj import ProjCore, meter2degree

from config import *


def get_match_data(date_str_in, frequency='daily'):
    """
    找到对应日期或者对应月份的匹配数据，并且读取所有数据
    :param date_str_in: str
    :param frequency: daily or monthly
    :return: pandas.DataFrame
    """
    match_dir = os.path.join(RESULT_DIR, 'MATCH')
    print("<<< {}".format(match_dir))
    if frequency == 'all':
        filenames = os.listdir(match_dir)
        files = [os.path.join(match_dir, filename) for filename in filenames]
    else:
        if frequency == 'monthly':
            date_str_in = date_str_in[:6]
        elif frequency == 'seasonly':
            date_str_in = date_str_in[:6]

        match_file_dict = defaultdict(list)

        for root, dirs, files in os.walk(match_dir):
            for name in files:
                if name[-3:].lower() != 'csv':
                    continue
                date_str = name.split('_')[1]
                if frequency == 'monthly':
                    date_str = date_str[:6]
                elif frequency == 'seasonly':
                    date_str = date_str[:6]
                match_file_dict[date_str].append(os.path.join(root, name))

        if date_str_in not in match_file_dict:
            print('此日期没有数据：{}'.format(date_str_in))

        if frequency == 'seasonly':
            files = list()
            for months_delta in range(3):
                datetime_start = datetime.strptime(date_str_in, "%Y%m")
                datetime_start += relativedelta(months=months_delta)
                date_str = datetime_start.strftime("%Y%m")
                if date_str in match_file_dict:
                    files.extend(match_file_dict[date_str])
        else:
            files = match_file_dict[date_str_in]

    data = None
    for match_file in files:
        print('<<< {}'.format(match_file))
        data_tmp = pd.read_csv(match_file)

        # 过滤无效数据
        data_tmp = data_tmp[
            (data_tmp.aod_y > 0) & (data_tmp.aod_y < 1.5) & (data_tmp.aod_x > 0) & (data_tmp.aod_x < 1.5)]

        # 范围过滤
        if LONGITUDE_RANGE and LATITUDE_RANGE:
            data_tmp = data_tmp[(data_tmp.lons_x > LONGITUDE_RANGE[0]) & (data_tmp.lons_x < LONGITUDE_RANGE[1]) &
                                (data_tmp.lats_x > LATITUDE_RANGE[0]) & (data_tmp.lats_x < LATITUDE_RANGE[1])]
        if data_tmp.empty:
            continue

        if data is None:
            data = data_tmp
        else:
            data = pd.concat((data, data_tmp), axis=0)
    return data


def plot_verification_picture(date_str, frequency='daily'):
    # 获取数据
    if frequency == 'monthly':
        date_str = date_str[:6]
    elif frequency == 'sensonly':
        date_str = date_str[:6]

    picture_dir = os.path.join(RESULT_DIR, 'PICTURE', AREA)
    if frequency == 'all':
        out_dir = os.path.join(RESULT_DIR, 'PICTURE', 'TIMESEIES')
        out_file = os.path.join(out_dir, 'regression_{}_{}_{}_{}.png'.format(AREA, frequency, _date_start, _date_end))
    else:
        out_dir = os.path.join(picture_dir, 'REGRESSION', frequency)
        out_file = os.path.join(out_dir, 'regression_{}_{}_{}.png'.format(AREA, frequency, date_str))
    # if os.path.isfile(out_file):
    #     return

    data = get_match_data(date_str, frequency=frequency)
    if data is None:
        print('没有获取到任何数据：{}'.format(date_str))
        return
    print('输入的样本数量: {}'.format(len(data)))

    # 获取x、y数据和信息
    x = data.aod_x
    y = data.aod_y

    count = len(x)
    if count < 50:
        print('数据量小于50: {}'.format(count))
        return

    # ====================== 绘制回归图 ===========================
    title = "{} {} AOD 550nm".format(date_str, AREA)

    x_range = [0, 1.5]
    y_range = [0, 1.5]
    x_label = '{}'.format(pair_x)
    y_label = '{}'.format(pair_y)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r = r_value  # 2020-11-04：散点图改为使用r，不使用r2
    r2 = r_value ** 2
    _rmse = rmse(x, y)
    x_mean = x.mean()
    y_mean = y.mean()
    bias_mean = (x - y).mean()
    print(slope, intercept, r_value, r, p_value, std_err, _rmse)
    annotate = {"left_top": ["Count :{:d}".format(count),
                             "Slope :{:0.2f}".format(slope),
                             "Intercept :{:0.2f}".format(intercept),
                             "R :{:0.2f}".format(r),
                             "RMSE :{:0.2f}".format(_rmse),
                             ]}
    if frequency == 'all':
        density = True
    else:
        density = True

    plot_regression(
        x=x,
        y=y,
        out_file=out_file,
        title=title,
        x_label=x_label,
        y_label=y_label,
        x_range=x_range,
        y_range=y_range,
        annotate=annotate,
        density=density
    )

    stats_data = {
        'date': date_str,
        'count': count,
        'slope': slope,
        'intercept': intercept,
        'r2': r2,
        'R': r,
        'RMSE': _rmse,
        'x_mean': x_mean,
        'y_mean': y_mean,
        'bias_mean': bias_mean,
    }
    return stats_data


def plot_verification_picture_map(date_str, frequency='daily'):
    # 获取数据
    if frequency == 'monthly':
        date_str = date_str[:6]
    elif frequency == 'sensonly':
        date_str = date_str[:6]

    out_dir = os.path.join(RESULT_DIR, 'PICTURE', 'TIMESEIES')
    file_out = os.path.join(out_dir, 'r2_map_{}_{}_{}.png'.format(AREA, frequency, date_str))
    # if os.path.isfile(file_out):
    #     print('already exist {}'.format(file_out))
    #     return

    data = get_match_data(date_str, frequency=frequency)
    if data is None:
        print('没有获取到任何数据：{}'.format(date_str))
        return
    print('输入的样本数量: {}'.format(len(data)))

    # 获取x、y数据和信息
    x = data.aod_x
    y = data.aod_y
    lats = data.lats_x.to_numpy()
    lons = data.lons_x.to_numpy()

    count = len(x)
    if count < 50:
        print('数据量小于50: {}'.format(count))
        return

    # 创建投影查找表
    print('创建投影查找表')
    res_degree = meter2degree(10000)  # 分辨率，10km
    projstr = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    proj = ProjCore(projstr, res_degree, unit="deg", pt_tl=(69.995, 55.995), pt_br=(139.995, 14.995))  # 角点也要放在格点中心位置
    result = proj.create_lut(lats=lats, lons=lons)

    # 投影
    print('投影')
    ii, jj = proj.lonslats2ij(lons, lats)

    # 投影新网格的经纬度信息
    print('投影新网格的经纬度信息')
    result['Longitude'] = proj.lons
    result['Latitude'] = proj.lats

    # 网格数据统计
    print('网格数据统计')
    data_dict = defaultdict(list)
    for i, j, x, y in zip(ii, jj, x, y):
        data_dict[(i, j)].append((x, y))

    # 新的网格的经纬度
    proj.grid_lonslats()
    lats_grid, lons_grid = proj.lats, proj.lons

    # 绘图
    print('绘图')
    lats_plot = list()
    lons_plot = list()
    r2_plot = list()
    r2_grid = np.full_like(lats_grid, 0, dtype=np.float)
    for (i, j), aod_x_y in data_dict.items():
        xs = list()
        ys = list()
        for x, y in aod_x_y:
            xs.append(x)
            ys.append(y)
        slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)
        lat = lats_grid[i, j]
        lon = lons_grid[i, j]
        if p_value < 0.05:
            lats_plot.append(lat)
            lons_plot.append(lon)
            r2_plot.append(r_value ** 2)
            r2_grid[i, j] = r_value ** 2
    if not lats_plot:
        print('没有数据： {}'.format(file_out))
        return

    if LONGITUDE_RANGE is not None:
        box = [LATITUDE_RANGE[0], LATITUDE_RANGE[1], LONGITUDE_RANGE[0],
               LONGITUDE_RANGE[1]]  # nlat, slat, wlon, elon:北（小），南（大），西（小），东（大）
    else:
        box = None
    title = ''
    vmin = 0
    vmax = 1
    markersize = 20
    lats_plot = np.array(lats_plot)
    lons_plot = np.array(lons_plot)
    r2_plot = np.array(r2_plot)
    plot_shanghai(lats_plot, lons_plot, r2_plot, file_out,
                  box=box, title=title, vmin=vmin, vmax=vmax, markersize=markersize)


def multi_plot_regression(date_start, date_end, frequency='daily'):
    datetime_start = datetime.strptime(date_start, "%Y%m%d")
    datetime_end = datetime.strptime(date_end, "%Y%m%d")
    stats_data = defaultdict(list)
    while datetime_start <= datetime_end:
        date_str = datetime_start.strftime("%Y%m%d")
        datas = plot_verification_picture(date_str, frequency)
        if datas:
            for key, value in datas.items():
                stats_data[key].append(value)

        if frequency == 'daily':
            datetime_start = datetime_start + relativedelta(days=1)
        elif frequency == 'monthly':
            datetime_start = datetime_start + relativedelta(months=1)
        elif frequency == 'seasonly':
            datetime_start = datetime_start + relativedelta(months=3)
    stats_data = pd.DataFrame(stats_data)
    out_dir = os.path.join(RESULT_DIR, 'STATS')
    make_sure_path_exists(out_dir)
    out_file = os.path.join(out_dir, '{}_{}_{}_{}.csv'.format(AREA, frequency, date_start, date_end))
    stats_data.to_csv(out_file, index=False)
    print('>>> {}'.format(out_file))


def multi_plot_map(date_start, date_end, frequency='daily'):
    datetime_start = datetime.strptime(date_start, "%Y%m%d")
    datetime_end = datetime.strptime(date_end, "%Y%m%d")
    while datetime_start <= datetime_end:
        date_str = datetime_start.strftime("%Y%m%d")
        plot_verification_picture_map(date_str, frequency)
        if frequency == 'daily':
            datetime_start = datetime_start + relativedelta(days=1)
        elif frequency == 'monthly':
            datetime_start = datetime_start + relativedelta(months=1)
        elif frequency == 'seasonly':
            datetime_start = datetime_start + relativedelta(months=3)
        else:
            break


# ############################################## timeseries ############################################
def get_stats_data(date_start, date_end, frequency='daily'):
    stats_dir = os.path.join(RESULT_DIR, 'STATS')
    stats_file = os.path.join(stats_dir, '{}_{}_{}_{}.csv'.format(AREA, frequency, date_start, date_end))
    print('<<< {}'.format(stats_file))
    data_tmp = pd.read_csv(stats_file)
    return data_tmp


def datestr2datetime(date_str):
    return datetime.strptime(str(date_str)[:8], "%Y%m%d")


def plot_timeseries_picture(date_start, date_end, frequency='daily'):
    data = get_stats_data(date_start, date_end, frequency)

    x = list()
    for i in data.date:
        x.append(datestr2datetime(i))

    out_dir = os.path.join(RESULT_DIR, 'PICTURE', 'TIMESEIES')

    # ================================ plot BIAS
    y = data.bias_mean
    y_label = 'BIAS'
    out_file = os.path.join(out_dir, 'timeseries_{}_{}_BIAS_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    title = '{} {}-{} {} And {} AOD BIAS'.format(AREA, date_start, date_end, pair_x, pair_y)
    y_range = [-0.5, 0.5]
    # y_range = [0, 1]
    plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True)

    # ================================ plot r
    y = data.R
    y_label = "R"
    out_file = os.path.join(out_dir, 'timeseries_{}_{}_R_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    title = '{} {}-{} {} And {} AOD R'.format(AREA, date_start, date_end, pair_x, pair_y)
    y_range = [-1, 1]
    plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True)

    # ================================ plot RMSE
    y = data.RMSE
    y_label = "RMSE"
    out_file = os.path.join(out_dir, 'timeseries_{}_{}_RMSE_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    title = '{} {}-{} {} And {} AOD RMSE'.format(AREA, date_start, date_end, pair_x, pair_y)
    y_range = [0, 1]
    plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True)

    # ++++++++++++++++++++++++++++++++ plot BIAS Hist
    y = data.bias_mean
    x_range = [-0.5, 0.5]
    # x_range = [0, 1]
    y_range = [0, 50]
    title = '{} {}-{} {} And {} AOD BIAS'.format(AREA, date_start, date_end, pair_x, pair_y)
    x_label = 'BIAS'
    out_file = os.path.join(out_dir, 'histogram_{}_{}_BIAS_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    plot_histogram(data=y, out_file=out_file, bins_count=20, title=title, x_label=x_label,
                   x_range=x_range, y_range=y_range, )


if __name__ == '__main__':
    # AREA = 'China'
    # AREA = 'ChangSanJiao'
    # AREA = 'ZhuSanJiao'
    # AREA = 'FenWei'
    # AREA = 'JingJinJi'

    # MATCH = 'AOD_FY3D_1KM_MODIS_3KM'
    # MATCH = 'AOD_FY3D_5KM_MODIS_10KM'
    # MATCH = 'AOD_FY3D_1KM_FY4A_4KM'

    MATCHES = ['AOD_FY3D_1KM_MODIS_3KM', 'AOD_FY3D_5KM_MODIS_10KM', 'AOD_FY3D_1KM_FY4A_4KM']
    MATCHES = ['AOD_FY3D_1KM_MODIS_3KM', 'AOD_FY3D_5KM_MODIS_10KM']
    for MATCH in MATCHES:

        if MATCH == 'AOD_FY3D_5KM_MODIS_10KM':
            RESULT_DIR = AOD_FY3D_5KM_MODIS_10KM_DIR
            pair_x = 'FY3D MERSI'
            pair_y = 'AQUA MODIS'
            _date_start = "20190101"
            _date_end = "20200531"
        elif MATCH == 'AOD_FY3D_1KM_MODIS_3KM':
            RESULT_DIR = AOD_FY3D_1KM_MODIS_3KM_DIR
            pair_x = 'FY3D MERSI'
            pair_y = 'AQUA MODIS'
            _date_start = "20190101"
            _date_end = "20200531"
        elif MATCH == 'AOD_FY3D_1KM_FY4A_4KM':
            RESULT_DIR = AOD_FY3D_1KM_FY4A_4KM_DIR
            pair_x = 'FY3D MERSI'
            pair_y = 'FY4A AGRI'
            _date_start = "20190101"
            _date_end = "20191231"
        else:
            _date_start = "20190101"
            _date_end = "20200531"

        # AREAs = ['China', 'ChangSanJiao', 'ZhuSanJiao', 'FenWei', 'JingJinJi']
        AREAs = ['ChangSanJiao', 'ZhuSanJiao', 'JingJinJi', 'FenWei']

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

            # multi_plot_regression(_date_start, _date_end, 'daily')
            # multi_plot_regression(_date_start, _date_end, 'monthly')
            # multi_plot_regression('20181201', _date_end, 'seasonly')
            # plot_verification_picture(_date_start, 'all')
            # plot_timeseries_picture(_date_start, _date_end, 'daily')

            # multi_plot_map('20181201', _date_end, 'seasonly')
            multi_plot_map('20190101', _date_end, 'all')
