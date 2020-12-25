#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-09-16 11:09
# @Author  : NingAnMe <ninganme@qq.com>
"""
绘制match结果的图像：regression，hist，timeseries，map_r2
"""
import argparse
from datetime import datetime
from dateutil.relativedelta import relativedelta
import os
from collections import defaultdict

import pandas as pd
import numpy as np
from scipy import stats

from lib.plot_stats import plot_regression, plot_timeseries, plot_histogram
from lib.scical import rmse
from lib.path import make_sure_path_exists
from lib.proj import ProjCore, meter2degree
from lib.province_mask import get_province_mask
from lib.proj_aod import proj_china

from aod_p02_plot_map_origin import plot_map_picture

from config import AOD_PICTURE_DIR, AOD_FY3D_1KM_FY4A_4KM_DIR, AOD_FY3D_1KM_MODIS_3KM_DIR, AOD_FY3D_5KM_MODIS_10KM_DIR
from config import get_areas, get_area_range

import warnings
warnings.filterwarnings('ignore')


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


def get_season(ym):
    season = {
        '201812': '2018 DJF',
        '201903': '2019 MAM',
        '201906': '2019 JJA',
        '201909': '2019 SON',
        '201912': '2019 DJF',
        '202003': '2020 MAM',
        '202006': '2020 JJA',
        '202009': '2020 SON',
        '202012': '2020 DJF',
        '202103': '2021 MAM',
        '202106': '2021 JJA',
        '202109': '2021 SON',
        '202112': '2021 DJF',
    }
    return season[ym]


def plot_verification_picture(date_str, date_end=None, frequency='daily'):
    # 获取数据
    if frequency == 'monthly':
        date_str = date_str[:6]
    elif frequency == 'seasonly':
        date_str = date_str[:6]

    picture_dir = os.path.join(AOD_PICTURE_DIR, 'STAT', MATCH)
    out_dir = os.path.join(picture_dir, 'REGRESSION', frequency)
    out_file = os.path.join(out_dir, 'regression_{}_{}_{}.png'.format(AREA, frequency, date_str))
    # if os.path.isfile(out_file):
    #     print('already exist {}'.format(out_file))
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
    if frequency == 'seasonly':
        season = get_season(date_str)
        title = '{} AOD (550nm) over {}'.format(season, AREA)
    elif frequency == 'all':
        title = '{}-{} AOD (550nm) over {}'.format(date_str, date_end, AREA)
    else:
        title = "{} AOD (550nm) over {}".format(date_str, AREA)

    x_range = [0, 1.5]
    y_range = [0, 1.5]
    x_label = '{}'.format(pair_x)
    y_label = '{}'.format(pair_y)
    x_interval = 0.3
    y_interval = 0.3

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
    if frequency != 'daily':
        density = True
    else:
        density = False

    plot_regression(
        x=x,
        y=y,
        out_file=out_file,
        title=title,
        x_label=x_label,
        y_label=y_label,
        x_range=x_range,
        y_range=y_range,
        x_interval=x_interval,
        y_interval=y_interval,
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


def plot_verification_picture_map(date_str, date_end=None, frequency='daily'):
    # 获取数据
    if frequency == 'monthly':
        date_str = date_str[:6]
    elif frequency == 'seasonly':
        date_str = date_str[:6]

    picture_dir = os.path.join(AOD_PICTURE_DIR, 'STAT', MATCH)
    out_dir = os.path.join(picture_dir, 'R_MAP', frequency)
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

    # 投影
    print('投影')
    ii, jj = proj.lonslats2ij(lons, lats)

    # 网格数据统计
    print('网格数据统计')
    data_dict = defaultdict(list)
    for i, j, x, y in zip(ii, jj, x, y):
        data_dict[(i, j)].append((x, y))

    # 新的网格的经纬度
    lons_grid, lats_grid = proj.grid_lonslats()

    # 绘图
    print('绘图')
    lats_plot = list()
    lons_plot = list()
    r2_plot = list()
    r2_grid = np.full_like(lats_grid, np.nan, dtype=np.float)
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
            r2_plot.append(np.abs(r_value))
            r2_grid[i, j] = np.abs(r_value)
    if not lats_plot:
        print('没有数据： {}'.format(file_out))
        return

    if frequency == 'seasonly':
        season = get_season(date_str)
        title = '{} R over {}'.format(season, AREA)
    elif frequency == 'all':
        title = '{}-{} R over {}'.format(date_str, date_end, AREA)
    else:
        title = "{} R over {}".format(date_str, AREA)
    vmin = 0
    vmax = 1

    ticks = np.arange(0, 1.1, 0.2)
    if AREA == 'China':
        mksize = 10
        nanhai = True
    else:
        mksize = 20
        nanhai = False

    data, lons, lats = proj_china(r2_grid, lons_grid, lats_grid, data_min=0, data_max=1)

    areas = get_areas(AREA)
    mask = get_province_mask(areas)

    valid = np.logical_and.reduce((data > vmin, data < vmax, mask))

    data_mask = data[valid]
    lons_mask = lons[valid]
    lats_mask = lats[valid]

    count = len(data_mask)
    print('count == {}'.format(count))

    box = [LATITUDE_RANGE[1], LATITUDE_RANGE[0], LONGITUDE_RANGE[0],
           LONGITUDE_RANGE[1]]  # nlat, slat, wlon, elon:北（大），南（小），西（小），东（大）

    plot_map_picture(data_mask, lons_mask, lats_mask, title=title, vmin=vmin, vmax=vmax,
                     areas=areas, box=box, ticks=ticks, file_out=file_out,
                     mksize=mksize, nanhai=nanhai)


def multi_plot_regression(date_start, date_end, frequency='daily'):
    datetime_start = datetime.strptime(date_start, "%Y%m%d")
    datetime_end = datetime.strptime(date_end, "%Y%m%d")
    stats_data = defaultdict(list)
    while datetime_start <= datetime_end:
        date_str = datetime_start.strftime("%Y%m%d")
        datas = plot_verification_picture(date_str, date_end=date_end, frequency=frequency)
        if datas:
            for key, value in datas.items():
                stats_data[key].append(value)

        if frequency == 'daily':
            datetime_start = datetime_start + relativedelta(days=1)
        elif frequency == 'monthly':
            datetime_start = datetime_start + relativedelta(months=1)
        elif frequency == 'seasonly':
            datetime_start = datetime_start + relativedelta(months=3)
        elif frequency == 'yearly':
            datetime_start = datetime_start + relativedelta(years=1)
        else:
            break
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
        plot_verification_picture_map(date_str, date_end=date_end, frequency=frequency)
        if frequency == 'daily':
            datetime_start = datetime_start + relativedelta(days=1)
        elif frequency == 'monthly':
            datetime_start = datetime_start + relativedelta(months=1)
        elif frequency == 'seasonly':
            datetime_start = datetime_start + relativedelta(months=3)
        elif frequency == 'yearly':
            datetime_start = datetime_start + relativedelta(years=1)
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
    try:
        data = get_stats_data(date_start, date_end, frequency)
    except pd.errors.EmptyDataError as why:
        print(why)
        return

    x = list()
    for i in data.date:
        x.append(datestr2datetime(i))

    picture_dir = os.path.join(AOD_PICTURE_DIR, 'STAT', MATCH)
    out_dir = os.path.join(picture_dir, 'TIMESERIES')

    # ================================ plot BIAS
    y = data.bias_mean
    y_label = 'BIAS'
    out_file = os.path.join(out_dir, 'timeseries_{}_{}_BIAS_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    title = '{}-{} BIAS AOD over {}'.format(date_start, date_end, AREA)
    y_range = [-0.5, 0.5]
    y_interval = 0.1
    # if not os.path.isfile(out_file):
    #     plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True)
    # else:
    #     print('already exist {}'.format(out_file))
    #     return
    plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True,
                    y_interval=y_interval,
                    ymd_start=date_start, ymd_end=date_end)

    # ================================ plot r
    y = data.R
    y_label = "R"
    out_file = os.path.join(out_dir, 'timeseries_{}_{}_R_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    title = '{}-{} R over {}'.format(date_start, date_end, AREA)
    y_range = [0, 1]
    # if not os.path.isfile(out_file):
    #     plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True)
    # else:
    #     print('already exist {}'.format(out_file))
    #     return
    plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True,
                    ymd_start=date_start, ymd_end=date_end)

    # ================================ plot RMSE
    y = data.RMSE
    y_label = "RMSE"
    out_file = os.path.join(out_dir, 'timeseries_{}_{}_RMSE_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    title = '{}-{} RMSE over {}'.format(date_start, date_end, AREA)
    y_range = [0, 1]
    # if not os.path.isfile(out_file):
    #     plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True)
    # else:
    #     print('already exist {}'.format(out_file))
    #     return
    plot_timeseries(x, y, out_file=out_file, title=title, y_label=y_label, y_range=y_range, plot_month=True,
                    ymd_start=date_start, ymd_end=date_end)

    # ++++++++++++++++++++++++++++++++ plot BIAS Hist
    y = data.bias_mean
    x_range = [-0.5, 0.5]
    y_range = [0, 50]
    title = '{}-{} Hist BIAS AOD over {}'.format(date_start, date_end, AREA)
    x_label = 'BIAS'
    y_label = 'Count'
    out_file = os.path.join(out_dir, 'histogram_{}_{}_BIAS_{}_{}.png'.format(AREA, frequency, date_start, date_end))
    # if not os.path.isfile(out_file):
    #     plot_histogram(data=y, out_file=out_file, bins_count=20, title=title, x_label=x_label,
    #                    x_range=x_range, y_range=y_range, )
    # else:
    #     print('already exist {}'.format(out_file))
    #     return
    plot_histogram(data=y, out_file=out_file, bins_count='auto', title=title, x_label=x_label, y_label=y_label,
                   x_range=x_range, y_range=y_range, )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='help')
    parser.add_argument('--matchType', '-m', help="匹配对：'AOD_FY3D_1KM_MODIS_3KM', 'AOD_FY3D_5KM_MODIS_10KM', 'AOD_FY3D_1KM_FY4A_4KM'", required=True)
    parser.add_argument('--dateType', help='时间类型（日、月、季、年、全部）：Daily、Monthly、Seasonly、Yearly、 All', required=False)
    parser.add_argument('--plotType', help='绘图类型（回归图、R分布图、时间序列图）：regression、map、timeseries', required=True)
    parser.add_argument('--areaType', help='地区类型：China、YRD、PRD、FWP、BTH', required=True)
    args = parser.parse_args()

    MATCH = args.matchType
    if args.dateType is not None:
        dateType = args.dateType.lower()
    else:
        dateType = None
    plotType = args.plotType.lower()
    areaType = args.areaType

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
        raise ValueError(MATCH)

    AREAs = ['China', 'YRD', 'PRD', 'FWP', 'BTH']

    if areaType is not None:
        assert areaType in AREAs
        AREAs = [areaType]

    for AREA in AREAs:
        LONGITUDE_RANGE, LATITUDE_RANGE = get_area_range(AREA)

        if plotType == 'regression':
            if dateType == 'daily':
                multi_plot_regression(_date_start, _date_end, 'daily')
            elif dateType == 'monthly':
                multi_plot_regression(_date_start, _date_end, 'monthly')
            elif dateType == 'seasonly':
                multi_plot_regression('20181201', _date_end, 'seasonly')
            elif dateType == 'yearly':
                multi_plot_regression(_date_start, _date_end, 'yearly')
            elif dateType == 'all':
                multi_plot_regression(_date_start, _date_end, 'all')
            else:
                raise ValueError(dateType)
        elif plotType == 'map':
            if dateType == 'monthly':
                multi_plot_map(_date_start, _date_end, 'monthly')
            elif dateType == 'seasonly':
                multi_plot_map('20181201', _date_end, 'seasonly')
            elif dateType == 'yearly':
                multi_plot_map(_date_start, _date_end, 'yearly')
            elif dateType == 'all':
                multi_plot_map(_date_start, _date_end, 'all')
            else:
                raise ValueError(dateType)
        elif plotType == 'timeseries':
            plot_timeseries_picture(_date_start, _date_end, 'daily')
        else:
            raise ValueError(plotType)
