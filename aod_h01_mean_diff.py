#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-10-29 10:17
# @Author  : NingAnMe <ninganme@qq.com>
"""
计算AOD数据月平均和季度平均
"""
from datetime import datetime
from collections import defaultdict
from scipy.interpolate import griddata

from lib.path import *
from lib.aod import *
from lib.hdf5 import write_hdf5_and_compress, get_hdf5_data
from lib.proj import ProjCore, meter2degree

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


def mean_modis_10kmto5km(frequency='monthly'):
    file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_MODIS_10KM_DIR):
        for name in files:
            if name[-3:].lower() != 'hdf':
                continue

            date_str_file = name[10:17]
            date_ = datetime.strptime(date_str_file, "%Y%j")
            if date_ > datetime(2020, 6, 1):
                continue
            date_str = date_.strftime("%Y%m")

            if frequency == 'monthly':
                file_dict[date_str[:6]].append(os.path.join(root, name))
            elif frequency == 'seasonly':
                season = get_season(date_str[:6])
                if season is not None:
                    file_dict[season].append(os.path.join(root, name))
            elif frequency == 'all':
                file_dict['all'].append(os.path.join(root, name))
            else:
                continue

    # 创建投影查找表
    print('创建投影查找表')
    res_degree = 0.05  # 分辨率，5km
    projstr = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
    proj = ProjCore(projstr, res_degree, unit="deg", pt_tl=(69.995, 55.995), pt_br=(139.995, 14.995))  # 角点也要放在格点中心位置

    for d_, files in file_dict.items():
        aod_sum = np.zeros((proj.row, proj.col), dtype=np.float)
        aod_count = np.zeros_like(aod_sum, dtype=np.float)
        filename_out = 'AOD_MODIS_5KM_MEAN_{}_{}.HDF'.format(frequency, d_)
        dir_out = os.path.join(AOD_MEAN_DIR, 'AOD_MEAN_MODIS_5KM', frequency.upper())
        file_out = os.path.join(dir_out, filename_out)
        # if os.path.isfile(file_out):
        #     print('already exist {}'.format(file_out))
        #     continue

        print('<<< {}'.format(d_))
        for file_ in files:
            print('<<< {}'.format(file_))

            modis_10km = AodModis(file_)
            aod = modis_10km.get_aod()
            lons, lats = modis_10km.get_lon_lat()
            print(np.nanmin(aod), np.nanmax(aod), np.nanmean(aod))
            if aod is None:
                print('aod 为 None: {}'.format(file_))
                continue

            # 投影
            print('投影')
            ii, jj = proj.lonslats2ij(lons, lats)
            valid = np.logical_and.reduce((ii >= 0, ii < proj.row,
                                           jj >= 0, jj < proj.col,
                                           aod > 0, aod < 1.5,
                                           ))
            if valid.sum() == 0:
                print('valid.size == 0, continue')
                continue
            print('valid.sum() == {}'.format(valid.sum()))

            ii = ii[valid]
            jj = jj[valid]
            aod = aod[valid]
            aod_sum[ii, jj] += aod
            aod_count[ii, jj] += 1

            print(np.nanmin(aod), np.nanmax(aod), np.nanmean(aod))
            print(np.nanmin(aod_sum), np.nanmax(aod_sum), np.nanmean(aod_sum))
        if aod_sum is not None and aod_count is not None:
            aod_mean = get_mean(aod_sum, aod_count)
        else:
            continue

        # 新的网格的经纬度
        lons_grid, lats_grid = proj.grid_lonslats()
        make_sure_path_exists(dir_out)
        print((aod_mean != -999).sum())
        data_write = {
            'aod_mean': aod_mean,
            # 'aod_sum': aod_sum,
            # 'aod_count': aod_count,
            'lon': lons_grid,
            'lat': lats_grid
        }
        write_hdf5_and_compress(data_write, file_out)
        break


def mean_modis_10km(frequency='monthly'):
    file_dict = defaultdict(list)
    for root, dirs, files in os.walk(AOD_MODIS_5KM_DIR):
        for name in files:
            if name[-3:].lower() != 'tif':
                continue
            date_str = name.split('-')[1][:6]
            if frequency == 'monthly':
                file_dict[date_str[:6]].append(os.path.join(root, name))
            elif frequency == 'seasonly':
                season = get_season(date_str[:6])
                if season is not None:
                    file_dict[season].append(os.path.join(root, name))
            elif frequency == 'all':
                file_dict['all'].append(os.path.join(root, name))
            else:
                continue

    for d_, files in file_dict.items():
        aod_sum = None
        aod_count = None
        filename_out = 'AOD_MODIS_10KM_MEAN_{}_{}.HDF'.format(frequency, d_)
        dir_out = os.path.join(AOD_MEAN_DIR, 'AOD_MEAN_MODIS_10KM', frequency.upper())
        file_out = os.path.join(dir_out, filename_out)
        # if os.path.isfile(file_out):
        #     print('already exist {}'.format(file_out))
        #     continue

        print('<<< {}'.format(d_))
        for file_ in files:
            print('<<< {}'.format(file_))
            loader = Dataset(file_)
            aod = loader.get_data(1)  # 无效值赋值为-999

            if aod is None:
                print('aod 为 None: {}'.format(file_))
                continue

            valid = np.logical_and(aod > 0, aod < 1500)
            aod[valid] = aod[valid] / 1000.
            aod[~valid] = -999
            # print(np.nanmin(aod), np.nanmax(aod), np.nanmean(aod))
            aod_sum, aod_count = get_sum_count(aod_sum, aod_count, aod)
        if aod_sum is not None and aod_count is not None:
            aod_mean = get_mean(aod_sum, aod_count)
        else:
            continue

        make_sure_path_exists(dir_out)
        loader = Dataset(files[0])
        lon, lat = loader.get_lon_lat()
        data_write = {
            'aod_mean': aod_mean,
            'lon': lon,
            'lat': lat
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
            elif frequency == 'all':
                file_dict['all'].append(os.path.join(root, name))
            else:
                continue

    for d_, files in file_dict.items():
        aod_sum = None
        aod_count = None
        filename_out = 'AOD_FY3D_5KM_MEAN_{}_{}.HDF'.format(frequency, d_)
        dir_out = os.path.join(AOD_MEAN_DIR, 'AOD_MEAN_FY3D_5KM', frequency.upper())
        file_out = os.path.join(dir_out, filename_out)
        if os.path.isfile(file_out):
            print('already exist {}'.format(file_out))
            continue

        print('<<< {}'.format(d_))
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


def diff_mersi_modis(frequency='monthly'):
    dir_out = os.path.join(AOD_MEAN_DIR, 'AOD_DIFF_MERSI_MODIS', frequency.upper())
    print("dir_out <<< : {}".format(dir_out))

    dir_out1 = os.path.join(AOD_MEAN_DIR, 'AOD_MEAN_FY3D_5KM', frequency.upper())
    print("dir_out1 <<< : {}".format(dir_out1))
    filenames = os.listdir(dir_out1)
    files1 = [os.path.join(dir_out1, x) for x in filenames]

    dir_out2 = os.path.join(AOD_MEAN_DIR, 'AOD_MEAN_MODIS_10KM', frequency.upper())
    print("dir_out2 <<< : {}".format(dir_out2))
    filenames = os.listdir(dir_out2)
    files2 = [os.path.join(dir_out2, x) for x in filenames]

    for aod_file1 in files1:
        d_1 = os.path.basename(aod_file1).split('_')[5][:6]
        for aod_file2 in files2:
            d_2 = os.path.basename(aod_file2).split('_')[5][:6]
            if d_1 != d_2:
                continue
            else:
                d_ = d_1
            print('aod_file1 <<< : {}'.format(aod_file1))
            print('aod_file2 <<< : {}'.format(aod_file2))

            aod1 = get_hdf5_data(aod_file1, 'aod_mean', 1, 0, [0, 1.5], np.nan)
            lons1 = get_hdf5_data(aod_file1, 'lon', 1, 0, [-180, 180])
            lats1 = get_hdf5_data(aod_file1, 'lat', 1, 0, [-90, 90])

            aod2 = get_hdf5_data(aod_file2, 'aod_mean', 1, 0, [0, 1.5], np.nan)
            lons2 = get_hdf5_data(aod_file2, 'lon', 1, 0, [-180, 180])
            lats2 = get_hdf5_data(aod_file2, 'lat', 1, 0, [-90, 90])

            # 创建投影查找表
            print('创建投影查找表')
            res_degree = meter2degree(10000)  # 分辨率，10km
            projstr = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
            proj = ProjCore(projstr, res_degree, unit="deg", pt_tl=(69.995, 55.995),
                            pt_br=(139.995, 14.995))  # 角点也要放在格点中心位置

            result1 = proj.create_lut(lats=lats1, lons=lons1)
            result2 = proj.create_lut(lats=lats2, lons=lons2)

            print('投影')
            # 创建投影之后的经纬度
            lons_proj, lats_proj = proj.grid_lonslats()

            # 投影，并且获取投影之后的结果
            aod1_proj = np.full_like(lons_proj, np.nan)  # 可以直接使用proj的row，col进行优化
            aod1_proj[result1['prj_i'], result1['prj_j']] = aod1[result1['pre_i'], result1['pre_j']]

            aod2_proj = np.full_like(lons_proj, np.nan)  # 可以直接使用proj的row，col进行优化
            aod2_proj[result2['prj_i'], result2['prj_j']] = aod2[result2['pre_i'], result2['pre_j']]

            aod_dif = np.full_like(lons_proj, -999)
            match_idx = np.logical_and(np.isfinite(aod1_proj), np.isfinite(aod2_proj))
            aod_dif[match_idx] = aod1_proj[match_idx] - aod2_proj[match_idx]

            result = {
                'aod_mean': aod_dif,
                'lon': lons_proj,
                'lat': lats_proj
            }
            filename_out = 'AOD_DIFF_MERSI_MODIS_{}_{}.HDF'.format(frequency, d_)
            file_out = os.path.join(dir_out, filename_out)
            make_sure_path_exists(dir_out)
            write_hdf5_and_compress(result, file_out)


if __name__ == '__main__':
    # mean_fy3d_5km('monthly')
    # mean_fy3d_5km('seasonly')
    # mean_fy3d_5km('all')
    # mean_modis_10km('monthly')
    # mean_modis_10km('seasonly')
    # mean_modis_10km('all')
    # diff_mersi_modis('monthly')
    # diff_mersi_modis('seasonly')
    # diff_mersi_modis('all')
    mean_modis_10kmto5km('monthly')
    mean_modis_10kmto5km('seasonly')
    mean_modis_10kmto5km('all')
