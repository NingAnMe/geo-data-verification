#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-09-21 15:38
# @Author  : NingAnMe <ninganme@qq.com>
import os

LONGITUDE_RANGE_China = [70, 140]
LATITUDE_RANGE_China = [15, 55]

LONGITUDE_RANGE_ChangSanJiao = [114, 123]
LATITUDE_RANGE_ChangSanJiao = [27, 36]

LONGITUDE_RANGE_JingJinJi = [113, 120]
LATITUDE_RANGE_JingJinJi = [36, 43]

LONGITUDE_RANGE_ZhuSanJiao = [110, 116]
LATITUDE_RANGE_ZhuSanJiao = [20, 25]

LONGITUDE_RANGE_FenWei = [105, 115]
LATITUDE_RANGE_FenWei = [32, 40]


# PROJ_DATA = '/DISK/DATA02/PROJECT/SourceData/ShangHai/AOD'
PROJ_DATA = '/home/aodo3/FY3D_AEROSOL_DATA'  # 数据文件根目录


# 结果
AOD_MEAN_DIR = os.path.join(PROJ_DATA, 'MEAN')  # 月、季度平均数据文件夹
AOD_PICTURE_DIR = os.path.join(PROJ_DATA, 'PICTURE')  # 分布图文件夹

# 原数据
AOD_FY3D_1KM_DIR = os.path.join(PROJ_DATA, '2FY3D_MERSI_L2_AOD-ORBT-1000M')  # FY3D数据的文件夹路径
GEO_FY3D_1KM_DIR = os.path.join(PROJ_DATA, '2FY3D-MERSI-L1-GEO1K')  # FY3D数据的文件夹路径
AOD_FY3D_5KM_DIR = os.path.join(PROJ_DATA, '1FY3D-MERSI-L2-AOD-DAILY-5000M')  # FY3D数据的文件夹路径

# 原数据
AOD_MODIS_3KM_DIR = os.path.join(PROJ_DATA, '3MYD04_L2_MODIS_AOD_10KM')  # MODIS数据的文件夹路径
AOD_MODIS_10KM_DIR = os.path.join(PROJ_DATA, '3MYD04_L2_MODIS_AOD_10KM')  # MODIS数据的文件夹路径

# 原数据
AOD_FY4A_4KM_DIR = os.path.join(PROJ_DATA, '4FY4A_AGRI_L2_AOD_4000M')  # FY4A数据的文件夹路径

# 中间文件
AOD_FY3D_1KM_MODIS_3KM_DIR = os.path.join(PROJ_DATA, 'MATCH_FY3D_1KM_MODIS_3KM')
AOD_FY3D_5KM_MODIS_10KM_DIR = os.path.join(PROJ_DATA, 'MATCH_FY3D_5KM_MODIS_10KM')
AOD_FY3D_1KM_FY4A_4KM_DIR = os.path.join(PROJ_DATA, 'MATCH_FY3D_1KM_FY4A_4KM')
AOD_FY3D_5KM_FY4A_4KM_DIR = os.path.join(PROJ_DATA, 'MATCH_FY3D_5KM_FY4A_4KM')

AOD_MODIS_5KM_DIR = os.path.join(PROJ_DATA, '3MYD04_L2_MODIS_AOD_5KM')  # MODIS数据的文件夹路径
AOD_MODIS_100KM_DIR = os.path.join(PROJ_DATA, '5MYD08_L3_MODIS_AOD_100KM')  # MODIS数据的文件夹路径
