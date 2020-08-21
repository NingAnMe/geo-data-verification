#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : NingAnMe <ninganme@qq.com>
from datetime import datetime
import os
import pickle

import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from lib.aeronet import Aeronet

data_path = '/DATA/PROJECT/SourceData/Aeronet/AOD/AOD20/ALL_POINTS'
all_site_pkl = 'all_site.pkl'
filter_site_pkl = 'filter_site.pkl'


def aeronet_site_info():
    files = os.listdir(data_path)
    files.sort()

    site_info = {}
    site_info_filter = {}
    for f in files:
        file_ = os.path.join(data_path, f)
        aeronet = Aeronet(file_)
        name, lon, lat = aeronet.get_site_name_lon_lat()

        dates = aeronet.get_datetime()
        datetime_start = datetime_end = dates[0]
        datetime_max = datetime(2020, 5, 1)
        for d in dates:
            if d < datetime_start:
                datetime_start = d
            if d > datetime_end:
                datetime_end = d
        print("{:20} {} {}".format(name, datetime_start, datetime_end))

        site_info[name] = (lon, lat, datetime_start, datetime_end)
        if datetime_end > datetime_max:
            print(datetime_end)
            site_info_filter[name] = (lon, lat, datetime_start, datetime_end)

    print(f'site_info: {len(site_info)}')
    if not os.path.isfile(all_site_pkl):
        with open(all_site_pkl, 'wb') as fp:
            pickle.dump(site_info, fp)

    print(f'site_info_filter: {len(site_info_filter)}')
    if not os.path.isfile(filter_site_pkl):
        with open(filter_site_pkl, 'wb') as fp:
            pickle.dump(site_info_filter, fp)


def station_datetime_start_end():
    site_info_out = {
        "name": list(),
        "lon": list(),
        "lat": list(),
        "dt_s": list(),
        "dt_e": list(),
    }
    with open(all_site_pkl, 'rb') as fp:
        site_info = pickle.load(fp)
        for k, v in site_info.items():
            site_info_out["name"].append(k)
            site_info_out["lon"].append(v[0])
            site_info_out["lat"].append(v[1])
            site_info_out["dt_s"].append(v[2])
            site_info_out["dt_e"].append(v[3])
        site_info_out = pd.DataFrame(site_info_out)
        site_info_out.to_csv("site_info_lev15.csv")


def __get_site_info():
    with open(filter_site_pkl, 'rb') as fp:
        return pickle.load(fp)


def plot_site_map():
    site_info = __get_site_info()
    print(len(site_info))
    print(site_info)
    lons = list()
    lats = list()
    for lon, lat, _, _ in site_info.values():
        lons.append(lon)
        lats.append(lat)

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()

    plt.scatter(lons, lats,
                color='blue', marker='.',
                transform=ccrs.PlateCarree(),
                )
    plt.savefig('site_map.png')


if __name__ == '__main__':
    # aeronet_site_info()
    plot_site_map()
    station_datetime_start_end()
