#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-07-30 12:36
# @Author  : NingAnMe <ninganme@qq.com>
from datetime import datetime
import os
import h5py
import numpy as np


class AodFy3d:

    def __init__(self, in_file, geo_file=None):
        self.in_file = in_file
        self.geo_file = geo_file
        self.filename = os.path.basename(in_file)
        ymdhm = "".join(os.path.splitext(self.filename)[0].split('_')[4:6])
        self.dt = datetime.strptime(ymdhm, '%Y%m%d%H%M')

    @classmethod
    def get_hdf5_data(cls, hdf5_file, data_name, slope=None, intercept=None, valid_range=None):
        with h5py.File(hdf5_file, 'r') as hdf:
            dataset = hdf.get(data_name)
            if slope is None:
                slope = dataset.attrs['Slope']
            if intercept is None:
                intercept = dataset.attrs['Intercept']
            if valid_range is None:
                valid_range = dataset.attrs['valid_range']
            data = dataset[:].astype(np.float)
            data[np.logical_or(data < valid_range[0], data > valid_range[1])] = -999
            data = data * slope + intercept
            return data

    def get_lon_lat(self):
        lons = self.get_hdf5_data(self.geo_file, 'Longitude', 1, 0, (-180, 180))
        lats = self.get_hdf5_data(self.geo_file, 'Latitude', 1, 0, (-90, 90))
        return lons, lats

    def get_aod(self):
        return self.get_hdf5_data(self.in_file, 'Optical_Depth_Land_And_Ocean', 1, 0, (0, 1000))
