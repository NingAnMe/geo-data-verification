#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-07-30 12:36
# @Author  : NingAnMe <ninganme@qq.com>
import os
import h5py
from pyhdf.SD import SD, SDC
import numpy as np


class AodFy3d1km:

    def __init__(self, in_file, geo_file=None):
        self.in_file = in_file
        self.geo_file = geo_file
        self.filename = os.path.basename(in_file)
        # ymdhm = "".join(os.path.splitext(self.filename)[0].split('_')[4:6])
        # self.dt = datetime.strptime(ymdhm, '%Y%m%d%H%M')

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
            valid_index = np.logical_or(data < valid_range[0], data > valid_range[1])
            data = data * slope + intercept
            data[valid_index] = -999
            return data

    def get_lon_lat(self):
        lons = self.get_hdf5_data(self.geo_file, '/Geolocation/Longitude')
        lats = self.get_hdf5_data(self.geo_file, '/Geolocation/Latitude')
        return lons, lats

    def get_aod(self):
        return self.get_hdf5_data(self.in_file, 'AOT_Land', 0.001, 0, (0, 1500))


class AodFy3d5km:

    def __init__(self, in_file, geo_file=None):
        self.in_file = in_file
        self.geo_file = geo_file
        self.filename = os.path.basename(in_file)
        # ymdhm = "".join(os.path.splitext(self.filename)[0].split('_')[4:6])
        # self.dt = datetime.strptime(ymdhm, '%Y%m%d%H%M')

    @classmethod
    def get_hdf5_data(cls, hdf5_file, data_name, slope=None, intercept=None, valid_range=None):
        with h5py.File(hdf5_file, 'r') as hdf:
            dataset = hdf.get(data_name)
            if dataset is not None:
                if slope is None:
                    slope = dataset.attrs['Slope']
                if intercept is None:
                    intercept = dataset.attrs['Intercept']
                if valid_range is None:
                    valid_range = dataset.attrs['valid_range']
                data = dataset[:].astype(np.float)
                valid_index = np.logical_or(data < valid_range[0], data > valid_range[1])
                data = data * slope + intercept
                data[valid_index] = -999
                return data

    @staticmethod
    def get_lon_lat():
        laty = np.arange(90, -90, -0.05)
        lonx = np.arange(-180, 180, 0.05)
        lons, lats = np.meshgrid(lonx, laty)
        return lons, lats

    def get_aod(self):
        return self.get_hdf5_data(self.in_file, 'AOT_Land_Mean', 0.001, 0, (0, 1500))


class AodFy4a4km:

    def __init__(self, in_file, geo_file='lib/FY4A_4KM_GEO.NC'):
        self.in_file = in_file
        self.geo_file = geo_file
        self.filename = os.path.basename(in_file)
        # ymdhm = "".join(os.path.splitext(self.filename)[0].split('_')[4:6])
        # self.dt = datetime.strptime(ymdhm, '%Y%m%d%H%M')

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
            invalid_index = np.logical_or(data < valid_range[0], data > valid_range[1])
            data = data * slope + intercept
            data[invalid_index] = -999
            return data

    def get_lon_lat(self):
        laty = self.get_hdf5_data(self.geo_file, 'lat', 1, 0)
        lonx = self.get_hdf5_data(self.geo_file, 'lon', 1, 0)
        lons, lats = np.meshgrid(lonx, laty)
        return lons, lats

    def get_aod(self):
        return self.get_hdf5_data(self.in_file, 'AOD', 1, 0, (0, 1.5))[1]


class AodModis:
    def __init__(self, in_file, geo_file=None):
        self.in_file = in_file
        self.geo_file = geo_file
        self.filename = os.path.basename(in_file)

    @classmethod
    def get_hdf4_data(cls, hdf4_file, data_name, slope=None, intercept=None, valid_range=None):
        hdf = SD(hdf4_file, SDC.READ)
        dataset = hdf.select(data_name)
        if dataset is not None:
            attrs = dataset.attributes()
            if slope is None:
                slope = attrs['scale_factor']
            if intercept is None:
                intercept = attrs['add_offset']
            if valid_range is None:
                valid_range = attrs['valid_range']
            data = dataset.get().astype(np.float)
            valid_index = np.logical_or(data < valid_range[0], data > valid_range[1])
            data = data * slope + intercept
            data[valid_index] = -999
            return data

    def get_lon_lat(self):
        lons = self.get_hdf4_data(self.in_file, 'Longitude')
        lats = self.get_hdf4_data(self.in_file, 'Latitude')
        return lons, lats

    def get_aod(self):
        return self.get_hdf4_data(self.in_file, 'Optical_Depth_Land_And_Ocean', valid_range=(0, 1500))


class AodModis100km:
    def __init__(self, in_file, geo_file=None):
        self.in_file = in_file
        self.geo_file = geo_file
        self.filename = os.path.basename(in_file)

    @classmethod
    def get_hdf4_data(cls, hdf4_file, data_name, slope=None, intercept=None, valid_range=None):
        hdf = SD(hdf4_file, SDC.READ)
        dataset = hdf.select(data_name)
        if dataset is not None:
            attrs = dataset.attributes()
            if slope is None:
                slope = attrs['scale_factor']
            if intercept is None:
                intercept = attrs['add_offset']
            if valid_range is None:
                valid_range = attrs['valid_range']
            data = dataset.get().astype(np.float)
            valid_index = np.logical_or(data < valid_range[0], data > valid_range[1])
            data = data * slope + intercept
            data[valid_index] = -999
            return data

    @staticmethod
    def get_lon_lat():
        laty = np.arange(90, -90, -1)
        lonx = np.arange(-180, 180, 1)
        lons, lats = np.meshgrid(lonx, laty)
        return lons, lats

    @staticmethod
    def get_lon_lat_5km():
        laty = np.arange(90, -90, -0.05)
        lonx = np.arange(-180, 180, 0.05)
        lons, lats = np.meshgrid(lonx, laty)
        return lons, lats

    def get_aod(self):
        return self.get_hdf4_data(self.in_file, 'Aerosol_Optical_Depth_Land_Mean_Mean', valid_range=(0, 1500))
