#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-07-30 16:51
# @Author  : NingAnMe <ninganme@qq.com>
from datetime import datetime
import os
import pickle

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from lib.aeronet import Aeronet


data_path = '/DATA/PROJECT/SourceData/Aeronet/AOD/AOD20/ALL_POINTS'

filenames = os.listdir(data_path)
aeronet_file = os.path.join(data_path, filenames[0])
print(aeronet_file)

aeronet = Aeronet(aeronet_file)

dts = aeronet.get_datetime()
dts = pd.to_datetime(dts, unit='s')
dts_delta = dts - datetime(1994, 3, 17)
dts_delta = np.abs(dts_delta)
print(np.argmin(dts_delta))
print(dts_delta[np.argmin(dts_delta)])
print(dts[np.argmin(dts_delta)])
print(type(dts))
for dt in dts.items():
    print(dt)
    index, timestamp = dt
    break
