#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-09-22 17:05
# @Author  : NingAnMe <ninganme@qq.com>
import os


def make_sure_path_exists(path):
    if not os.path.isdir(path):
        os.makedirs(path)
        print(f'创建文件夹：{path}')