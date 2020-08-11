#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Time    : 2020-07-24 12:55
# @Author  : NingAnMe <ninganme@qq.com>
import os
import matplotlib.pyplot as plt

from lib.plot_core import get_ds_font, PlotAx
# from lib import dv_map

ORG_NAME = ''
FONT0 = get_ds_font("OpenSans-Regular.ttf")
FONT_MONO = get_ds_font("DroidSansMono.ttf")

RED = '#f63240'
BLUE = '#1c56fb'
GRAY = '#c0c0c0'
EDGE_GRAY = '#303030'

LINE_WIDTH = 0.5

TICKER_FONT = FONT0.copy()
TICKER_FONT.set_size(11)

TITLE_FONT = FONT0.copy()
TITLE_FONT.set_size(14)

LABEL_FONT = FONT0.copy()
LABEL_FONT.set_size(12)

BOTTOM_FONT = FONT0.copy()
BOTTOM_FONT.set_size(13)

REGRESSION_ANNOTATE_SIZE = 13
REGRESSION_ANNOTATE_COLOR = 'red'


def make_sure_path_exists(path):
    if not os.path.isdir(path):
        os.makedirs(path)
        print(f'创建文件夹：{path}')


def plot_regression(
        x, y, w=None,
        out_file=None,
        title=None,
        x_label=None,
        y_label=None,
        x_range=None,
        y_range=None,
        x_interval=None,
        y_interval=None,
        annotate=None,
        ymd_start=None,
        ymd_end=None,
        ymd=None,
        density=False, ):
    # style_file = os.path.join('plot_regression.mplstyle')
    # plt.style.use(style_file)
    figsize = (5, 5)
    dpi = 100
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax1 = plt.subplot2grid((1, 1), (0, 0))

    plot_ax = PlotAx()

    # ##### 画散点
    marker_size = 5
    if density:
        marker = 'o'
        alpha = 0.8
        zorder = 90
        plot_ax.plot_density_scatter(ax1, x, y, marker=marker, alpha=alpha, marker_size=marker_size,
                                     zorder=zorder)
    else:
        alpha = 0.8  # 透明度
        marker = "o"  # 形状
        color = "b"  # 颜色
        ax1.scatter(x, y, s=marker_size, marker=marker, c=color, lw=0, alpha=alpha, zorder=90)

    # ##### 画回归线
    color = 'r'
    linewidth = 1.2
    zorder = 80
    plot_ax.plot_regression_line(ax1, x, y, w, x_range, color=color, linewidth=linewidth,
                                 zorder=zorder)

    # ##### 画对角线
    color = '#808080'
    linewidth = 1.2
    zorder = 70
    plot_ax.plot_diagonal_line(ax1, x, y, x_range, y_range, color, linewidth, zorder=zorder)

    # ##### 格式化图片
    format_kwargs = {}

    if x_range is not None:
        format_kwargs['x_axis_min'] = x_range[0]
        format_kwargs['x_axis_max'] = x_range[1]
    if y_range is not None:
        format_kwargs['y_axis_min'] = y_range[0]
        format_kwargs['y_axis_max'] = y_range[1]
    if x_label is not None:
        format_kwargs['x_label'] = x_label
    if y_label is not None:
        format_kwargs['y_label'] = y_label
    if x_interval is not None:
        x_major_count = (x_range[1] - x_range[0]) / x_interval + 1
        format_kwargs['x_major_count'] = x_major_count
        if x_major_count <= 11:
            x_minor_count = 4
        else:
            x_minor_count = 1
        format_kwargs['x_minor_count'] = x_minor_count
    if y_interval is not None:
        y_major_count = (y_range[1] - y_range[0]) / y_interval + 1
        format_kwargs['y_major_count'] = y_major_count
        if y_major_count <= 11:
            y_minor_count = 4
        else:
            y_minor_count = 1
        format_kwargs['y_minor_count'] = y_minor_count
    if annotate is not None:
        format_kwargs['annotate'] = annotate
        format_kwargs['annotate_color'] = REGRESSION_ANNOTATE_COLOR
        format_kwargs['annotate_font_size'] = REGRESSION_ANNOTATE_SIZE

    plot_ax.format_ax(ax1, **format_kwargs)

    # ##### 标题 底部文字 LOGO
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.15, top=0.85)

    if '\n' in title:
        title_y = 0.96
    else:
        title_y = 0.94
    fig.suptitle(title, y=title_y, ha='center', fontproperties=TITLE_FONT, fontsize=10)

    bottom_text = ''
    bottom_text_l = 0.7
    bottom_text_b = 0.02
    if ymd_start and ymd_end:
        bottom_text = bottom_text + '%s-%s' % (ymd_start, ymd_end)
    elif ymd:
        bottom_text = bottom_text + '%s' % ymd
    if ORG_NAME is not None:
        bottom_text = bottom_text + '   ' + ORG_NAME
    if bottom_text:
        fig.text(bottom_text_l, bottom_text_b, bottom_text, fontproperties=BOTTOM_FONT)

    # # ##### LOGO
    # plot_fig = PlotFigure()
    # img_width = 0.1
    # img_high = float(figsize[0]) / float(figsize[1]) * img_width
    # img_size = (img_width, img_high)
    # if LOGO_LEFT:
    #     plot_fig.add_image(fig, img_size, LOGO_LEFT, position='LB')
    # if LOGO_RIGHT:
    #     plot_fig.add_image(fig, img_size, LOGO_RIGHT, position='RB')

    # ##### 输出图片
    make_sure_path_exists(os.path.dirname(out_file))
    fig.savefig(out_file, dpi=dpi)
    fig.clear()
    plt.close()
    print('>>> {}'.format(out_file))


# def plot_map_project(
#         latitude,
#         longitude,
#         value,
#         out_file,
#         box=None,
#         title=None,
#         ptype=None,
#         vmin=None,
#         vmax=None,
#         marker=None,
#         markersize=None):
#     if box is not None:
#         box = box
#
#     if title is not None:
#         title = title
#     else:
#         title = "Map"
#
#     if vmin is not None:
#         vmin = vmin
#
#     if vmax is not None:
#         vmax = vmax
#
#     if marker is not None:
#         marker = marker
#     else:
#         marker = 's'
#
#     if markersize is not None:
#         markersize = markersize
#     else:
#         markersize = 5
#
#     p = dv_map.dv_map()
#
#     p.easyplot(latitude, longitude, value, vmin=vmin, vmax=vmax, box=box,
#                ptype=ptype, markersize=markersize, marker=marker)
#
#     p.title = title
#     p.savefig(out_file)
#     print('>>> {}'.format(out_file))
