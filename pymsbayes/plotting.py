#! /usr/bin/env python

import os
import sys
import math
import string

from pymsbayes.utils.stats import (get_freqs, Partition, IntegerPartition,
        ValidationProbabilities, root_mean_square_error)
from pymsbayes.utils.probability import (almost_equal,
        get_probability_from_bayes_factor)
from pymsbayes.utils.functions import frange
from pymsbayes.utils.parsing import DMCSimulationResults, spreadsheet_iter
from pymsbayes.utils import GLOBAL_RNG
from pymsbayes.utils.messaging import get_logger

_LOG = get_logger(__name__)

MATPLOTLIB_AVAILABLE = False
try:
    import matplotlib
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    _LOG.warning('matplotlib could not be imported; '
            'plotting functionality not supported')

class ScatterData(object):
    def __init__(self, x, y,
            marker = 'o',
            markerfacecolor = 'none',
            markeredgecolor = '0.35',
            markeredgewidth = 0.7,
            linestyle = '',
            zorder = 100,
            **kwargs):
        self.x = x
        self.y = y
        self.marker = marker
        self.markerfacecolor = markerfacecolor
        self.markeredgecolor = markeredgecolor
        self.markeredgewidth = markeredgewidth
        self.linestyle = linestyle
        self.zorder = zorder
        self.kwargs = kwargs

    def plot(self, ax):
        l = ax.plot(self.x, self.y)
        plt.setp(l,
                marker = self.marker,
                linestyle = self.linestyle,
                markerfacecolor = self.markerfacecolor,
                markeredgecolor = self.markeredgecolor,
                markeredgewidth = self.markeredgewidth,
                zorder = self.zorder,
                **self.kwargs)
        return l

class HistData(object):
    def __init__(self, x,
            normed = True,
            bins = 10,
            range = None,
            cumulative = False,
            histtype = 'bar',
            align = 'mid',
            orientation = 'vertical',
            rwidth = None,
            log = False,
            color = None,
            edgecolor = '0.5',
            facecolor = '0.5',
            fill = True,
            hatch = None,
            label = None,
            linestyle = None,
            linewidth = None,
            zorder = 10,
            **kwargs):
        self.x = x
        self.normed = normed
        self.bins = bins
        self.range = range
        self.cumulative = cumulative
        self.histtype = histtype
        self.align = align
        self.orientation = orientation
        self.rwidth = rwidth
        self.log = log
        self.color = color
        self.edgecolor = edgecolor
        self.facecolor = facecolor
        self.fill = fill
        self.hatch = hatch
        self.label = label
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.zorder = zorder
        self.kwargs = kwargs

    def plot(self, ax):
        n, bins, patches = ax.hist(self.x,
                normed = self.normed,
                bins = self.bins,
                range = self.range,
                cumulative = self.cumulative,
                histtype = self.histtype,
                align = self.align,
                orientation = self.orientation,
                rwidth = self.rwidth,
                log = self.log,
                color = self.color,
                edgecolor = self.edgecolor,
                facecolor = self.facecolor,
                fill = self.fill,
                hatch = self.hatch,
                label = self.label,
                linestyle = self.linestyle,
                linewidth = self.linewidth,
                zorder = self.zorder,
                **self.kwargs)
        return n, bins, patches

class VerticalLine(object):
    def __init__(self, x,
            ymin = 0,
            ymax = 1,
            color = '0.5',
            label = None,
            linestyle = '--',
            linewidth = 1.0,
            zorder = 0,
            **kwargs):
        self.x = x
        self.ymin = ymin
        self.ymax = ymax
        self.color = color
        self.label = label
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.zorder = zorder
        self.kwargs = kwargs

    def plot(self, ax):
        ax.axvline(x = self.x, 
                ymin = self.ymin,
                ymax = self.ymax, 
                color = self.color,
                label = self.label,
                linestyle = self.linestyle,
                linewidth = self.linewidth,
                **self.kwargs)

class HorizontalLine(object):
    def __init__(self, y,
            xmin = 0,
            xmax = 1,
            color = '0.5',
            label = None,
            linestyle = '--',
            linewidth = 1.0,
            zorder = 0,
            **kwargs):
        self.y = y
        self.xmin = xmin
        self.xmax = xmax
        self.color = color
        self.label = label
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.zorder = zorder
        self.kwargs = kwargs

    def plot(self, ax):
        ax.axhline(y = self.y, 
                xmin = self.xmin,
                xmax = self.xmax, 
                color = self.color,
                label = self.label,
                linestyle = self.linestyle,
                linewidth = self.linewidth,
                **self.kwargs)

class Ticks(object):
    def __init__(self,
            ticks,
            minor = False,
            labels = None,
            **kwargs):
        self.ticks = ticks
        self.labels = labels
        self.minor = minor
        self.kwargs = kwargs

    def set_x(self, ax):
        ax.set_xticks(ticks = self.ticks, minor = self.minor)
        if self.labels:
            ax.set_xticklabels(labels = self.labels,
                    **self.kwargs)

    def set_y(self, ax):
        ax.set_yticks(ticks = self.ticks, minor = self.minor)
        if self.labels:
            ax.set_yticklabels(labels = self.labels,
                    **self.kwargs)

class TextObj(object):
    def __init__(self, x, y, s,
            rotation = 'horizontal',
            horizontalalignment = 'center',
            verticalalignment = 'bottom',
            weight = 'normal',
            style = 'normal',
            size = 14,
            **kwargs):
        self.x = x
        self.y = y
        self.s = s
        self.rotation = rotation
        self.horizontalalignment = horizontalalignment
        self.verticalalignment = verticalalignment
        self.weight = weight
        self.style = style
        self.size = size
        self.kwargs = kwargs

    def add_to_figure(self, fig):
        fig.text(x = self.x, y = self.y, s = self.s,
                horizontalalignment = self.horizontalalignment,
                verticalalignment = self.verticalalignment,
                rotation = self.rotation,
                weight = self.weight,
                style = self.style,
                size = self.size,
                **self.kwargs)

class ScatterPlot(object):
    def __init__(self, scatter_data_list = [],
            hist_data_list = [],
            vertical_lines = [],
            horizontal_lines = [],
            plot_label = None,
            x_label = None,
            y_label = None,
            left_text = None,
            center_text = None,
            right_text = None,
            position = (1,1,1),
            xlim = (None, None),
            ylim = (None, None),
            xticks_obj = None,
            yticks_obj = None,
            identity_line = False,
            tab = 0.06):
        self.scatter_data_list = list(scatter_data_list)
        self.hist_data_list = list(hist_data_list)
        self.vertical_lines = list(vertical_lines)
        self.horizontal_lines = list(horizontal_lines)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(*position)
        self._plot_label = plot_label
        self._x_label = x_label
        self._y_label = y_label
        self._left_text = left_text
        self._center_text = center_text
        self._right_text = right_text
        self.text_objects = {'plot_label': None,
                             'left': None,
                             'center': None,
                             'right': None,
                             'x_label': None,
                             'y_label': None}
        self.plot_label_size = 14.0
        self.plot_label_weight = 'bold'
        self.plot_label_style = 'normal'
        self.left_text_size = 14.0
        self.left_text_weight = 'normal'
        self.left_text_style = 'normal'
        self.center_text_size = 14.0
        self.center_text_weight = 'normal'
        self.center_text_style = 'normal'
        self.right_text_size = 14.0
        self.right_text_weight = 'normal'
        self.right_text_style = 'normal'
        self.shared_x_ax = None
        self.shared_y_ax = None
        self.xlim_left = xlim[0]
        self.xlim_right = xlim[1]
        self.ylim_bottom = ylim[0]
        self.ylim_top = ylim[1]
        self.xticks_obj = xticks_obj
        self.yticks_obj = yticks_obj
        self.identity_line = identity_line
        self.identity_color = '0.5'
        self.identity_style = '-'
        self.identity_width = 1.0
        self.tab = tab
        self._plot()
        self._reset_text_objects()

    def clear(self):
        self.ax.clear()

    def reset_plot(self):
        self._new_instance()
        self._plot()
        self._reset_text_objects()

    def _new_instance(self):
        geo = self.get_geometry()
        self.ax = self.fig.add_subplot(*geo, sharex = self.shared_x_ax,
                sharey = self.shared_y_ax)

    def _plot(self):
        # using `plot` method rather than `scatter`, because the `collections`
        # attribute created by `scatter` seems difficult to adjust after the
        # initial creation. Whereas `lines` are easy to manipulate.
        for d in self.scatter_data_list:
            self._plot_scatter_data(d)
        for h in self.hist_data_list:
            self._plot_hist_data(h)
        for v in self.vertical_lines:
            self._plot_v_line(v)
        for h in self.horizontal_lines:
            self._plot_h_line(h)
        self.ax.set_xlim(left = self.xlim_left, right = self.xlim_right)
        self.ax.set_ylim(bottom = self.ylim_bottom, top = self.ylim_top)
        if self.xticks_obj:
            self.xticks_obj.set_x(self.ax)
        if self.yticks_obj:
            self.yticks_obj.set_y(self.ax)
        if self.identity_line:
            mn = self.get_minimum()
            mx = self.get_maximum()
            l = self.ax.plot([mn, mx], [mn, mx]) 
            plt.setp(l,
                    color = self.identity_color,
                    linestyle = self.identity_style,
                    linewidth = self.identity_width,
                    marker = '',
                    zorder = 0)
        self.ax.set_xlim(left = self.xlim_left, right = self.xlim_right)
        self.ax.set_ylim(bottom = self.ylim_bottom, top = self.ylim_top)

    def _plot_scatter_data(self, d):
        l = d.plot(self.ax)

    def _plot_hist_data(self, h):
        n, bins, patches = h.plot(self.ax)

    def append_plot(self):
        self._plot()
        self.adjust_text_objects()

    def get_top_text_baseline(self):
        ymin, ymax = self.ax.get_ylim()
        return ymax + (math.fabs(ymax - ymin) * 0.01)

    def get_tab_indent(self):
        xmin, xmax = self.ax.get_xlim()
        return xmin + (math.fabs(xmax - xmin) * self.tab)

    def get_x_center(self):
        xmin, xmax = self.ax.get_xlim()
        return xmin + (math.fabs(xmax - xmin) * 0.5)

    def get_origin(self):
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        return xmin, ymin

    def get_minimum(self):
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        return min([xmin, ymin])

    def get_maximum(self):
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        return max([xmax, ymax])

    def remove_text_object(self, text_obj):
        if text_obj and (text_obj in self.ax.texts):
            self.ax.texts.remove(text_obj)

    def _reset_text_objects(self):
        self.set_plot_label()
        self.set_left_text()
        self.set_center_text()
        self.set_right_text()
        self.set_xlabel()
        self.set_ylabel()

    def set_plot_label(self, label = None, fontdict = None, withdash = False,
            **kwargs):
        if label:
            self._plot_label = label
        if self._plot_label is None:
            return
        self.remove_text_object(self.text_objects['plot_label'])
        verticalalignment = 'bottom'
        horizontalalignment = 'left'
        self.text_objects['plot_label'] = self.ax.text(x = 0, y = 0,
                s = self._plot_label,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)
        self._adjust_plot_label()

    def set_left_text(self, left_text = None, fontdict = None, withdash = False,
            **kwargs):
        if left_text:
            self._left_text = left_text
        if self._left_text is None:
            return
        self.remove_text_object(self.text_objects['left'])
        verticalalignment = 'bottom'
        horizontalalignment = 'left'
        self.text_objects['left'] = self.ax.text(x = 0, y = 0,
                s = self._left_text,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)
        self._adjust_left_text()

    def set_center_text(self, center_text = None, fontdict = None, withdash = False,
            **kwargs):
        if center_text:
            self._center_text = center_text
        if self._center_text is None:
            return
        self.remove_text_object(self.text_objects['center'])
        verticalalignment = 'bottom'
        horizontalalignment = 'center'
        self.text_objects['center'] = self.ax.text(x = 0, y = 0,
                s = self._center_text,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)
        self._adjust_center_text()

    def set_right_text(self, right_text = None, fontdict = None, withdash = False,
            **kwargs):
        if right_text:
            self._right_text = right_text
        if self._right_text is None:
            return
        self.remove_text_object(self.text_objects['right'])
        verticalalignment = 'bottom'
        horizontalalignment = 'right'
        self.text_objects['right'] = self.ax.text(x = 0, y = 0,
                s = self._right_text,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)
        self._adjust_right_text()

    def adjust_text_objects(self):
        self._adjust_plot_label()
        self._adjust_left_text()
        self._adjust_center_text()
        self._adjust_right_text()

    def _adjust_plot_label(self):
        if not self.text_objects.get('plot_label', None):
            return
        x = self.get_xmin()
        y = self.get_top_text_baseline()
        self.text_objects['plot_label'].set_position((x, y))
        self.text_objects['plot_label'].set_size(self.plot_label_size)
        self.text_objects['plot_label'].set_style(self.plot_label_style)
        self.text_objects['plot_label'].set_weight(self.plot_label_weight)

    def _adjust_left_text(self):
        if not self.text_objects.get('left', None):
            return
        x = self.get_tab_indent()
        y = self.get_top_text_baseline()
        self.text_objects['left'].set_position((x, y))
        self.text_objects['left'].set_size(self.left_text_size)
        self.text_objects['left'].set_style(self.left_text_style)
        self.text_objects['left'].set_weight(self.left_text_weight)

    def _adjust_center_text(self):
        if not self.text_objects.get('center', None):
            return
        x = self.get_x_center()
        y = self.get_top_text_baseline()
        self.text_objects['center'].set_position((x, y))
        self.text_objects['center'].set_size(self.center_text_size)
        self.text_objects['center'].set_style(self.center_text_style)
        self.text_objects['center'].set_weight(self.center_text_weight)

    def _adjust_right_text(self):
        if not self.text_objects.get('right', None):
            return
        x = self.get_xmax()
        y = self.get_top_text_baseline()
        self.text_objects['right'].set_position((x, y))
        self.text_objects['right'].set_size(self.right_text_size)
        self.text_objects['right'].set_style(self.right_text_style)
        self.text_objects['right'].set_weight(self.right_text_weight)

    def set_xlabel(self, xlabel = None, fontdict = None, labelpad = None,
            **kwargs):
        if xlabel:
            self._x_label = xlabel
        if self._x_label is None:
            return
        self.remove_text_object(self.text_objects['x_label'])
        self.text_objects['x_label'] = self.ax.set_xlabel(
                xlabel = self._x_label,
                fontdict = fontdict,
                labelpad = labelpad,
                **kwargs)

    def set_ylabel(self, ylabel = None, fontdict = None, labelpad = None,
            **kwargs):
        if ylabel:
            self._y_label = ylabel
        if self._y_label is None:
            return
        self.remove_text_object(self.text_objects['y_label'])
        self.text_objects['y_label'] = self.ax.set_ylabel(
                ylabel = self._y_label,
                fontdict = fontdict,
                labelpad = labelpad,
                multialignment = 'center',
                **kwargs)

    def set_xlim(self, left = None, right = None, emit = True, auto = False,
            **kwargs):
        self.xlim_left = left
        self.xlim_right = right
        self.ax.set_xlim(left = left, right = right, emit = emit, auto = auto,
                **kwargs)

    def set_ylim(self, bottom = None, top = None, emit = True, auto = False,
            **kwargs):
        self.ylim_bottom = bottom
        self.ylim_top = top
        self.ax.set_ylim(bottom = bottom, top = top, emit = emit, auto = auto,
                **kwargs)

    def get_xlim(self):
        return self.ax.get_xlim()

    def get_xmin(self):
        return self.get_xlim()[0]

    def get_xmax(self):
        return self.get_xlim()[-1]

    def get_ylim(self):
        return self.ax.get_ylim()

    def get_ymin(self):
        return self.get_ylim()[0]

    def get_ymax(self):
        return self.get_ylim()[-1]

    def get_geometry(self):
        return self.ax.get_geometry()

    def change_geometry(self, numrows, numcols, num):
        self.ax.change_geometry(numrows = numrows, numcols = numcols, num = num)
        self.reset_plot()

    def set_figure(self, fig):
        self.fig = fig
        self.ax.set_figure(fig)

    def get_figure(self):
        return self.ax.get_figure()

    def is_last_row(self):
        return self.ax.is_last_row()

    def is_first_col(self):
        return self.ax.is_first_col()

    def _plot_v_line(self, v):
        v.plot(self.ax)

    def add_v_line(self, vertical_line_object):
        self.vertical_lines.append(vertical_line_object)
        self._plot_v_line(vertical_line_object)

    def _plot_h_line(self, h):
        h.plot(self.ax)

    def add_h_line(self, horizontal_line_object):
        self.horizontal_lines.append(horizontal_line_object)
        self._plot_h_line(horizontal_line_object)

    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)


class PlotGrid(object):
    valid_label_schemas = ['uppercase', 'lowercase', 'numbers']

    def __init__(self, subplots = [],
            num_columns = 2,
            share_x = False,
            share_y = False,
            label_schema = 'uppercase',
            title = None,
            title_top = True,
            y_title = None,
            y_title_position = 0.001,
            height = 6.0,
            width = 8.0,
            auto_height = True):
        self.num_columns = num_columns
        self._set_label_schema(label_schema)
        self.share_x = share_x
        self.share_y = share_y
        self.plot_label_size = 14.0
        self.plot_label_weight = 'bold'
        self.plot_label_style = 'normal'
        self.plot_label_suffix = ''
        self.title = title
        self.title_top = title_top
        self._width = width
        self._height = height
        self.fig = plt.figure(figsize = self.size)
        self.subplots = subplots
        self.auto_height = auto_height
        self.y_title = None
        if y_title:
            self.y_title = TextObj(x = y_title_position,
                    y = 0.5,
                    s = y_title,
                    rotation = 'vertical',
                    horizontalalignment = 'left',
                    verticalalignment = 'center',
                    size = 14)
        self.perimeter_padding = 0.25
        self.padding_between_vertical = 0.8
        self.padding_between_horizontal = None
        self.margin_left = 0
        self.margin_right = 1
        self.margin_bottom = 0
        self.margin_top = 0.975
        self.auto_adjust_margins = True
        for sp in self.subplots:
            f = sp.get_figure()
            if f != self.fig:
                plt.close(f)
                sp.set_figure(self.fig)
        self.reset_figure()

    def _get_size(self):
        return (self._width, self._height)

    def _set_size(self, width_height_tup):
        self._width, self._height = width_height_tup
        self.fig.set_size_inches(self._get_size())

    size = property(_get_size, _set_size)

    def _get_width(self):
        return self._width

    def _set_width(self, width):
        self._width = width
        self.size = self.size

    width = property(_get_width, _set_width)

    def _get_height(self):
        return self._height

    def _set_height(self, height):
        self._height = height
        self.size = self.size

    height = property(_get_height, _set_height)

    def _get_label_schema(self):
        return self._label_schema

    def _set_label_schema(self, schema):
        if schema:
            schema = schema.lower()
        self._label_schema = schema
        if self._label_schema and (
                self._label_schema not in self.valid_label_schemas):
            raise ValueError('invalid label schema {0}; valid options:'
                    '\n\t{1}'.format(self._label_schema,
                            ', '.join(self.valid_label_schemas)))

    label_schema = property(_get_label_schema, _set_label_schema)

    def get_num_rows(self):
        return int(math.ceil(len(self.subplots) / float(self.num_columns)))

    def get_plot_labels(self):
        if not self.label_schema:
            return None
        l = len(self.subplots)
        s = self.plot_label_suffix
        if self.label_schema == 'uppercase':
            return [c + s for c in string.ascii_uppercase[0: l]]
        elif self.label_schema == 'lowercase':
            return [c + s for c in string.ascii_lowercase[0: l]]
        elif self.label_schema == 'numbers':
            return [str(x) + s for x in range(1, l + 1)]
        else:
            raise Exception('invalid label schema {0!r}'.format(
                    self.label_schema))

    def reset_figure(self):
        nrows = self.get_num_rows()
        ncols = self.num_columns
        if self.auto_height:
            self.height = (self.width / (ncols * 1.1)) * nrows
        self.fig = plt.figure(figsize = (self.width, self.height))
        plot_labels = self.get_plot_labels()
        for i, subplot in enumerate(self.subplots):
            f = subplot.get_figure()
            if f != self.fig:
                plt.close(f)
                subplot.set_figure(self.fig)
            subplot.change_geometry(numrows = nrows, numcols = ncols,
                    num = i + 1)
            subplot.plot_label_size = self.plot_label_size
            subplot.plot_label_weight = self.plot_label_weight
            subplot.plot_label_style = self.plot_label_style
            label = plot_labels[i]
            subplot.set_plot_label(label = label)
            if self.share_x:
                subplot.shared_x_ax = self.subplots[0].ax
            if self.share_y:
                subplot.shared_y_ax = self.subplots[0].ax
            subplot.ax.tick_params(labelbottom = True)
            subplot.ax.tick_params(labelleft = True)
            if self.share_x and (not subplot.is_last_row()):
                subplot.ax.tick_params(labelbottom = False)
            if self.share_y and (not subplot.is_first_col()):
                subplot.ax.tick_params(labelleft = False)
        rect = [0, 0, 1, 0.975]
        if self.title:
            if self.title_top:
                self.fig.suptitle(self.title,
                        verticalalignment = 'top',
                        horizontalalignment = 'center',
                        y = 0.999)
                if self.auto_adjust_margins:
                    self.margin_top -= 0.55
            else:
                self.fig.suptitle(self.title,
                        verticalalignment = 'bottom',
                        horizontalalignment = 'center',
                        y = 0.001)
                if self.auto_adjust_margins:
                    self.margin_bottom += 0.06
        if self.y_title:
            self.y_title.add_to_figure(self.fig)
            if self.auto_adjust_margins:
                self.margin_left += 0.02
        if self.subplots:
            rect = (self.margin_left, self.margin_bottom, self.margin_right,
                    self.margin_top)
            self.fig.tight_layout(pad = self.perimeter_padding,
                    h_pad = self.padding_between_vertical,
                    w_pad = self.padding_between_horizontal,
                    rect = rect) # available space on figure

    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)

class PowerPlotGrid(object):
    valid_variables = ['psi', 'omega', 'tau_exclusion']
    variable_symbols = {'psi': r'\Psi',
                        'omega': r'\Omega',
                        'tau_exclusion': r'\tau'}
    def __init__(self,
            observed_config_to_estimates,
            variable = 'psi',
            variable_symbol = None,
            num_columns = 2,
            x_title = None,
            y_title = 'Density',
            width = 8,
            height = 9,
            auto_height = False,
            auto_adjust_margins = False,
            margin_left = 0.025,
            margin_bottom = 0.025,
            margin_right = 1,
            margin_top = 0.98,
            padding_between_horizontal = 0.5,
            padding_between_vertical = 1.3,
            tab = 0.08):
        self.config_estimates_tups = sorted(
            [(c, list(e)) for c, e in observed_config_to_estimates.iteritems()],
            key = lambda x : x[0].tau.mean)
        self.variable = variable.lower()
        self.variable_symbol = variable_symbol
        if not self.variable_symbol:
            self.variable_symbol = self.variable_symbols[self.variable]
        if not self.variable in self.valid_variables:
            raise ValueError('{0!r} is not a valid variable; valid options:'
                    '\n\t{1}'.format(self.variable,
                            ', '.join(self.valid_variables)))
        self.num_taxon_pairs = self.config_estimates_tups[0][0].npairs
        for c, e in self.config_estimates_tups:
            if c.npairs != self.num_taxon_pairs:
                raise ValueError('configs have differing number of taxa')
        self.num_columns = num_columns
        self.subplots = []
        self.plot_grid = None
        self.x_title = x_title
        if self.x_title is None:
            if self.variable == 'psi' or self.variable == 'omega':
                self.x_title = r'$\hat{{{0}}}$'.format(self.variable_symbol)
            elif self.variable == 'tau_exclusion':
                self.x_title = r'Number of true ${0}$ excluded'.format(self.variable_symbol)
            else:
                raise Exception('unexpected variable {0!r}'.format(
                        self.variable))
        self.y_title = y_title
        self.width = width
        self.height = height
        self.auto_height = auto_height
        self.auto_adjust_margins = auto_adjust_margins
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_bottom = margin_bottom
        self.margin_top = margin_top
        self.padding_between_horizontal = padding_between_horizontal
        self.padding_between_vertical = padding_between_vertical
        self.bins = 20
        if self.variable == 'psi':
            self.bins = range(1, self.num_taxon_pairs + 2)
        elif self.variable == 'tau_exclusion':
            mx = [max(e) for (c, e) in self.config_estimates_tups]
            mx = max(mx)
            self.bins = range(0, mx + 2)
        self.tab = tab
        self.vertical_lines = []
        if self.variable == 'omega':
            self.vertical_lines.append(VerticalLine(
                    x = 0.01,
                    color = '0.25'))
        self.populate_subplots()

    def populate_subplots(self):
        for cfg, estimates in self.config_estimates_tups:
            dist = r'$\tau \sim {0}$'.format(str(cfg.tau))
            prob = None
            if self.variable == 'psi':
                p = estimates.count(1) / float(len(estimates))
                prob = r'$p(\hat{{{0}}} = 1) = {1}$'.format(
                        self.variable_symbol,
                        p)
            elif self.variable == 'omega':
                c = len([e for e in estimates if e < 0.01])
                p = c / float(len(estimates))
                prob = r'$p(\hat{{{0}}} < 0.01) = {1}$'.format(
                        self.variable_symbol,
                        p)
            elif self.variable == 'tau_exclusion':
                matplotlib.rc('mathtext',**{'fontset': 'stix'})
                c = len([e for e in estimates if e > 0])
                p = c /float(len(estimates))
                prob = (r'$p(\mathbf{{\tau}} \, \notin \, \hat{{M}}) '
                        r'= {0}$'.format(p))
            if len([e for e in estimates if ((e > 0.00000001) or (e < -0.00000001))]) < 1:
                self.bins = list(frange(0, 0.1, 20, include_end_point = True))
            hd = HistData(x = estimates,
                    normed = True,
                    bins = self.bins,
                    histtype = 'bar',
                    align = 'mid',
                    orientation = 'vertical',
                    zorder = 0)
            # freqs = get_freqs(estimates)
            # s = ScatterPlot()
            # f, bins, patches = hd.plot(s.ax)
            # for i, v in enumerate(f):
            #     assert almost_equal(v, freqs.get(i + 1, 0))
            xticks_obj = None
            if (self.variable == 'psi') or (self.variable == 'tau_exclusion'):
                tick_labels = []
                for x in self.bins[0:-1]:
                    if x % 2:
                        tick_labels.append(str(x))
                    else:
                        tick_labels.append('')
                xticks_obj = Ticks(ticks = self.bins,
                        labels = tick_labels,
                        horizontalalignment = 'left')
            hist = ScatterPlot(hist_data_list = [hd],
                    vertical_lines = self.vertical_lines,
                    left_text = dist,
                    right_text = prob,
                    xticks_obj = xticks_obj,
                    tab = self.tab)
            if self.variable == 'omega':
                hist.left_text_size = 12.0
                hist.right_text_size = 12.0
                xticks = [i for i in hist.ax.get_xticks()]
                xtick_labels = [i for i in xticks]
                if len(xtick_labels) >= 10:
                    for i in range(1, len(xtick_labels), 2):
                        xtick_labels[i] = ''
                xticks_obj = Ticks(ticks = xticks,
                        labels = xtick_labels,
                        horizontalalignment = 'center')
                hist.xticks_obj = xticks_obj
            if self.variable == 'psi':
                hist.set_xlim(left = (self.bins[0]), right = (self.bins[-1]))
            self.subplots.append(hist)

    def create_grid(self):
        if len(self.subplots) < 2:
            self.num_columns = 1
        share_x = True
        if self.variable == 'omega':
            share_x = False
        self.plot_grid = PlotGrid(subplots = self.subplots,
                num_columns = self.num_columns,
                share_x = share_x,
                share_y = False,
                label_schema = 'uppercase',
                title = self.x_title,
                title_top = False,
                y_title = self.y_title,
                width = self.width,
                height = self.height,
                auto_height = self.auto_height)
        self.plot_grid.auto_adjust_margins = self.auto_adjust_margins
        self.plot_grid.margin_left = self.margin_left
        self.plot_grid.margin_bottom = self.margin_bottom
        self.plot_grid.margin_right = self.margin_right
        self.plot_grid.margin_top = self.margin_top
        self.plot_grid.padding_between_horizontal = \
                self.padding_between_horizontal
        self.plot_grid.padding_between_vertical = self.padding_between_vertical
        self.plot_grid.reset_figure()
        return self.plot_grid

class AccuracyValidationPlotGrid(object):
    def __init__(self,
            validation_result_obj,
            num_subsample = 5000,
            rng = None,
            omega_symbol = r'\Omega',
            psi_symbol = r'\Psi',
            mean_time_symbol = r'E(\tau)',
            math_font = None,
            width = 8,
            height = 5,
            auto_height = False,
            auto_adjust_margins = False,
            margin_left = 0.0,
            margin_bottom = 0.0,
            margin_right = 1,
            margin_top = 0.975,
            padding_between_horizontal = 0.5,
            padding_between_vertical = 1.1,
            tab = 0.08):
        self.omega_symbol = omega_symbol
        self.psi_symbol = psi_symbol
        self.mean_time_symbol = mean_time_symbol
        self.math_font = math_font
        self.width = width
        self.height = height
        self.auto_height = auto_height
        self.auto_adjust_margins = auto_adjust_margins
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_bottom = margin_bottom
        self.margin_top = margin_top
        self.padding_between_horizontal = padding_between_horizontal
        self.padding_between_vertical = padding_between_vertical
        self.tab = tab
        self.num_columns = 3
        if not rng:
            rng = GLOBAL_RNG
        self.v = validation_result_obj.get_random_subsample(num_subsample, rng)
        self.v.psi.jitter_true_and_modes(sigma = 0.02, rng = rng)
        self.subplots = []
        self.plot_grid = None
        self.populate_subplots()

    def populate_subplots(self):
        if self.math_font:
            matplotlib.rc('mathtext',**{'fontset': self.math_font})
        else:
            matplotlib.rc('text',**{'usetex': True})
        self.subplots = []

        x = self.v.psi.true_jitter
        y = self.v.psi.mode_jitter
        rmse = root_mean_square_error(x, y)
        sd_psi = ScatterData(
                x = x,
                y = y)
        sp_psi = ScatterPlot(scatter_data_list = [sd_psi],
                x_label = r'True ${0}$'.format(self.psi_symbol),
                y_label = '\\textit{{\\textbf{{Unadjusted}}}}\n$\\hat{{{0}}}$ (mode)'.format(self.psi_symbol),
                right_text = r'$RMSE = {0:.2f}$'.format(rmse),
                identity_line = True,
                tab = self.tab)
        sp_psi.right_text_size = 10.0

        x = self.v.psi.true_jitter
        y = self.v.psi.mode_glm_jitter
        rmse = root_mean_square_error(x, y)
        sd_psi_glm = ScatterData(
                x = x,
                y = y)
        sp_psi_glm = ScatterPlot(scatter_data_list = [sd_psi_glm],
                x_label = r'True ${0}$'.format(self.psi_symbol),
                y_label = '\\textit{{\\textbf{{GLM-adjusted}}}}\n$\\hat{{{0}}}$ (mode)'.format(self.psi_symbol),
                right_text = r'$RMSE = {0:.2f}$'.format(rmse),
                identity_line = True,
                tab = self.tab)
        sp_psi_glm.right_text_size = 10.0

        x = self.v.omega.true
        y = self.v.omega.median
        rmse = root_mean_square_error(x, y)
        sd_omega = ScatterData(
                x = x,
                y = y)
        sp_omega = ScatterPlot(scatter_data_list = [sd_omega],
                x_label = r'True ${0}$'.format(self.omega_symbol),
                y_label = r'$\hat{{{0}}}$ (median)'.format(self.omega_symbol),
                right_text = r'$RMSE = {0:.2f}$'.format(rmse),
                identity_line = True,
                tab = self.tab)
        sp_omega.right_text_size = 10.0

        x = self.v.omega.true
        y = self.v.omega.mode_glm
        rmse = root_mean_square_error(x, y)
        sd_omega_glm = ScatterData(
                x = x,
                y = y)
        sp_omega_glm = ScatterPlot(scatter_data_list = [sd_omega_glm],
                x_label = r'True ${0}$'.format(self.omega_symbol),
                y_label = r'$\hat{{{0}}}$ (mode)'.format(self.omega_symbol),
                right_text = r'$RMSE = {0:.2f}$'.format(rmse),
                identity_line = True,
                tab = self.tab)
        sp_omega_glm.right_text_size = 10.0

        x = self.v.tau.true
        y = self.v.tau.median
        rmse = root_mean_square_error(x, y)
        sd_tau = ScatterData(
                x = x,
                y = y)
        sp_tau = ScatterPlot(scatter_data_list = [sd_tau],
                x_label = r'True ${0}$'.format(self.mean_time_symbol),
                y_label = r'$\hat{{{0}}}$ (median)'.format(self.mean_time_symbol),
                right_text = r'$RMSE = {0:.2f}$'.format(rmse),
                identity_line = True,
                tab = self.tab)
        sp_tau.right_text_size = 10.0

        x = self.v.tau.true
        y = self.v.tau.mode_glm
        rmse = root_mean_square_error(x, y)
        sd_tau_glm = ScatterData(
                x = x,
                y = y)
        sp_tau_glm = ScatterPlot(scatter_data_list = [sd_tau_glm],
                x_label = r'True ${0}$'.format(self.mean_time_symbol),
                y_label = r'$\hat{{{0}}}$ (mode)'.format(self.mean_time_symbol),
                right_text = r'$RMSE = {0:.2f}$'.format(rmse),
                identity_line = True,
                tab = self.tab)
        sp_tau_glm.right_text_size = 10.0
        self.subplots = [sp_psi, sp_omega, sp_tau, sp_psi_glm, sp_omega_glm,
                sp_tau_glm]
        for sp in self.subplots:
            mx = max(sp.scatter_data_list[0].x + sp.scatter_data_list[0].y)
            mn = min(sp.scatter_data_list[0].x + sp.scatter_data_list[0].y)
            buff = (mx - mn) * 0.04
            lm = (mn - buff, mx + buff)
            sp.set_xlim(left = lm[0], right = lm[1])
            sp.set_ylim(bottom = lm[0], top = lm[1])
            xticks = [i for i in sp.ax.get_xticks()]
            yticks = [i for i in sp.ax.get_yticks()]
            xtick_labels = [i for i in xticks]
            if len(xtick_labels) >= 10:
                for i in range(1, len(xtick_labels), 2):
                    xtick_labels[i] = ''
            xticks_obj = Ticks(ticks = xticks,
                    labels = xtick_labels,
                    horizontalalignment = 'center',
                    size = 10.0)
            yticks_obj = Ticks(ticks = yticks,
                    labels = yticks,
                    size = 10.0)
            sp.xticks_obj = xticks_obj
            sp.yticks_obj = yticks_obj

    def create_grid(self):
        if len(self.subplots) < 2:
            self.num_columns = 1
        self.plot_grid = PlotGrid(subplots = self.subplots,
                num_columns = self.num_columns,
                share_x = False,
                share_y = False,
                label_schema = 'uppercase',
                width = self.width,
                height = self.height,
                auto_height = self.auto_height)
        self.plot_grid.auto_adjust_margins = self.auto_adjust_margins
        self.plot_grid.margin_left = self.margin_left
        self.plot_grid.margin_bottom = self.margin_bottom
        self.plot_grid.margin_right = self.margin_right
        self.plot_grid.margin_top = self.margin_top
        self.plot_grid.padding_between_horizontal = \
                self.padding_between_horizontal
        self.plot_grid.padding_between_vertical = self.padding_between_vertical
        self.plot_grid.reset_figure()
        return self.plot_grid

class ProbabilityValidationPlotGrid(object):
    def __init__(self,
            psi_validation_probs,
            psi_validation_probs_glm,
            omega_validation_probs,
            omega_validation_probs_glm,
            omega_symbol = r'\Omega',
            psi_symbol = r'\Psi',
            math_font = None,
            width = 8,
            height = 6,
            auto_height = False,
            auto_adjust_margins = False,
            margin_left = 0.0,
            margin_bottom = 0.0,
            margin_right = 1,
            margin_top = 0.985,
            padding_between_horizontal = 0.5,
            padding_between_vertical = 1.0,
            tab = 0.08):
        self.omega_symbol = omega_symbol
        self.psi_symbol = psi_symbol
        self.math_font = math_font
        self.width = width
        self.height = height
        self.auto_height = auto_height
        self.auto_adjust_margins = auto_adjust_margins
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_bottom = margin_bottom
        self.margin_top = margin_top
        self.padding_between_horizontal = padding_between_horizontal
        self.padding_between_vertical = padding_between_vertical
        self.tab = tab
        self.num_columns = 2
        self.psi = psi_validation_probs
        self.psi_glm = psi_validation_probs_glm
        self.omega = omega_validation_probs
        self.omega_glm = omega_validation_probs_glm
        self.subplots = []
        self.plot_grid = None
        self.populate_subplots()

    def populate_subplots(self):
        if self.math_font:
            matplotlib.rc('mathtext',**{'fontset': self.math_font})
        else:
            matplotlib.rc('text',**{'usetex': True})
        self.subplots = []
        sd_psi = ScatterData(
                x = self.psi.estimated_probs,
                y = self.psi.true_probs)
        sd_psi_glm = ScatterData(
                x = self.psi_glm.estimated_probs,
                y = self.psi_glm.true_probs)
        sd_omega = ScatterData(
                x = self.omega.estimated_probs,
                y = self.omega.true_probs)
        sd_omega_glm = ScatterData(
                x = self.omega_glm.estimated_probs,
                y = self.omega_glm.true_probs)
        lm = (0.0, 1.0)
        sp_psi = ScatterPlot(scatter_data_list = [sd_psi],
                y_label = '\\textit{\\textbf{Unadjusted}}\nTrue probability',
                xlim = lm,
                ylim = lm,
                identity_line = True,
                tab = self.tab)
        sp_psi_glm = ScatterPlot(scatter_data_list = [sd_psi_glm],
                x_label = r'Estimated $p({0} = 1 \, | \, B_{{\epsilon}}(S*))$'.format(self.psi_symbol),
                y_label = '\\textit{\\textbf{GLM-adjusted}}\nTrue probability',
                xlim = lm,
                ylim = lm,
                identity_line = True,
                tab = self.tab)
        sp_omega = ScatterPlot(scatter_data_list = [sd_omega],
                xlim = lm,
                ylim = lm,
                identity_line = True,
                tab = self.tab)
        sp_omega_glm = ScatterPlot(scatter_data_list = [sd_omega_glm],
                x_label = r'Estimated $p({0} < 0.01 \, | \, B_{{\epsilon}}(S*))$'.format(self.omega_symbol),
                xlim = lm,
                ylim = lm,
                identity_line = True,
                tab = self.tab)
        self.subplots.extend([sp_psi, sp_omega, sp_psi_glm, sp_omega_glm])

    def create_grid(self):
        if len(self.subplots) < 2:
            self.num_columns = 1
        self.plot_grid = PlotGrid(subplots = self.subplots,
                num_columns = self.num_columns,
                share_x = True,
                share_y = True,
                label_schema = 'uppercase',
                width = self.width,
                height = self.height,
                auto_height = self.auto_height)
        self.plot_grid.auto_adjust_margins = self.auto_adjust_margins
        self.plot_grid.margin_left = self.margin_left
        self.plot_grid.margin_bottom = self.margin_bottom
        self.plot_grid.margin_right = self.margin_right
        self.plot_grid.margin_top = self.margin_top
        self.plot_grid.padding_between_horizontal = \
                self.padding_between_horizontal
        self.plot_grid.padding_between_vertical = self.padding_between_vertical
        self.plot_grid.reset_figure()
        return self.plot_grid

class ProbabilityPowerPlotGrid(object):
    valid_variables = ['psi', 'omega', 'tau_exclusion']
    variable_symbols = {'psi': r'\Psi',
                        'omega': r'\Omega',
                        'tau_exclusion': r'\tau'}
    valid_div_model_priors = ['psi', 'dpp', 'uniform']
    def __init__(self,
            observed_config_to_estimates,
            variable = 'psi',
            variable_symbol = None,
            div_model_prior = 'psi',
            bayes_factor = 10,
            bayes_factor_prob = None,
            cfg_to_prob_of_bf_exclusion = None,
            bayes_factor_line_color = '0.25',
            dpp_concentration_mean = None,
            num_columns = 2,
            x_title = None,
            y_title = 'Density',
            width = 8,
            height = 9,
            auto_height = False,
            auto_adjust_margins = False,
            margin_left = 0.025,
            margin_bottom = 0.025,
            margin_right = 1,
            margin_top = 0.98,
            padding_between_horizontal = 0.5,
            padding_between_vertical = 1.0,
            text_size = 12.0,
            tab = 0.08):
        self.config_estimates_tups = sorted(
            [(c, list(e)) for c, e in observed_config_to_estimates.iteritems()],
            key = lambda x : x[0].tau.mean)
        self.variable = variable.lower()
        self.variable_symbol = variable_symbol
        if not self.variable_symbol:
            self.variable_symbol = self.variable_symbols[self.variable]
        if not self.variable in self.valid_variables:
            raise ValueError('{0!r} is not a valid variable; valid options:'
                    '\n\t{1}'.format(self.variable,
                            ', '.join(self.valid_variables)))
        self.num_taxon_pairs = self.config_estimates_tups[0][0].npairs
        for c, e in self.config_estimates_tups:
            if c.npairs != self.num_taxon_pairs:
                raise ValueError('configs have differing number of taxa')
        self.div_model_prior = div_model_prior.lower()
        if not self.div_model_prior in self.valid_div_model_priors:
            raise ValueError('{0!r} is not a valid div model prior; valid '
                    'options:\n\t{1}'.format(self.div_model_prior,
                            ', '.join(self.valid_div_model_priors)))
        self.dpp_concentration_mean = dpp_concentration_mean
        if (self.div_model_prior == 'dpp') and (
                not self.dpp_concentration_mean):
            raise ValueError('if the div model prior is {0!r}, you need '
                    'to specify the `dpp_concentration_mean`.'.format(
                            self.div_model_prior))
        self.one_div_prior = 0.0
        if self.div_model_prior == 'psi':
            ip = IntegerPartition([0] * self.num_taxon_pairs)
            self.one_div_prior = ip.psi_uniform_prior_probability()
        elif self.div_model_prior == 'uniform':
            ip = IntegerPartition([0] * self.num_taxon_pairs)
            self.one_div_prior = ip.uniform_prior_probability()
        elif self.div_model_prior == 'dpp':
            p = Partition([0] * self.num_taxon_pairs)
            self.one_div_prior = p.dirichlet_process_prior_probability(
                    alpha = self.dpp_concentration_mean)
        else:
            raise Exception('unexpected div model prior {0!r}'.format(
                    self.div_model_prior))
        self.bayes_factor = float(bayes_factor)
        self.cfg_to_prob_of_bf_exclusion = cfg_to_prob_of_bf_exclusion
        self.bayes_factor_line_color = str(bayes_factor_line_color)
        self.bayes_factor_prob = bayes_factor_prob
        if not self.bayes_factor_prob:
            if self.variable == 'psi':
                self.bayes_factor_prob = get_probability_from_bayes_factor(
                        self.bayes_factor, self.one_div_prior)
        self.num_columns = num_columns
        self.subplots = []
        self.plot_grid = None
        self.x_title = x_title
        if self.x_title is None:
            if self.variable == 'psi':
                self.x_title = r'Estimated $p({0} = 1 \, | \, B_{{\epsilon}}(S*))$'.format(self.variable_symbol)
            elif self.variable == 'omega':
                self.x_title = r'Estimated $p({0} < 0.01 \, | \, B_{{\epsilon}}(S*))$'.format(self.variable_symbol)
            elif self.variable == 'tau_exclusion':
                matplotlib.rc('mathtext',**{'fontset': 'stix'})
                self.x_title = (r'Estimated $p(\mathbf{\tau} \, \notin \, '
                        r'M \, | \, B_{\epsilon}(S*))$')
            else:
                raise Exception('unexpected variable {0!r}'.format(
                        self.variable))
        self.y_title = y_title
        self.width = width
        self.height = height
        self.auto_height = auto_height
        self.auto_adjust_margins = auto_adjust_margins
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_bottom = margin_bottom
        self.margin_top = margin_top
        self.padding_between_horizontal = padding_between_horizontal
        self.padding_between_vertical = padding_between_vertical
        self.bins = list(frange(0, 1, 20, include_end_point = True))
        self.tab = tab
        self.text_size = text_size
        self.bf_line = VerticalLine(x = self.bayes_factor_prob,
                color = self.bayes_factor_line_color)
        self.populate_subplots()

    def populate_subplots(self):
        for cfg, estimates in self.config_estimates_tups:
            dist = r'$\tau \sim {0}$'.format(str(cfg.tau))
            p = None
            prob = None
            if self.bayes_factor_prob is not None:
                count = len([e for e in estimates if e > self.bayes_factor_prob])
                p = count / float(len(estimates))
            if p is not None:
                if self.variable == 'psi':
                    prob = (r'$p(BF_{{{0} = 1, {0} \neq 1}} > {1}) = '
                            '{2}$'.format(self.variable_symbol,
                                    int(self.bayes_factor),
                                    p))
                elif self.variable == 'omega':
                    prob = (r'$p(BF_{{{0} < 0.01, {0} \geq 0.01}} > {1}) '
                            '= {2}$'.format(self.variable_symbol,
                                    int(self.bayes_factor),
                                    p))
                elif self.variable == 'tau_exclusion':
                    if self.cfg_to_prob_of_bf_exclusion:
                        p = self.cfg_to_prob_of_bf_exclusion[cfg]
                    prob = (r'$p(BF_{{\mathbf{{\tau}} \, \notin \, '
                            r'M, \mathbf{{\tau}} \, \notin \, '
                            r'M}} > {0}) = {1}$'.format(
                                    int(self.bayes_factor), p))
            hd = HistData(x = estimates,
                    normed = True,
                    bins = self.bins,
                    histtype = 'bar',
                    align = 'mid',
                    orientation = 'vertical',
                    zorder = 0)
            # b = list(frange(0, 1, 20, include_end_point = True))
            # counts = [0] * 20
            # for i in range(1, len(b)):
            #     for e in estimates:
            #         if (e >= b[i - 1]) and (e < b[i]):
            #             counts[i-1] += 1
            # densities = [((c * 20) / float(len(estimates))) for c in counts]
            # s = ScatterPlot()
            # d, bins, patches = hd.plot(s.ax)
            # for i, v in enumerate(d):
            #     assert almost_equal(v, densities[i])
            tick_labels = []
            for i, x in enumerate(self.bins[0:-1]):
                if not (i - 1) % 4:
                    tick_labels.append(str(x))
                else:
                    tick_labels.append('')
            xticks_obj = Ticks(ticks = self.bins,
                    labels = tick_labels,
                    horizontalalignment = 'center',
                    size = 10.0)
            hist = ScatterPlot(hist_data_list = [hd],
                    vertical_lines = [self.bf_line],
                    left_text = dist,
                    right_text = prob,
                    xticks_obj = xticks_obj,
                    tab = self.tab)
            hist.left_text_size = self.text_size
            hist.right_text_size = self.text_size
            hist.set_xlim(left = (self.bins[0]), right = (self.bins[-1]))
            self.subplots.append(hist)

    def create_grid(self):
        if len(self.subplots) < 2:
            self.num_columns = 1
        self.plot_grid = PlotGrid(subplots = self.subplots,
                num_columns = self.num_columns,
                share_x = True,
                share_y = False,
                label_schema = 'uppercase',
                title = self.x_title,
                title_top = False,
                y_title = self.y_title,
                width = self.width,
                height = self.height,
                auto_height = self.auto_height)
        self.plot_grid.auto_adjust_margins = self.auto_adjust_margins
        self.plot_grid.margin_left = self.margin_left
        self.plot_grid.margin_bottom = self.margin_bottom
        self.plot_grid.margin_right = self.margin_right
        self.plot_grid.margin_top = self.margin_top
        self.plot_grid.padding_between_horizontal = \
                self.padding_between_horizontal
        self.plot_grid.padding_between_vertical = self.padding_between_vertical
        self.plot_grid.reset_figure()
        return self.plot_grid

class AccuracyPowerPlotGrid(object):
    def __init__(self,
            observed_config_to_estimates,
            num_columns = 2,
            variable_symbol = r'\Omega',
            y_title_position = 0.006,
            width = 8,
            height = 9,
            auto_height = False,
            auto_adjust_margins = False,
            margin_left = 0.025,
            margin_bottom = 0.02,
            margin_right = 1,
            margin_top = 0.975,
            padding_between_horizontal = 0.5,
            padding_between_vertical = 1.4,
            tab = 0.08):
        self.config_estimates_tups = sorted(
            [(c, (list(e['x']), list(e['y']))) for c, e in observed_config_to_estimates.iteritems()],
            key = lambda x : x[0].tau.mean)
        self.num_columns = num_columns
        self.subplots = []
        self.plot_grid = None
        self.variable_symbol = variable_symbol
        self.x_title = r'True ${0}$'.format(self.variable_symbol)
        self.y_title = r'$\hat{{{0}}}$'.format(self.variable_symbol)
        self.y_title_position = y_title_position
        self.width = width
        self.height = height
        self.auto_height = auto_height
        self.auto_adjust_margins = auto_adjust_margins
        self.margin_left = margin_left
        self.margin_right = margin_right
        self.margin_bottom = margin_bottom
        self.margin_top = margin_top
        self.padding_between_horizontal = padding_between_horizontal
        self.padding_between_vertical = padding_between_vertical
        self.tab = tab
        self.populate_subplots()

    def populate_subplots(self):
        for cfg, estimates in self.config_estimates_tups:
            dist = r'$\tau \sim {0}$'.format(str(cfg.tau))
            assert len(estimates) == 2
            x, y = estimates[0], estimates[1]
            rmse = root_mean_square_error(x, y)
            assert len(x) == len(y)
            c = len([1 for i in range(len(x)) if y[i] < x[i]])
            p = c / float(len(x))
            prob = r'$p(\hat{{{0}}} < {0}) = {1}$'.format(
                    self.variable_symbol,
                    p)
            rmse_str = r'$RMSE = {0:.2f}$'.format(rmse)
            mx = max(x + y)
            mn = min(x + y)
            buff = (mx - mn) * 0.04
            xlim = (mn - buff, mx + buff)
            ylim = xlim
            sd = ScatterData(x = x, y = y,
                    marker = 'o',
                    markerfacecolor = 'none',
                    markeredgecolor = '0.35',
                    markeredgewidth = 0.7,
                    linestyle = '',
                    zorder = 100)
            sp = ScatterPlot(scatter_data_list = [sd],
                    left_text = dist,
                    right_text = prob,
                    xlim = xlim,
                    ylim = ylim,
                    identity_line = True,
                    tab = self.tab)
            xticks = [i for i in sp.ax.get_xticks()]
            yticks = [i for i in sp.ax.get_yticks()]
            xtick_labels = [i for i in xticks]
            if len(xtick_labels) >= 10:
                for i in range(1, len(xtick_labels), 2):
                    xtick_labels[i] = ''
            xticks_obj = Ticks(ticks = xticks,
                    labels = xtick_labels,
                    horizontalalignment = 'center',
                    size = 10.0)
            yticks_obj = Ticks(ticks = yticks,
                    labels = yticks,
                    size = 10.0)
            sp.xticks_obj = xticks_obj
            sp.yticks_obj = yticks_obj
            sp.left_text_size = 12.0
            sp.right_text_size = 12.0
            self.subplots.append(sp)

    def create_grid(self):
        if len(self.subplots) < 2:
            self.num_columns = 1
        self.plot_grid = PlotGrid(subplots = self.subplots,
                num_columns = self.num_columns,
                share_x = False,
                share_y = False,
                label_schema = 'uppercase',
                title = self.x_title,
                title_top = False,
                y_title = self.y_title,
                y_title_position = self.y_title_position,
                width = self.width,
                height = self.height,
                auto_height = self.auto_height)
        self.plot_grid.auto_adjust_margins = self.auto_adjust_margins
        self.plot_grid.margin_left = self.margin_left
        self.plot_grid.margin_bottom = self.margin_bottom
        self.plot_grid.margin_right = self.margin_right
        self.plot_grid.margin_top = self.margin_top
        self.plot_grid.padding_between_horizontal = \
                self.padding_between_horizontal
        self.plot_grid.padding_between_vertical = self.padding_between_vertical
        self.plot_grid.reset_figure()
        return self.plot_grid

class SaturationPlotGrid(object):
    def __init__(self, d,
            x_key = 'PRI.t',
            y_keys = None,
            y_labels = {},
            num_columns = 2,
            vertical_line_positions = [],
            label_schema = 'uppercase'):
        if not d.has_key(x_key):
            raise ValueError('`x_key` {0!r} is not in data dict'.format(x_key))
        self.d = d
        self.x_key = x_key
        self.y_keys = y_keys
        if not self.y_keys:
            self.y_keys = d.keys()
            self.y_keys.pop(self.x_key)
        for y in self.y_keys:
            if not self.d.has_key(y):
                raise ValueError('y key {0!r} is not in data dict'.format(y))
        self.y_labels = y_labels
        if not self.y_labels:
            self.y_labels = {'pi': r'$\pi$',
                             'pi.net': r'$\pi_{net}$',
                             'wattTheta': r'$\theta_W$',
                             'tajD.denom': r'$SD(\pi - \theta_W)$'}
        self.num_columns = num_columns
        self.v_lines = []
        for x in vertical_line_positions:
            self.v_lines.append(VerticalLine(x))
        self.label_schema = label_schema
        self.subplots = []
        self.plot_grid = None
        self.x_title = r'Divergence time $\tau$ in $4N_C$ generations'
        self.populate_subplots()

    def populate_subplots(self):
        for i, y in enumerate(self.y_keys):
            y_lab = self.y_labels.get(y, y)
            scatter_data = ScatterData(x = self.d[self.x_key], y = self.d[y])
            sp = ScatterPlot(scatter_data_list = [scatter_data],
                    vertical_lines = self.v_lines,
                    y_label = y_lab)
            self.subplots.append(sp)

    def create_grid(self):
        if len(self.subplots) < 2:
            self.num_columns = 1
        self.plot_grid = PlotGrid(subplots = self.subplots,
                num_columns = self.num_columns,
                share_x = True,
                share_y = False,
                label_schema = self.label_schema,
                title = self.x_title,
                title_top = False)
        return self.plot_grid

class SimResult(object):
    def __init__(self):
        self.true = []
        self.mode = []
        self.mode_glm = []
        self.median = []
        self.prob = []
        self.prob_glm = []
        self.validation_probs = None
        self.validation_probs_glm = None
        self.true_jitter = None
        self.mode_jitter = None
        self.mode_glm_jitter = None

    def jitter_true_and_modes(self, sigma = 0.005, rng = None):
        if not rng:
            rng = GLOBAL_RNG
        self.true_jitter = []
        self.mode_jitter = []
        self.mode_glm_jitter = []
        for i in range(len(self.true)):
            self.true_jitter.append(self.true[i] + rng.normalvariate(0.0,
                    sigma))
            self.mode_jitter.append(self.mode[i] + rng.normalvariate(0.0,
                    sigma))
            self.mode_glm_jitter.append(self.mode_glm[i] + rng.normalvariate(
                    0.0, sigma))

    def get_random_subsample(self, n, rng = None):
        l = len(self.true)
        if l <= n:
            raise ValueError('subsample is as big or bigger than population')
        if not rng:
            rng = GLOBAL_RNG
        indices = rng.sample(range(l), n)
        sr = SimResult()
        if self.true:
            sr.true = [self.true[i] for i in indices]
        if self.mode:
            sr.mode = [self.mode[i] for i in indices]
        if self.mode_glm:
            sr.mode_glm = [self.mode_glm[i] for i in indices]
        if self.median:
            sr.median = [self.median[i] for i in indices]
        if self.prob:
            sr.prob = [self.prob[i] for i in indices]
        if self.prob_glm:
            sr.prob_glm = [self.prob_glm[i] for i in indices]
        sr.validation_probs = self.validation_probs
        sr.validation_probs_glm = self.validation_probs_glm
        return sr
        
class ValidationResult(object):
    def __init__(self, result_paths = [],
            psi_of_interest = 1,
            omega_threshold = 0.01,
            omega_symbol = r'\Omega',
            psi_symbol = r'\Psi',
            mean_time_symbol = r'E(\tau)',
            math_font = None):
        self.psi_of_interest = psi_of_interest
        self.omega_threshold = omega_threshold
        self.omega_symbol = omega_symbol
        self.psi_symbol = psi_symbol
        self.mean_time_symbol = mean_time_symbol
        self.math_font = math_font
        self.psi = SimResult()
        self.omega = SimResult()
        self.tau = SimResult()
        self.prob_plot = None
        self.accuracy_plot = None
        if result_paths:
            self._parse_result_summary(result_paths)
    
    def get_random_subsample(self, n, rng = None):
        if not rng:
            rng = GLOBAL_RNG
        vr = ValidationResult(psi_of_interest = self.psi_of_interest,
                omega_threshold = self.omega_threshold)
        vr.psi = self.psi.get_random_subsample(n, rng)
        vr.omega = self.omega.get_random_subsample(n, rng)
        vr.tau = self.tau.get_random_subsample(n, rng)
        return vr

    def _parse_result_summary(self, result_paths):
        psi_prob_str = 'psi_{0}_prob'.format(self.psi_of_interest)
        psi_prob_glm_str = psi_prob_str + '_glm'
        psi_true_prob_glm_triples = []
        omega_true_prob_glm_triples = []
        for d in spreadsheet_iter(result_paths):
            psi_true = int(d['psi_true'])
            omega_true = float(d['omega_true'])
            self.psi.true.append(psi_true)
            self.omega.true.append(omega_true)
            self.psi.mode.append(int(d['psi_mode']))
            self.psi.mode_glm.append(float(d['psi_mode_glm']))
            self.omega.mode.append(float(d['omega_mode']))
            self.omega.median.append(float(d['omega_median']))
            self.omega.mode_glm.append(float(d['omega_mode_glm']))
            self.tau.true.append(float(d['mean_tau_true']))
            self.tau.mode.append(float(d['mean_tau_mode']))
            self.tau.median.append(float(d['mean_tau_median']))
            self.tau.mode_glm.append(float(d['mean_tau_mode_glm']))
            psi_true_prob_glm_triples.append((psi_true,
                    float(d[psi_prob_str]),
                    float(d[psi_prob_glm_str])))
            omega_true_prob_glm_triples.append((omega_true,
                    float(d['omega_prob_less']),
                    float(d['omega_prob_less_glm'])))
        self.psi.validation_probs = ValidationProbabilities(
                ((t, p) for (t, p, g) in psi_true_prob_glm_triples),
                true_value_of_interest = self.psi_of_interest,
                value_to_true_comparison = '==')
        self.psi.validation_probs_glm = ValidationProbabilities(
                ((t, g) for (t, p, g) in psi_true_prob_glm_triples),
                true_value_of_interest = self.psi_of_interest,
                value_to_true_comparison = '==')
        self.omega.validation_probs = ValidationProbabilities(
                ((t, p) for (t, p, g) in omega_true_prob_glm_triples),
                true_value_of_interest = self.omega_threshold,
                value_to_true_comparison = '<')
        self.omega.validation_probs_glm = ValidationProbabilities(
                ((t, g) for (t, p, g) in omega_true_prob_glm_triples),
                true_value_of_interest = self.omega_threshold,
                value_to_true_comparison = '<')
        self.prob_plot = ProbabilityValidationPlotGrid(
                psi_validation_probs = self.psi.validation_probs,
                psi_validation_probs_glm = self.psi.validation_probs_glm,
                omega_validation_probs = self.omega.validation_probs,
                omega_validation_probs_glm = self.omega.validation_probs_glm,
                omega_symbol = self.omega_symbol,
                psi_symbol = self.psi_symbol,
                math_font = self.math_font)
        self.accuracy_plot = AccuracyValidationPlotGrid(self,
                omega_symbol = self.omega_symbol,
                psi_symbol = self.psi_symbol,
                mean_time_symbol = self.mean_time_symbol,
                math_font = self.math_font)
    
    def save_prob_plot(self, path):
        fig = self.prob_plot.create_grid()
        fig.savefig(path)

    def save_accuracy_plot(self,path):
        fig = self.accuracy_plot.create_grid()
        fig.savefig(path)

def plot_validation_results(info_path, observed_indices = None,
        prior_indices = None,
        omega_symbol = r'\Omega',
        psi_symbol = r'\Psi',
        mean_time_symbol = r'E(\tau)',
        math_font = None):
    results = DMCSimulationResults(info_path)
    result_dir = os.path.dirname(info_path)
    plot_dir = os.path.join(result_dir, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    for obs_idx, obs_cfg in results.observed_index_to_config.iteritems():
        if observed_indices:
            if not obs_idx in observed_indices:
                continue
        for p_idx, p_cfg in results.prior_index_to_config.iteritems():
            if prior_indices:
                if not p_idx in prior_indices:
                    continue
            obs_name = os.path.splitext(os.path.basename(obs_cfg))[0]
            p_name = os.path.splitext(os.path.basename(p_cfg))[0]
            prob_plot_name = '_'.join([obs_name, p_name, 'mc_behavior']) + '.pdf'
            acc_plot_name = '_'.join([obs_name, p_name, 'accuracy']) + '.pdf'
            prob_plot_path = os.path.join(plot_dir, prob_plot_name)
            acc_plot_path = os.path.join(plot_dir, acc_plot_name)
            result_path = results.get_result_summary_path(obs_idx, p_idx)
            vr = ValidationResult([result_path],
                    omega_symbol = omega_symbol,
                    psi_symbol = psi_symbol,
                    mean_time_symbol = mean_time_symbol,
                    math_font = math_font)
            vr.save_prob_plot(prob_plot_path)
            vr.save_accuracy_plot(acc_plot_path)

def plot_model_choice_validation_results(info_path):
    results = DMCSimulationResults(info_path)
    result_dir = os.path.dirname(info_path)
    plot_dir = os.path.join(result_dir, 'plots')
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    result_paths = []
    for obs_idx in results.observed_index_to_config.iterkeys():
        result_paths.append(results.get_result_summary_path(obs_idx,
                results.combined_prior_index))
    vr = ValidationResult(result_paths)
    vr.save_prob_plot(os.path.join(plot_dir, 'mc_behavior.pdf'))
    vr.save_accuracy_plot(os.path.join(plot_dir, 'accuracy.pdf'))

