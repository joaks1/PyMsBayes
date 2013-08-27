#! /usr/bin/env python

import os
import sys
import math
import string

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
            zorder = 1,
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

class ScatterPlot(object):
    def __init__(self, scatter_data_list = [],
            vertical_lines = [],
            horizontal_lines = [],
            plot_label = None,
            x_label = None,
            y_label = None,
            left_text = None,
            center_text = None,
            right_text = None,
            marker = 'o',
            marker_face_color = 'none',
            marker_edge_color = '0.35',
            marker_edge_width = 0.7,
            position = (1,1,1)):
        self.scatter_data_list = list(scatter_data_list)
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
        self.marker = marker
        self.marker_face_color = marker_face_color
        self.marker_edge_color = marker_edge_color
        self.marker_edge_width = marker_edge_width
        self.shared_x_ax = None
        self.shared_y_ax = None
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
        for v in self.vertical_lines:
            self._plot_v_line(v)
        for h in self.horizontal_lines:
            self._plot_h_line(h)

    def _plot_scatter_data(self, d):
        l = self.ax.plot(d.x, d.y)
        plt.setp(l,
                marker = d.marker,
                linestyle = d.linestyle,
                markerfacecolor = d.markerfacecolor,
                markeredgecolor = d.markeredgecolor,
                markeredgewidth = d.markeredgewidth,
                zorder = d.zorder,
                **d.kwargs)

    def append_plot(self):
        self._plot()
        self.adjust_text_objects()

    def get_top_text_baseline(self):
        ymin, ymax = self.ax.get_ylim()
        return ymax + (math.fabs(ymax - ymin) * 0.01)

    def get_tab_indent(self):
        xmin, xmax = self.ax.get_xlim()
        return xmin + (math.fabs(xmax - xmin) * 0.06)

    def get_x_center(self):
        xmin, xmax = self.ax.get_xlim()
        return xmin + (math.fabs(xmax - xmin) * 0.5)

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
                **kwargs)

    def set_xlim(self, left = None, right = None, emit = True, auto = False,
            **kwargs):
        self.ax.set_xlim(left = left, right = right, emit = emit, auto = auto,
                **kwargs)

    def set_ylim(self, bottom = None, top = None, emit = True, auto = False,
            **kwargs):
        self.ax.set_xlim(bottom = bottom, top = top, emit = emit, auto = auto,
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
        self.ax.axvline(x = v.x, 
                ymin = v.ymin,
                ymax = v.ymax, 
                color = v.color,
                label = v.label,
                linestyle = v.linestyle,
                linewidth = v.linewidth,
                **v.kwargs)

    def add_v_line(self, vertical_line_object):
        self.vertical_lines.append(vertical_line_object)
        self._plot_v_line(vertical_line_object)

    def _plot_h_line(self, h):
        self.ax.axhline(y = h.y, 
                xmin = h.xmin,
                xmax = h.xmax, 
                color = h.color,
                label = h.label,
                linestyle = h.linestyle,
                linewidth = h.linewidth,
                **h.kwargs)

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
            title_top = True):
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
        self.fig = plt.figure()
        self.subplots = subplots
        for sp in self.subplots:
            sp.set_figure(self.fig)
        self.reset_figure()

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
        self.fig = plt.figure()
        plot_labels = self.get_plot_labels()
        nrows = self.get_num_rows()
        ncols = self.num_columns
        for i, subplot in enumerate(self.subplots):
            if subplot.get_figure() != self.fig:
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
        rect = (0, 0, 1, 0.98)
        if self.title:
            if self.title_top:
                self.fig.suptitle(self.title,
                        verticalalignment = 'top',
                        horizontalalignment = 'center',
                        y = 0.999)
                rect = (0, 0, 1, 0.94)
            else:
                self.fig.suptitle(self.title,
                        verticalalignment = 'bottom',
                        horizontalalignment = 'center',
                        y = 0.001)
                rect = (0, 0.04, 1, 0.98)
        self.fig.tight_layout(pad = 0.25, # outside margin
                h_pad = 0.8, # vertical padding between subplots
                w_pad = None, # horizontal padding between subplots
                rect = rect) # available space on figure

    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)


def saturation_plot(d, x_key = 'PRI.t', y_keys = None, y_labels = {},
        num_columns = 2, vertical_line_positions = []):
    if not d.has_key(x_key):
        raise ValueError('`x_key` {0!r} is not in data dict'.format(x_key))
    if not y_keys:
        y_keys = d.keys()
        y_keys.pop(x_key)
    for y in y_keys:
        if not d.has_key(y):
            raise ValueError('y key {0!r} is not in data dict'.format(y_key))
    if not y_labels:
        y_labels = {'pi': r'$\pi$',
                       'pi.net': r'$\pi_{net}$',
                       'wattTheta': r'$\theta_W$',
                       'tajD.denom': r'$SD(\pi - \theta_W)$'}
    subplots = []
    v_lines = []
    for x in vertical_line_positions:
        v_lines.append(VerticalLine(x))
    for i, y in enumerate(y_keys):
        y_lab = y_labels.get(y, y)
        scatter_data = ScatterData(x = d[x_key], y = d[y])
        sp = ScatterPlot(scatter_data_list = [scatter_data],
                vertical_lines = v_lines,
                y_label = y_lab)
        subplots.append(sp)

    return PlotGrid(subplots = subplots,
            num_columns = num_columns,
            share_x = True,
            share_y = False,
            label_schema = 'uppercase',
            title = r'Divergence time $\tau$ in $4N_C$ generations',
            title_top = False)
    
