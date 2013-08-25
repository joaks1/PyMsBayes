#! /usr/bin/env python

import os
import sys

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


class ScatterPlot(object):
    def __init__(self, x = None, y = None,
            plot_label = None,
            x_label = None,
            y_label = None,
            left_text = None,
            center_text = None,
            right_text = None):
        self.x = x
        self.y = y
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1,1,1)
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
        self.reset_plot()

    def clear(self):
        self.ax.clear()

    def reset_plot(self):
        self.clear()
        self._plot()
        self._reset_text_objects()

    def _plot(self):
        self.ax.scatter(self.x, self.y)
    
    def append_plot(self):
        self._plot()
        self.adjust_text_objects()

    def get_top_text_baseline(self):
        ymin, ymax = self.ax.get_ylim()
        return ymax + (math.fabs(ymax - ymin) * 0.01)

    def get_tab_indent(self):
        xmin, xmax = self.ax.get_xlim()
        return xmin + (math.fabs(xmax - xmin) * 0.02)

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

    def set_plot_label(label = None, fontdict = None, withdash = False,
            **kwargs):
        if label:
            self._plot_label = label
        if self._plot_label is None:
            return
        self.remove_text_object(self.text_objects['plot_label'])
        verticalalignment = 'bottom'
        horizontalalignment = 'left'
        x = self.get_xmin()
        y = self.get_top_text_baseline()
        self.text_objects['plot_label'] = self.ax.text(x = x, y = y,
                s = self._plot_label,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)

    def set_left_text(left_text = None, fontdict = None, withdash = False,
            **kwargs):
        if left_text:
            self._left_text = left_text
        if self._left_text is None:
            return
        self.remove_text_object(self.text_objects['left_text'])
        verticalalignment = 'bottom'
        horizontalalignment = 'left'
        x = self.get_tab_indent()
        y = self.get_top_text_baseline()
        self.text_objects['left_text'] = self.ax.text(x = x, y = y,
                s = self._left_text,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)

    def set_center_text(center_text = None, fontdict = None, withdash = False,
            **kwargs):
        if center_text:
            self._center_text = center_text
        if self._center_text is None:
            return
        self.remove_text_object(self.text_objects['center_text'])
        verticalalignment = 'bottom'
        horizontalalignment = 'center'
        x = self.get_x_center()
        y = self.get_top_text_baseline()
        self.text_objects['center_text'] = self.ax.text(x = x, y = y,
                s = self._center_text,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)

    def set_right_text(right_text = None, fontdict = None, withdash = False,
            **kwargs):
        if right_text:
            self._right_text = right_text
        if self._right_text is None:
            return
        self.remove_text_object(self.text_objects['right_text'])
        verticalalignment = 'bottom'
        horizontalalignment = 'right'
        x = self.get_xmax()
        y = self.get_top_text_baseline()
        self.text_objects['right_text'] = self.ax.text(x = x, y = y,
                s = self._right_text,
                fontdict = fontdict,
                withdash = withdash, 
                verticalalignment = verticalalignment,
                horizontalalignment = horizontalalignment,
                **kwargs)

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

    def _adjust_left_text(self):
        if not self.text_objects.get('left_text', None):
            return
        x = self.get_tab_indent()
        y = self.get_top_text_baseline()
        self.text_objects['left_text'].set_position((x, y))

    def _adjust_center_text(self):
        if not self.text_objects.get('center_text', None):
            return
        x = self.get_x_center()
        y = self.get_top_text_baseline()
        self.text_objects['center_text'].set_position((x, y))

    def _adjust_right_text(self):
        if not self.text_objects.get('right_text', None):
            return
        x = self.get_xmax()
        y = self.get_top_text_baseline()
        self.text_objects['right_text'].set_position((x, y))

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

    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)

# class ScatterPlotGrid(object):
#     def __init__(self,
#             xy_pairs,
#             shared_x = False,
#             marker_face_color = 'none',
#             marker_edge_color = '0.4',
#             marker_edge_width = 0.7):
#         pass
#         self.xy_pairs

    # axis_plot_labels = {'pi': r'$\pi$',
    #                'pi.net': r'$\pi_{net}$',
    #                'wattTheta': r'$\theta_W$',
    #                'tajD.denom': r'$SD(\pi - \theta_W)$'}
    # fig = plt.figure()
    # ncols = 2
    # nrows = get_num_rows(len(stat_keys))
    # for i, k in enumerate(stat_keys):
    #     if fig.axes:
    #         ax = fig.add_subplot(nrows, ncols, i + 1, sharex = fig.axes[0])
    #     else:
    #         ax = fig.add_subplot(nrows, ncols, i + 1)
    #     ax.plot(stats_by_time['PRI.t'], stats_by_time[k])
    #     ax.set_ylabel(axis_plot_labels[k])
    # plt.setp([a.lines for a in fig.axes],
    #         marker = 'o',
    #         linestyle='',
    #         markerfacecolor = 'none',
    #         markeredgecolor = '0.4',
    #         markeredgewidth = 0.7)
    # fig.suptitle(r'Divergence time $\tau$ in $4N_C$ generations',
    #         verticalalignment = 'bottom',
    #         y = 0.001)
    # fig.tight_layout(pad = 0.25, # out side margin
    #                  h_pad = None, # height padding between subplots
    #                  w_pad = None, # width padding between subplots
    #                  rect = (0,0.05,1,1))
    # fig.savefig('plot.pdf')
