#! /usr/bin/env python

import os
import sys
import random
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import numpy
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.stats import gamma
import numpy as np

def get_rejection_plot_3d(maximum = 1.0,
        observed = (0.5, 0.5, 0.5),
        nsamples = 0,
        show_rejected = True,
        rng = None):
    if not rng:
        rng = random.Random()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot([observed[0]], [observed[1]], [observed[2]],
            marker = "o",
            zorder=100
            )
    ax.set_xlim((0.0, maximum))
    ax.set_ylim((0.0, maximum))
    ax.set_zlim((0.0, maximum))
    rect = [-0.12, 0, 1, 1.07]
    fig.tight_layout(pad = 0.25, rect = rect)
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    
    x = (0.1 * np.outer(np.cos(u), np.sin(v))) + observed[0]
    y = (0.1 * np.outer(np.sin(u), np.sin(v))) + observed[1]
    z = (0.1 * np.outer(np.ones(np.size(u)), np.cos(v))) + observed[2]
    if show_rejected:
        ax.plot_surface(x, y, z,  rstride=10, cstride=15, linewidth = 0.1,
                color='none',
                edgecolor='0.5',
                zorder=1000
                )
    if nsamples > 0:
        draws = [[], [], []]
        for i in range(nsamples):
            draws[0].append(rng.random())
            draws[1].append(rng.random())
            draws[2].append(rng.random())
        if not show_rejected:
            s = [[], [], []]
            for i in range(nsamples):
                e = math.sqrt(math.pow(draws[0][i] - observed[0], 2) +
                        math.pow(draws[1][i] - observed[1], 2) +
                        math.pow(draws[2][i] - observed[2], 2))
                if e < 0.1:
                    s[0].append(draws[0][i])
                    s[1].append(draws[1][i])
                    s[2].append(draws[2][i])
            draws = s
        ax.plot(draws[0], draws[1], draws[2],
                marker = 'o',
                markerfacecolor = 'black',
                markeredgecolor = 'black',
                markeredgewidth = 0.7,
                markersize = 0.8,
                linestyle = '',
                zorder = 500,
                )
    
    return ax, fig

def main_cli():
    maximum = 1.0
    seed = 111
    rng = random.Random()
    rng.seed(seed)
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = 0,
            show_rejected = False,
            rng = rng)
    fig.savefig('rejection-sampling-observed.png')
    rng.seed(seed)
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = 0,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-tolerance.png')
    rng.seed(seed)
    nsamples = 10
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-{0}.png'.format(nsamples))
    rng.seed(seed)
    nsamples = 100
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-{0}.png'.format(nsamples))
    rng.seed(seed)
    nsamples = 200
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-{0}.png'.format(nsamples))
    rng.seed(seed)
    nsamples = 500
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-{0}.png'.format(nsamples))
    rng.seed(seed)
    nsamples = 1000
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-{0}.png'.format(nsamples))
    rng.seed(seed)
    nsamples = 10000
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-{0}.png'.format(nsamples))
    rng.seed(seed)
    nsamples = 20000
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = True,
            rng = rng)
    fig.savefig('rejection-sampling-{0}.png'.format(nsamples))
    rng.seed(seed)
    nsamples = 20000
    ax, fig = get_rejection_plot_3d(maximum = maximum,
            observed = (0.5, 0.5, 0.5),
            nsamples = nsamples,
            show_rejected = False,
            rng = rng)
    fig.savefig('rejection-sampling-{0}-post.png'.format(nsamples))


if __name__ ==  '__main__':
    main_cli()

