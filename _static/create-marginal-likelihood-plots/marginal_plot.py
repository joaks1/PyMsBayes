#! /usr/bin/env python

import os
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import numpy
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy.stats import gamma

def get_bivariate_normal_and_uniform_densities(maximum = 1.0,
        mean = (0.15, 0.247),
        variance = (0.039, 0.026),
        covariance = 0.0,
        npoints = 100):
    a = numpy.linspace(0, maximum, npoints)
    b = numpy.linspace(0, maximum, npoints)
    X, Y = numpy.meshgrid(a, b)
    Z1 = get_bivariate_normal_density(X, Y,
            mean = mean,
            variance = variance,
            covariance = covariance)
    Z2 = (Z1 * 0.0) + (1.0 / (maximum ** 2))
    return X, Y, Z1, Z2

def get_bivariate_normal_density(x, y,
        mean = (2.0, 3.0),
        variance = (0.2, 0.2),
        covariance = 0.0):
    return matplotlib.mlab.bivariate_normal(x, y,
            sigmax = variance[0],
            sigmay = variance[1],
            mux = mean[0],
            muy = mean[1],
            sigmaxy = covariance)

def get_marginal_likelihood(x, y, z):
    max_x, max_y = 0.0, 0.0
    for i in x:
        max_x = max([max_x] + [max(i)])
    for i in y:
        max_y = max([max_y] + [max(i)])
    prior = 1.0 / (max_x * max_y)
    l = 0.0
    w = 0.0
    for i in range(len(z)):
        for j in range(len(z[0])):
            l += (z[i][j] * prior)
            w += prior
    return l/w

def get_marginal_likelihood_constrained(x, y, z):
    assert len(x) == len(y)
    max_x, max_y = 0.0, 0.0
    for i, a in enumerate(x):
        assert len(x[i]) == len(y[i])
        max_x = max([max_x] + [max(a)])
    for a in y:
        max_y = max([max_y] + [max(a)])
    assert max_x == max_y
    prior = 1.0 / max_x
    l = 0.0
    w = 0.0
    for i in range(len(z)):
        l += (z[i][i] * prior)
        w += prior
    return l/w

def get_marginal_plot_2d(maximum = 1.0,
        likelihood_shape = 50.0,
        likelihood_scale = 0.002,
        prior_shape = 3.0,
        prior_scale = 0.06,
        npoints = 500,
        include_uniform_prior = True,
        include_gamma_prior = True,
        include_likelihood = True,
        include_function_labels = True,
        parameter_symbol = r'\theta',
        linewidth = 2.0,
        prior_label_x = 0.5):
    x = numpy.linspace(0.0000001, maximum, npoints)
    likelihood = gamma(likelihood_shape, scale = likelihood_scale)
    y = [likelihood.pdf(i) for i in x]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    if include_likelihood:
        likelihood_line = ax.plot(x, y)
        plt.setp(likelihood_line,
                color = '0.3',
                linestyle = '-',
                linewidth = linewidth,
                marker = '',
                zorder = 200)
        max_idx = y.index(max(y)) + int(round(0.01 * npoints))
        label_target = (x[max_idx], y[max_idx])
        label_position = (label_target[0] + (0.08 * maximum), label_target[1])
        if include_function_labels:
            plt.annotate(r'$p(X \mid \, {0})$'.format(parameter_symbol),
                    xy = label_target,
                    arrowprops = dict(arrowstyle = '->'),
                    xytext = label_position,
                    size = 18.0)
    prior_x = prior_label_x * maximum
    u_density = 1.0 / maximum
    ymax = max(ax.get_ylim())
    prior_label_position = (prior_x, u_density + (0.1 * ymax))
    if include_uniform_prior:
        u = [u_density for i in range(len(x))]
        u_line = ax.plot(x, u)
        plt.setp(u_line,
                color = 'r',
                linestyle = '-',
                linewidth = linewidth,
                marker = '',
                zorder = 0)
        u_label_target = (prior_x + (0.04 * maximum), u_density)
        if include_function_labels:
            plt.annotate(r'$p({0})$'.format(parameter_symbol),
                    xy = u_label_target,
                    arrowprops = dict(arrowstyle = '->'),
                    xytext = prior_label_position,
                    size = 18.0)
                    # verticalalignment = 'bottom',
                    # horizontalalignment = 'center')
    if include_gamma_prior:
        g_prior = gamma(prior_shape, scale = prior_scale)
        g = [g_prior.pdf(i) for i in x]
        g_line = ax.plot(x, g)
        plt.setp(g_line,
                color = 'b',
                linestyle = '-',
                linewidth = linewidth,
                marker = '',
                zorder = 100)
        idx = g.index(max(g)) + int(round(0.1 * npoints))
        g_label_target = (x[idx], g[idx])
        if include_function_labels:
            plt.annotate('',
                    xy = g_label_target,
                    arrowprops = dict(arrowstyle = '->'),
                    xytext = prior_label_position,
                    size = 18.0)
                    # verticalalignment = 'center',
                    # horizontalalignment = 'center')
    ax.set_xlabel(r'${0}$'.format(parameter_symbol), size=18.0)
    ax.set_ylabel(r'Density', size=18.0)
    rect = [0, 0, 1, 1]
    fig.tight_layout(pad = 0.25, rect = rect)
    return ax, fig

def get_marginal_plot_3d(maximum = 1.0,
        mean = (0.15, 0.247),
        variance = (0.039, 0.026),
        covariance = 0.0,
        npoints = 100,
        include_prior = True,
        include_constrained_density = True,
        linewidth = 0.1):
    X, Y, Z1, Z2 = get_bivariate_normal_and_uniform_densities(maximum = maximum,
            mean = mean,
            variance = variance,
            covariance = covariance,
            npoints = npoints)
    ml_2p = get_marginal_likelihood(X, Y, Z1)
    ml_1p = get_marginal_likelihood_constrained(X, Y, Z1)
    sys.stdout.write('marginal likelihood of 2-parameter model: {0}\n'.format(ml_2p))
    sys.stdout.write('marginal likelihood of 1-parameter model: {0}\n'.format(ml_1p))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, linewidth=linewidth, antialiased=False, shade=True, cmap=cm.coolwarm, zorder=200)
    if include_prior:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        prior_d = 1.0 / (xmax * ymax)
        prior_d *= 2.0
        x_back_line = ax.plot([xmin, xmax], [ymax, ymax], [prior_d, prior_d])
        x_front_line = ax.plot([xmin, xmax], [ymin, ymin], [prior_d, prior_d], zorder=200)
        y_back_line = ax.plot([xmin, xmin], [ymin, ymax], [prior_d, prior_d], zorder=-10)
        y_front_line = ax.plot([xmax, xmax], [ymin, ymax], [prior_d, prior_d], zorder=200)
        plt.setp([x_back_line, y_back_line, x_front_line, y_front_line],
                color = 'r',
                linestyle = '--',
                linewidth = 1.0,
                marker = '')
    if include_constrained_density:
        a, b, c = [], [], []
        for i in range(len(X)):
            a.append(X[i][i])
            b.append(Y[i][i])
            c.append(Z1[i][i])
        identity_line = ax.plot(a, b, c)
        plt.setp(identity_line,
                color = 'w',
                linestyle = '-',
                linewidth = 0.75,
                marker = '',
                zorder = 100)
    ax.set_xlabel(r'$T_1$', size=14.0)
    ax.set_ylabel(r'$T_2$', size=14.0)
    ax.set_zlabel('Density', size=14.0)
    rect = [-0.12, 0, 1, 1.07]
    fig.tight_layout(pad = 0.25, rect = rect)
    return ax, fig

def main_cli():
    # maximum = 1.0
    # ax, fig = get_marginal_plot_3d(maximum = maximum,
    #         mean = (0.15, 0.247),
    #         variance = (0.039, 0.026),
    #         covariance=0.0,
    #         npoints = 100,
    #         include_prior = True,
    #         include_constrained_density = True,
    #         linewidth=0.1)
    # fig.savefig('marginal-plot-3d.png', dpi=300)

    # ax, fig = get_marginal_plot_3d(maximum = maximum,
    #         mean = (0.15, 0.247),
    #         variance = (0.039, 0.026),
    #         covariance=0.0,
    #         npoints = 100,
    #         include_prior = False,
    #         include_constrained_density = False)
    # fig.savefig('marginal-plot-3d-bare.png', dpi=300)

    # ax, fig = get_marginal_plot_3d(maximum = maximum,
    #         mean = (0.15, 0.247),
    #         variance = (0.039, 0.026),
    #         covariance=0.0,
    #         npoints = 100,
    #         include_prior = True,
    #         include_constrained_density = False)
    # fig.savefig('marginal-plot-3d-prior.png', dpi=300)

    # ax, fig = get_marginal_plot_2d(maximum = maximum,
    #     likelihood_shape = 50.0,
    #     likelihood_scale = 0.002,
    #     prior_shape = 3.0,
    #     prior_scale = 0.06,
    #     npoints = 500,
    #     include_uniform_prior = True,
    #     include_gamma_prior = True,
    #     linewidth = 2.0)
    # fig.savefig('marginal-plot-2d.png')

    # ax, fig = get_marginal_plot_2d(maximum = maximum,
    #     likelihood_shape = 50.0,
    #     likelihood_scale = 0.002,
    #     prior_shape = 3.0,
    #     prior_scale = 0.06,
    #     npoints = 500,
    #     include_uniform_prior = False,
    #     include_gamma_prior = False,
    #     linewidth = 2.0)
    # fig.savefig('marginal-plot-2d-no-priors.png')

    # ax, fig = get_marginal_plot_2d(maximum = maximum,
    #     likelihood_shape = 50.0,
    #     likelihood_scale = 0.002,
    #     prior_shape = 3.0,
    #     prior_scale = 0.06,
    #     npoints = 500,
    #     include_uniform_prior = True,
    #     include_gamma_prior = False,
    #     linewidth = 2.0)
    # fig.savefig('marginal-plot-2d-uniform-prior.png')
    #
    maximum = 10.0
    ax, fig = get_marginal_plot_3d(maximum = maximum,
            # mean = (2.46, 3.8),
            # variance = (0.58, 0.4314),
            mean = (2.46, 3.584),
            variance = (0.58, 0.4314),
            covariance=0.1,
            # covariance=0.0,
            npoints = 100,
            include_prior = True,
            include_constrained_density = True,
            linewidth=0.1)
    fig.savefig('../marginal-plot-3d.png', dpi=200)

    ax, fig = get_marginal_plot_2d(maximum = maximum,
        likelihood_shape = 8.0,
        likelihood_scale = 0.2,
        prior_shape = 2.0,
        prior_scale = 1.0,
        npoints = 500,
        include_uniform_prior = True,
        include_gamma_prior = True,
        include_likelihood = True,
        include_function_labels = True,
        parameter_symbol = r'T',
        linewidth = 2.0)
    fig.savefig('../marginal-plot-2d.png')


if __name__ ==  '__main__':
    main_cli()

