#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
@author: B. Font Garcia
@style: PEP8
@contact: b.fontgarcia@soton.ac.uk
"""

# Imports
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Change working directory to current script folder
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# File containig data
fname1 = '../fort.11'
fname2 = '../fort.12'
fname3 = '../fort.13'


def plot1D(data, file):
    """
    Generate a 1D plot using the matplotlib library given the arguments data = [x,y] 
    and the file name to save the plot in .pdf format (can be easly converted to .svg).
    """
    # Basic definitions
    plt.rcParams['text.usetex'] = True  # Set TeX interpreter
    ax = plt.gca()
    fig = plt.gcf()

    # Define data
    x = data[:, 0]
    y = data[:, 1]

    # Show lines
    plt.plot(x, y, color='blue', lw=1.5)
    # plt.show()

    # Set limits
    ax.set_xlim(min(x), max(x))
    ax.set_ylim(min(y), max(y)*1.2)

    # Make the plot square --------
    fwidth = fig.get_figwidth()
    fheight = fig.get_figheight()
    # get the axis size and position in relative coordinates
    # this gives a BBox object
    bb = ax.get_position()
    # calculate them into real world coordinates
    axwidth = fwidth*(bb.x1-bb.x0)
    axheight = fheight*(bb.y1-bb.y0)
    # if the axis is wider than tall, then it has to be narrowe
    if axwidth > axheight:
        # calculate the narrowing relative to the figure
        narrow_by = (axwidth-axheight)/fwidth
        # move bounding box edges inwards the same amount to give the correct width
        bb.x0 += narrow_by/2
        bb.x1 -= narrow_by/2
    # else if the axis is taller than wide, make it vertically smaller
    # works the same as above
    elif axheight > axwidth:
        shrink_by = (axheight-axwidth)/fheight
        bb.y0 += shrink_by/2
        bb.y1 -= shrink_by/2
    ax.set_position(bb)
    # --------

    # Edit frame, labels and legend
    ax.axhline(linewidth=1)
    ax.axvline(linewidth=1)
    # plt.xlabel(r'$x$')
    # plt.ylabel(r'$y$')
    # leg = plt.legend(loc='upper right')
    # leg.get_frame().set_edgecolor('black')

    # Show plot and save figure
    plt.savefig(file, transparent=True, bbox_inches='tight')
    plt.close()
    return


"""""""""""""""""""""""""""
Main function

"""""""""""""""""""""""""""


def main():

    global data, labels

#     h(x) = (f*g)(x)
#     f(x) : Function to convolve
#     g(x) : Kernel Gaussian function used for the convolution
#     h(x) : Function resulting from the convolution

    # Kernel function (g)
    data = np.genfromtxt(fname1)[1:, :]
    labels = np.genfromtxt(fname1, dtype=str)[0, :]
    x = data[:, 1]
    for i in range(2, np.size(data, 1)):
        plot1D(np.vstack((x, data[:, i])).T, labels[i]+'.pdf')

    # Convolution results (f, h)
    data = np.genfromtxt(fname2)[1:, :]
    labels = np.genfromtxt(fname2, dtype=str)[0, :]
    x = data[:, 1]
    for i in range(2, np.size(data, 1)):
        plot1D(np.vstack((x, data[:, i])).T, labels[i]+'.pdf')

    # Padded functions (f_pad, g_pad, h_pad)
    data = np.genfromtxt(fname3)[1:, :]
    labels = np.genfromtxt(fname3, dtype=str)[0, :]
    x = data[:, 1]
    for i in range(2, np.size(data, 1)):
        plot1D(np.vstack((x, data[:, i])).T, labels[i]+'.pdf')

if __name__ == '__main__':
    main()
