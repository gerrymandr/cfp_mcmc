import os
import sys
import matplotlib.pyplot as plt
import numpy as np


def get_x_bounds(xs, nums):
    nonzero_indices = np.argwhere(nums)
    return xs[nonzero_indices[0, 0]], xs[nonzero_indices[-2, 0]]


def plot(xs, nums, x_bounds, width, vline, line_color='red', fig=None):
    if not fig:
        fig = plt.figure()
    plt.bar(xs, nums, width=width, color='#aaaaaa')
    plt.plot(xs, nums, color='black')
    plt.xlim(*x_bounds)
    plt.ylim(0, 1.1 * nums[5:-5].max())
    plt.axvline(vline, color=line_color)
    ax = plt.gca()
    from matplotlib import ticker

    def my_formatter_fun(x, *args):
        return "%.0E" % x if x > 100 else x
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(my_formatter_fun))
    return fig


def load_hist_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        nbins, start_val, bin_width = lines[1].split(', ')
        nums = str.split(lines[3], ', ')
        nums = [int(x) for x in nums[:-1]]
    return int(nbins), float(start_val), float(bin_width), nums


def plot_from_file(filename, vline_location):
    nbins, start_val, bin_width, nums = load_hist_data(filename)
    xs = np.arange(nbins+1) * bin_width + start_val
    fig = plot(xs, np.array(nums), get_x_bounds(xs, nums), bin_width, float(vline_location))
    return fig


if __name__ == '__main__':
    # args:
    # 1: file
    # 2: red_line_loc
    # 3: output filename
    fig = plot_from_file(sys.argv[1], sys.argv[2])
    plt.savefig(os.path.expanduser(sys.argv[3]))

    exit(0)
