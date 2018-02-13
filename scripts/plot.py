import sys
import matplotlib.pyplot as plt
import numpy as np


def plot(xs, nums, width, vline):
    plt.bar(xs, nums, width=width)
    plt.axvline(vline, color='r')
    plt.show()


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
    plot(xs, nums, bin_width, float(vline_location))


if __name__ == '__main__':
    # args:
    # 1: file
    # 2: red_line_loc
    plot_from_file(sys.argv[1], sys.argv[2])

    exit(0)
