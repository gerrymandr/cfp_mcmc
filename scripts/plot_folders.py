from collections import defaultdict
from glob import glob
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re


def get_x_bounds(xs, nums):
    nonzero_indices = np.argwhere(nums)
    return xs[nonzero_indices[0, 0]], xs[nonzero_indices[-2, 0]]


def plot(xs, nums, x_bounds, width, vline, line_color, fig=None):
    if not fig:
        fig = plt.figure()
    plt.bar(xs, nums, width=width, color='#aaaaaa')
    plt.plot(xs, nums, color='black')
    plt.xlim(*x_bounds)
    plt.ylim(0, 1.1 * nums[5:-5].max())
    plt.axvline(vline, color=line_color)
    return fig


def generate_plot(dfs, initial_values, metric, race, district, subtract_1E5, out_folder):
    fig = plt.figure()
    values = dfs[metric][race][district]
    if subtract_1E5:
        values = values - dfs[metric][race][district + '_1E5'].values
    initial_value = initial_values[metric][race][district]
    width = values.index[1] - values.index[0]

    # inefficient to recompute this within the loop, but whatever.
    x_bounds = []
    for r, districts_dicts in dfs[metric].items():
        # same race as our current one
        if r[:3] == race[:3]:
            for d, df in districts_dicts.items():
                x_bounds.append(get_x_bounds(df.index.values, df[d].values))
                # add in initial value as potential max bound (skip lower bound by making it absurd)
                if d in initial_values[metric][r]:
                    x_bounds.append((99999999, initial_values[metric][r][d]))
    mins, maxes = zip(*x_bounds)
    x_bounds = min(mins), max(maxes)

    line_color = 'red' if district.endswith('Rep') else 'purple'

    plot(values.index.values, values[district].values, x_bounds, width,
         initial_value, line_color, fig)
    fn = '_'.join((metric, race, district)) + '_sub1E5' if subtract_1E5 else '' + '.png'
    plt.title(''.join((race, district)))
    plt.savefig(os.path.join(out_folder, fn))
    plt.close()


def load_hist_data(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        nbins, start_val, bin_width = lines[1].split(', ')
        nums = str.split(lines[3], ', ')
        nums = [int(x) for x in nums[:-1]]
    return int(nbins), float(start_val), float(bin_width), nums


def get_attribute_files(dirname, races, district_ids):
    folders = []
    for race in races:
        folders.extend(glob(os.path.join(dirname, "{}*-*".format(race))))
    files = {'mm': defaultdict(dict),
             'eg': defaultdict(dict),
             'smoothed': defaultdict(dict)}
    initial_values = {'mm': defaultdict(dict),
                      'eg': defaultdict(dict),
                      'smoothed': defaultdict(dict)}
    for folder in folders:
        race, district_id = folder.split('-')
        race = os.path.basename(os.path.splitext(race)[0])
        district_id = os.path.splitext(district_id)[0]
        if district_id in district_ids:
            with open(os.path.join(folder, 'out_err.log')) as f:
                log_data = f.read()
            files['mm'][race][district_id] = os.path.join(folder,
                                                'median_mean_hist.txt')
            files['mm'][race][district_id + '_1E5'] = os.path.join(folder,
                                                'median_mean_hist_1E5.txt')
            initial_values['mm'][race][district_id] = float(
                re.search(r"initial.*median.*(\d\.\d+)", log_data).group(1))
            files['eg'][race][district_id] = os.path.join(folder,
                                                'eg_hist.txt')
            files['eg'][race][district_id + '_1E5'] = os.path.join(folder,
                                                'eg_hist_1E5.txt')
            initial_values['eg'][race][district_id] = float(
                re.search(r"initial.*gap.*(\d\.\d+)", log_data).group(1))
            files['smoothed'][race][district_id] = os.path.join(folder,
                                                'smoothed_seats_hist.txt')
            files['smoothed'][race][district_id + '_1E5'] = os.path.join(folder,
                                                'smoothed_seats_hist_1E5.txt')
            initial_values['smoothed'][race][district_id] = float(
                re.search(r"initial\ssmoothed.*(\d\.\d+)", log_data).group(1))
    return files, initial_values


def load_files_into_dfs(files_dict):
    dfs = defaultdict(dict)
    for metric, races_dicts in files_dict.items():
        races_dfs = dfs.get(metric, {})
        for race, district_dicts in races_dicts.items():
            districts_dfs = races_dfs.get(race, {})
            for district, filename in district_dicts.items():
                nbins, start_val, bin_width, bin_values = load_hist_data(filename)
                xs = np.arange(nbins+1) * bin_width + start_val
                df = pd.DataFrame({district: bin_values}, index=xs)
                districts_dfs[district] = df
            races_dfs[race] = districts_dfs
        dfs[metric] = races_dfs
    return dfs


def plot_from_file(filename, vline_location):
    nbins, start_val, bin_width, nums = load_hist_data(filename)
    xs = np.arange(nbins+1) * bin_width + start_val
    plot(xs, nums, bin_width, float(vline_location))


if __name__ == '__main__':
    # args:
    # 1: dirname
    # 2: output folder (defaults to . if not provided)
    races = 'cong', 'sen', 'pres'
    district_ids = 'CurrentRep', 'NewRep', 'Feb_11_Wes_units'
    out_folder = sys.argv[2] if len(sys.argv) > 2 else '.'
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)
    files, initial_values = get_attribute_files(sys.argv[1], races, district_ids)
    dfs = load_files_into_dfs(files)
    for metric, races_dicts in dfs.items():
        for race, dfs_dicts in races_dicts.items():
            for district in (k for k in dfs_dicts.keys() if not k.endswith('_1E5')):
                for subtract_1E5 in (True, False):
                    generate_plot(dfs, initial_values, metric, race, district, subtract_1E5, out_folder)
    exit(0)
