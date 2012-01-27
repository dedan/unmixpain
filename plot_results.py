#!/usr/bin/env python
# encoding: utf-8

from os import path
import pylab as plt
import numpy as np
import pickle
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import mpl
import matplotlib.gridspec as gridspec

data_folder = '/Users/dedan/projects/fu/data/unmixpain/'
plot_format = 'pdf'
out_folder = path.join(data_folder, 'out')
res = pickle.load(open(path.join(out_folder, 'results.pickle')))

subsamples = {'sub1':
              ['2010-03-19-u1', '2010-03-19-u2', '2010-03-19-u4'],
              'sub2':
              ['2010-12-02-u1', '2010-12-02-u2'],
              'sub3':
              ['2010-12-13-u3', '2010-12-13-u4', '2010-12-14-u4', '2010-12-14-u5']}

# first plot the results for each unit seperate
for key, result in res.items():

    fig_res = plt.figure()
    gs = gridspec.GridSpec(2, 2)
    gs.update(hspace=0.5)
    p_res_mech = fig_res.add_subplot(gs[0, 0])
    p_res_temp = fig_res.add_subplot(gs[0, 1])


    # plot the ISIs over time (x-axis normalized!)
    for i, isi in enumerate(result['isis']):
        if len(isi) > 3:
            bla = result['flevels'][i] / np.max(result['flevels'])
            c = cm.jet(bla, 1)
            p_res_mech.plot(np.array(range(len(isi))) / float(len(isi)-1), isi, '.-', color=c)
    p_res_mech.set_title('relative time vs. ISI')
    p_res_mech.set_xlabel('normalized time')
    p_res_mech.set_ylabel('ISI')

    # plot ISIs over temperature
    min_temp = result['temp_range'][0]
    max_temp = result['temp_range'][1]
    for i, temp_t in enumerate(result['temp_t']):
        c = cm.jet((temp_t - min_temp) / (max_temp - min_temp), 1)
        pp = p_res_temp.plot(i, result['temp_isis'][i], '.-', color=c)
    axins1 = inset_axes(p_res_temp, width="50%", height="5%", loc=1)
    norm = mpl.colors.Normalize(vmin=min_temp, vmax=max_temp)
    ticks = [round(min_temp, 2)+0.01, round(max_temp, 2)-0.01]
    mpl.colorbar.ColorbarBase(axins1, norm=norm, cmap=cm.jet,
                              orientation="horizontal",
                              ticks=ticks)
    p_res_temp.set_title('ISI over time \n (temperature coded by color)')
    p_res_temp.set_xlabel('time')
    p_res_temp.set_ylabel('ISI')
    fig_res.savefig(path.join(out_folder, 'fig_res_%s.%s') % (key, plot_format))

plt.show()


for subsample, units in subsamples.items():

    # plot CV and firing rate in relation to stimulation force level
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2)
    gs.update(hspace=0.5, wspace=0.5)
    p_cv = fig.add_subplot(gs[0, 0])
    p_cv.set_title('CV against force level')
    p_cv.set_xlabel('force')
    p_cv.set_ylabel('CV')
    p_rate = plt.subplot(gs[0, 1])
    p_rate.set_title('rate against force level')
    p_rate.set_xlabel('force')
    p_rate.set_ylabel('rate')
    p_temp = fig.add_subplot(gs[1, 0])
    p_temp.set_title('relative time vs. ISI')
    p_temp.set_xlabel('relative time')
    p_temp.set_ylabel('ISI')
    p_legend = plt.subplot(gs[1, 1])
    p_legend.get_xaxis().set_visible(False)
    p_legend.get_yaxis().set_visible(False)

    min_temp = np.min([res[unit]['temp_range'] for unit in units])
    max_temp = np.max([res[unit]['temp_range'] for unit in units])
    norm = mpl.colors.Normalize(vmin=min_temp, vmax=max_temp)

    for unit in units:

        flevels = res[unit]['flevels']
        p_cv.plot(flevels, res[unit]['cvs'], '*-')
        p_rate.plot(flevels, res[unit]['rates'], '*-')
        p_legend.plot(0, 0, '.-', label=path.basename(unit))

        isi_len = len(res[unit]['temp_t'])
        x_range = np.array(range(isi_len)) / float(isi_len-1)
        for i, temp_t in enumerate(res[unit]['temp_t']):
            c = cm.jet((temp_t - min_temp) / (max_temp - min_temp), 1)
            pp = p_temp.plot(x_range[i], res[unit]['temp_isis'][i], '.-', color=c)
    axins1 = inset_axes(p_temp, width="50%", height="5%", loc=1)
    mpl.colorbar.ColorbarBase(axins1, norm=norm, cmap=cm.jet,
                              orientation="horizontal",
                              ticks=[round(min_temp, 2)+0.01, round(max_temp, 2)])
    p_legend.legend()
    plt.savefig(path.join(out_folder, 'rates_vs_force_%s.%s' % (subsample, plot_format)))

    plt.figure()
    # take only the last 6 mechanical stimulations because this is what we 
    # have for all units in common. some also have 7 or 8 force levels
    rates = np.array([res[unit]['rates'][-6:] for unit in units])
    flevels = np.array([res[unit]['flevels'][-6:] for unit in units])
    plt.plot(np.mean(flevels, axis=0), np.mean(rates, axis=0))
    plt.title('mean force level vs. mean rate')
    plt.xlabel('mean force')
    plt.ylabel('mean rate')
    plt.savefig(path.join(out_folder, 'rates_vs_force_mean_%s.%s' % (subsample, plot_format)))

plt.show()
