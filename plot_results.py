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
out_folder = path.join(data_folder, 'out')
res = pickle.load(open(path.join(out_folder, 'results.pickle')))

subsamples = {}
# TODO: '2010-09-10-u3' falls Du es einfach auf die gleiche Base-line setzen kannst, wie die anderen drei,
# dann sieht es so aus, als ob die Amplitude(n) den anderen entsprechen müsste(n).
# Falls aber nicht einfach, dann eben nicht machen..)
subsamples['sub1'] = ['2010-03-19-u1', '2010-03-19-u2', '2010-03-19-u4', '2010-09-10-u3']
subsamples['sub2'] = ['2010-12-02-u1', '2010-12-02-u2']
subsamples['sub3'] = ['2010-12-13-u3', '2010-12-13-u4', '2010-12-14-u4', '2010-12-14-u5']


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
    for i, temp_t in enumerate(result['temp_t']):
        min_temp = result['temp_range'][0]
        max_temp = result['temp_range'][1]
        c = cm.jet((temp_t - min_temp) / (max_temp - min_temp), 1)
        pp = p_res_temp.plot(i, result['temp_isis'][i], '.-', color=c)
    axins1 = inset_axes(p_res_temp, width="50%", height="5%", loc=1)
    norm = mpl.colors.Normalize(vmin=min_temp, vmax=max_temp)
    mpl.colorbar.ColorbarBase(axins1, norm=norm, cmap=cm.jet,
                              orientation="horizontal",
                              ticks=[round(min_temp, 2)+0.01, round(max_temp, 2)])
    p_res_temp.set_title('ISI over time \n (temperature coded by color)')
    p_res_temp.set_xlabel('time')
    p_res_temp.set_ylabel('ISI')
    fig_res.savefig(path.join(out_folder, 'fig_res_%s.png') % key)

plt.show()


# plot CV and firing rate in relation to stimulation force level
fig = plt.figure()
gs = gridspec.GridSpec(2, 2)
p1 = fig.add_subplot(gs[0, 0])
p1.set_title('CV against force level')
p2 = plt.subplot(gs[0, 1])
p2.set_title('rate against force level')
p3 = plt.subplot(gs[1, :])

for key, result in res.items():
    flevels = result['flevels'] # / np.max(result['flevels'])
    p1.plot(flevels, result['cvs'], '*-')
    p2.plot(flevels, result['rates'], '*-')
    p3.plot(0, 0, '.-', label=path.basename(key))

    p3.legend()
    plt.savefig(path.join(out_folder, 'rates_vs_force_new.png'))
plt.show()
