#!/usr/bin/env python
# encoding: utf-8

import pickle, json
import glob
import logging
from nexioplus import NexIOplus
import pylab as plt
import numpy as np
import quantities as pq
from math import floor
from os import path
import matplotlib.gridspec as gridspec

# logger setup
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger()

data_folder = '/Users/dedan/projects/fu/data/unmixpain/'
plot_format = 'pdf'
out_folder = path.join(data_folder, 'out')

# downsample to save working memory, but !!! if the factor is two high
# the stepper activity detection fails for small force levels.
downsample_factor = 10

# threshold to detect activity of the stepper-motor
stepper_thresh = 0.02

# width of the kernel used for smoothing the stepper signal
win_width = 100

# width of the kernel for the firing rate estimate in seconds
kernel_delta_t = 0.1

def cv(train):
    '''compute the coefficient of variation'''
    isi = np.diff(win)
    return np.std(isi)/np.mean(isi)


nexlist = glob.glob(data_folder + '*.nex')

res = {}
rate_estimates = {}
wins = []

for i, nexname in enumerate(nexlist):
# for i, nexname in enumerate([nexlist[1]]):

    print ' '
    logger.info('read in: %s' % nexname)
    block = NexIOplus(filename=nexname, downsample=downsample_factor).read()
    nex_base = path.basename(nexname)[:-4]

    res[nex_base] = {'flevels': [], 'cvs': [], 'rates': [], 'isis': []}

    # prepare all the figures
    fig_signal = plt.figure()
    gs = gridspec.GridSpec(4, 1)
    gs.update(hspace=0.5)
    p_temp = fig_signal.add_subplot(gs[0, 0])
    p_temp.set_title(path.basename(nexname))
    p_mech = fig_signal.add_subplot(gs[1, 0])
    p_spik = fig_signal.add_subplot(gs[2, 0])
    p_esti = fig_signal.add_subplot(gs[3, 0])

    rate_estimates[nex_base] = {'data': []}
    for segment in block.segments:

        train = segment.spiketrains[0]
        if len(segment.analogsignals) == 3:
            spikes = []
            stepper = segment.analogsignals[0]
            force = segment.analogsignals[1]
            temp = segment.analogsignals[2]
        else:
            spikes = segment.analogsignals[0]
            stepper = segment.analogsignals[1]
            force = segment.analogsignals[2]
            temp = segment.analogsignals[3]

        rate = stepper.sampling_rate
        length = len(stepper)
        start = stepper.t_start

        res[nex_base]['train'] = train.magnitude
        if 'temp' in segment.annotations:
            # find the spikes that occured during temp stimulation
            temp_spikes = train[train > temp.t_start]
            # extract the temperature for each of the spikes
            idx = [int(floor((spike - start) * rate)) for spike in temp_spikes]
            res[nex_base]['temp_t'] = temp[idx[1:]].magnitude
            res[nex_base]['temp_range'] = [np.min(temp).magnitude, np.max(temp).magnitude]
            res[nex_base]['temp_isis'] = np.diff(temp_spikes).magnitude

        # plot the analog signals
        x_range = range(int(floor(start*rate)), int(floor(start*rate))+length)
        p_mech.plot(x_range, force, 'g')
        p_temp.plot(x_range, temp, 'r')

        # plot the spiketrain
        if np.any(train):
            tmp = np.zeros(length)
            for spike in train:
                tmp[int(floor((spike - start) * rate))] = 1
            p_spik.plot(x_range, tmp, 'k', linewidth=0.3)

            kernel = np.ones(rate * kernel_delta_t) / kernel_delta_t
            rate_estimate = np.convolve(tmp, kernel, mode='same')
            p_esti.plot(x_range, rate_estimate, 'b')

        # compute features in stimulation windows
        for epoch in segment.epochs:
            x1, x2 = epoch.time, epoch.time + epoch.duration
            x1_t, x2_t = (x1 / rate) + start, (x2 / rate) + start
            flevel = np.max(force[x1:x2]) - np.min(force)
            res[nex_base]['flevels'].append(flevel.magnitude.item())
            # extract spiketrain during stimulation
            win = train[(train > x1_t) & (train < x2_t)]

            wins.append(win)
            res[nex_base]['isis'].append(np.diff(win).magnitude)
            res[nex_base]['cvs'].append(cv(win).magnitude.item())
            res[nex_base]['rates'].append(float(len(win)) / (x2 - x1))
            rate_estimates[nex_base]['data'].append(rate_estimate[x1:x2])
            rate_estimates[nex_base]['start'] = x1
            rate_estimates[nex_base]['end'] = x2

    # annotate axis
    ticks = p_temp.get_xticks()
    p_temp.set_xticks(ticks)
    p_temp.set_xticklabels(ticks/rate, rotation=10)
    p_esti.set_xlabel('time')
    p_esti.set_ylabel('firing rate')
    p_esti.set_xticklabels([])
    p_mech.set_xticklabels([])
    p_spik.set_xticklabels([])
    p_spik.set_xlim(p_mech.get_xlim())
    p_esti.set_xlim(p_mech.get_xlim())
    fig_signal.savefig(path.join(out_folder, 'fig_signal_%s.%s') % 
                       (path.basename(nexname), plot_format))
    plt.show()

subsamples = json.load(open(path.join(data_folder, 'subsamples.json')))
for subsample, units in subsamples.items():
    p_mean_est = plt.figure()
    subs = []
    print subsample
    epoch_lengths = np.array([[len(epo) for epo in rate_estimates[unit]['data']][-6:]
                                for unit in units])
    min_epoch_lengths = np.min(epoch_lengths, axis=0)
    epochs = []
    for unit in units:
        tmp = []
        for i, epo in enumerate(rate_estimates[unit]['data'][-6:]):
            tmp.append(np.array(epo[-min_epoch_lengths[i]:]))
        epochs.append(tmp)
    blam = np.mean(epochs, axis=0)
    for i, mepo in enumerate(blam):
        sub = p_mean_est.add_subplot(1, 6, i+1)
        subs.append(sub)
        sub.plot(mepo)
    for i, sub in enumerate(subs):
        sub.set_ylim(subs[-1].get_ylim())
        sub.set_xticklabels([])
        if not i == 0:
            sub.set_yticklabels([])
    p_mean_est.suptitle('%s:  %s' % (subsample, ', '.join(units)))
    p_mean_est.savefig(path.join(out_folder, 'fig_meanest_%s.%s') 
                       % (subsample, plot_format))
    plt.show()
pickle.dump(res, open(path.join(out_folder, 'results.pickle'), 'w'))











