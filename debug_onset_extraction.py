#!/usr/bin/env python
# encoding: utf-8

import pickle, json
import glob
import logging
import nexioplus
import pylab as plt
import numpy as np
import quantities as pq
from math import floor
from os import path
import matplotlib.gridspec as gridspec
reload(nexioplus)

# logger setup
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger()

data_folder = '/Users/dedan/projects/fu/data/unmixpain/'
data_folder = '/Users/dedan/projects/fu/unmixpain/tests/test_data/'
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


nexlist = glob.glob(data_folder + '*.nex')

res = {}
rate_estimates = {}
wins = []

# for i, nexname in enumerate(nexlist):
for i, nexname in enumerate([nexlist[0]]):

    print ' '
    logger.info('read in: %s' % nexname)
    block = nexioplus.NexIOplus(filename=nexname, downsample=downsample_factor).read()
    nex_base = path.basename(nexname)[:-4]

    fig = plt.figure()

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

        # plot the analog signals
        x_range = np.arange(start*rate, start*rate+length)

        plt.plot(x_range, force, 'g')
        plt.plot(x_range, stepper, 'b')

        # compute features in stimulation windows
        for epoch in segment.epochs:
            x1, x2 = epoch.time, epoch.time + epoch.duration
            print 'times: ', x1, x2

            plt.plot([x1+start*rate, x2+start*rate], [0, 0], '*r')
plt.show()
