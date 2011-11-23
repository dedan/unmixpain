
import glob
import logging
from nexioplus import NexIOplus
import pylab as plt
import numpy as np
import quantities as pq
from math import floor
import os

# logger setup
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger()

data_folder = '/Users/dedan/projects/fu/data/'
downsample_factor = 10
boring_thresh = 1
stepper_thresh = 0.1
win_width = 100


def extract_onsets(signal, threshold, win_width):
    """find onsets in stepper signal"""
    s = np.convolve(abs(signal), np.ones(win_width) / win_width)
    x = 0
    onoffs = []
    for i, y in enumerate(s):
        if y > threshold:
            if x == 0:
                x = i
        elif not x == 0:
            onoffs.append((x, i))
            x = 0
    # sort them differently because my algo is so stupid
    onoffs = [(onoffs[i][1], onoffs[i+1][0]) for i in range(len(onoffs)-1)]
    # don't use the times in between
    return onoffs[::2]

def is_boring_segment(signals, threshold):
    """docstring for is_boring_segment"""
    no_temp_change = np.max(signals[2]) - np.min(signals[2]) < threshold
    no_stepper = np.max(signals[0]) - np.min(signals[0]) < threshold
    return no_temp_change and no_stepper

def cv(train):
    """compute the coefficient of variation"""
    isi = np.diff(win)
    return np.std(isi)/np.mean(isi)



nexlist = glob.glob(data_folder + '*.nex')

res = {}
wins = []

# actual data processing
# for nexname in [nexlist[2]]:
for i, nexname in enumerate(nexlist):

    print " "
    logger.info('read in: %s' % nexname)
    block = NexIOplus(filename=nexname, downsample=downsample_factor).read()

    res[nexname] = {"flevels": [], "cvs": []}

    plt.figure()
    for segment in block.segments:

        if is_boring_segment(segment.analogsignals, boring_thresh):
            continue

        stepper = segment.analogsignals[0]
        force = segment.analogsignals[1]
        temp = segment.analogsignals[2]

        rate = stepper.sampling_rate * 1/pq.s
        length = len(stepper)
        start = stepper.t_start * pq.s

        x_range = range(int(floor(start*rate)), int(floor(start*rate))+length)
        # plt.plot(x_range, stepper, 'b')
        plt.plot(x_range, force, 'g')
        plt.plot(x_range, temp, 'r')

        onoffs = extract_onsets(stepper, stepper_thresh, win_width)
        for x1, x2 in onoffs:
            flevel = np.max(force[x1:x2]) - np.min(force)
            res[nexname]["flevels"].append(flevel)
            l = plt.plot([x1 + start*rate, x2 + start*rate], [0, 0], 'v')
            l[0].set_markersize(10)

        train = segment.spiketrains[0]
        if np.any(train):
            tmp = np.zeros(length)
            for spike in train:
                tmp[int(floor((spike - start) * rate))] = 1
            plt.plot(x_range, tmp, 'k', linewidth=0.3)

        onoffs_time = [(x1/rate, x2/rate) for x1, x2 in onoffs]
        if onoffs_time:
            for x1, x2 in onoffs_time:
                win = train[(train > x1+start) & (train < x2+start)]
                wins.append(win)
                res[nexname]["cvs"].append(cv(win))
    ticks = plt.xticks()
    plt.xticks(ticks[0], ticks[0]/rate, rotation=25)
    plt.show()
    # plt.savefig("fig%d.pdf" % (i+1))

plt.figure()
plt.axes([0.1, 0.1, 0.5, 0.8])
plt.title("CV against force level")
for key, bla in res.items():
    plt.plot(bla['flevels'], bla['cvs'], '*', label=os.path.basename(key))
plt.legend(loc=(1, 0))
# plt.savefig("cv_vs_force.pdf")
plt.show()

# isi over spike-number
plt.figure()
for win in wins:
    plt.plot(np.diff(win), '.')
plt.show()

