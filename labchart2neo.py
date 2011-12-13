
import glob
import logging
from nexioplus import NexIOplus
import pylab as plt
import numpy as np
import quantities as pq
from math import floor
from os import path
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

# logger setup
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger()

data_folder = '/Users/dedan/projects/fu/data/'
out_folder = path.join(data_folder, 'out')

# flags for what to plot
plot_stepper = False
plot_section_markers = False

# downsample to save working memory, but !!! if the factor is two high
# the stepper activity detection fails for small force levels.
downsample_factor = 10

# threshold to detect activity of the stepper-motor
stepper_thresh = 0.02

# width of the kernel used for smoothing the stepper signal
win_width = 100


def extract_onsets(signal, threshold, win_width):
    '''find onsets in stepper signal (stupid but workin version)

    Simple looks whether a certain threshold is exceeded. Signal is first
    convolved with a rectangular kernel to be less susceptible to small
    spikes.
    '''
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


def cv(train):
    '''compute the coefficient of variation'''
    isi = np.diff(win)
    return np.std(isi)/np.mean(isi)


nexlist = glob.glob(data_folder + '*.nex')

res = {}
wins = []

# for nexname in [nexlist[2]]:
for i, nexname in enumerate(nexlist):

    print ' '
    logger.info('read in: %s' % nexname)
    block = NexIOplus(filename=nexname, downsample=downsample_factor).read()

    res[nexname] = {'flevels': [], 'cvs': [], 'rates': [], 'isis': []}

    plt.figure()
    gs = gridspec.GridSpec(2, 2)
    gs.update(hspace=0.5)
    plt.subplot(gs[0,:])
    plt.title(path.basename(nexname))
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

        x_range = range(int(floor(start*rate)), int(floor(start*rate))+length)
        plt.plot(x_range, force, 'g')
        plt.plot(x_range, temp, 'r')
        if plot_stepper:
            plt.plot(x_range, stepper, 'b')

        # extract windows in which stimulation took place
        onoffs = extract_onsets(stepper, stepper_thresh, win_width)
        for x1, x2 in onoffs:
            flevel = np.max(force[x1:x2]) - np.min(force)
            res[nexname]['flevels'].append(flevel)
            if plot_section_markers:
                start_r = start * rate
                l = plt.plot([x1 + start_r, x2 + start_r], [0, 0], 'v')
                l[0].set_markersize(10)

        # plot the spiketrain of segment
        if np.any(train):
            tmp = np.zeros(length)
            for spike in train:
                tmp[int(floor((spike - start) * rate))] = 1
            plt.plot(x_range, tmp, 'k', linewidth=0.3)

        # extract spiketrain during stimulation
        onoffs_time = [(x1/rate, x2/rate) for x1, x2 in onoffs]
        for x1, x2 in onoffs_time:
            win = train[(train > x1+start) & (train < x2+start)]

            wins.append(win)
            res[nexname]['isis'].append(np.diff(win))
            res[nexname]['cvs'].append(cv(win))
            res[nexname]['rates'].append(len(win) / (x2 - x1))

    # annotate x-axis
    ticks = plt.xticks()
    plt.xticks(ticks[0], ticks[0]/rate, rotation=25)

    # plot the ISIs over time (x-axis normalized!)
    plt.subplot(gs[1, 0])
    for i, isi in enumerate(res[nexname]['isis']):
        if len(isi) > 3:
            bla = res[nexname]['flevels'][i] / np.max(res[nexname]['flevels'])
            c = cm.jet(bla, 1)
            plt.plot(np.array(range(len(isi))) / float(len(isi)-1), isi, '.-', color=c)
    plt.title('relative time vs. ISI')
    plt.xlabel('normalized time')
    plt.ylabel('ISI')

    # plot ISIs over temperature
    # TODO should I extract the temperature stimulation onset?
    plt.subplot(gs[1, 1])
    # find the spikes that occured during temp stimulation
    temp_spikes = train[train > temp.t_start]
    # extract the temperature for each of the spikes
    idx = [int(floor((spike - start) * rate)) for spike in temp_spikes]
    temp_t = temp[idx]
    plt.plot(temp_t[:-1], np.diff(temp_spikes), '.')
    plt.title('temperature vs. ISI')
    plt.xlabel('temperature')

    plt.savefig(path.join(out_folder, 'fig_%s.png') % path.basename(nexname))
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
    flevels = result['flevels'] / np.max(result['flevels'])
    p1.plot(flevels, result['cvs'], '*-')
    p2.plot(flevels, result['rates'], '*-')
    p3.plot(0, 0, '.-', label=path.basename(key))

p3.legend()
plt.savefig(path.join(out_folder, 'rates_vs_force_new.png'))
plt.show()
