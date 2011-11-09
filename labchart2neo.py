
import glob
import logging
from nexioplus import NexIOplus
import pylab as plt
import numpy as np
import quantities as pq
import math

# logger setup
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger()

data_folder = '/Users/dedan/projects/fu/data/'
downsample_factor = 100

nexlist = glob.glob(data_folder + '*.nex')

# actual data processing
for i, nexname in enumerate(nexlist):

    logger.info('read in: %s' % nexname)
    block = NexIOplus(filename=nexname, downsample=downsample_factor).read()

    plt.figure()
    for j, segments in enumerate(block.segments):

        plt.subplot(3, 3, j+1)
        rate = segments.analogsignals[0].sampling_rate
        length = len(segments.analogsignals[0])
        start = segments.analogsignals[0].t_start * pq.s

        for signal in segments.analogsignals:
            plt.plot(signal)

        train = segments.spiketrains[0]
        if np.any(train):
            tmp = np.zeros(length)
            isi = np.diff(train)
            print "CV: %f" % (np.std(isi)/np.mean(isi))

            for spike in train:
                tmp[int(math.floor((spike-start)*rate))] = 1

            plt.plot(tmp)


    plt.savefig("fig%d.pdf" % (i+1))

