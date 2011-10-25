
import glob
import logging
from nexioplus import NexIOplus
import pylab as plt

# logger setup
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger()

data_folder = '/Users/dedan/projects/fu/data/'
downsample_factor = 1000

nexlist = glob.glob(data_folder + '*.nex')
plt.figure()

# actual data processing
for i, nexname in enumerate(nexlist):
    plt.subplot(len(nexlist), 1, i)

    print 'read in: %s' % nexname
    seg = NexIOplus(filename=nexname, downsample=downsample_factor).read()
    if len(seg.spiketrains) > 1:
        logger.error('more then one unit was defined in nexfile')

    for signal in seg.analogsignals:
        plt.plot(signal)
plt.show()


