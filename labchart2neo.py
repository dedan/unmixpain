#!/usr/bin/env python
# encoding: utf-8

"""
labchart2neo.py

We get the data for this project as labcharts .adicht files. They contain
4 channels:

1. electrode recording
2. the signal for the stepper motor (mechanical stimulation)
3. the mechanical stimulation signal
4. the heat stimulation signal

The files we get are already sorted for a unit and the only thing we are
using labchart for is to export the spiketrain of the first channel via

    Spike Histogram -> Export to NeuroExplorer

and the remaining 3 channels via

    File -> Export to a mat file

This script here is used to read in the `.nex` and the `.mat` file and
to convert it into the neo-datastructure:
http://packages.python.org/neo/classes.html

For this I wrote a class that inherits from NeuroExplorerIO which only reads
in the next files and extented its functionality to add the information from
the mat file

Created by Stephan Gabler on 2011-10-20.
"""

# about the nexfile
#
# when I read in a nexfile, I get a segment.
# the segment contains events. those events correspond to the blue lines
# in the labchart view of the file. the blue lines is when recording was
# switched on and off
# TODO is this right?

import sys
import os
import glob
import logging
import numpy as np
import scipy.io as sio
from neo.io.neuroexplorerio import NeuroExplorerIO
from neo.core import AnalogSignal

# logger setup
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%m-%d %H:%M')
logger = logging.getLogger()


class NexIOplus(NeuroExplorerIO):
    """Read from .nex file and add the info from .mat file

        this class inherits from the neo nexfile reader and adds
        functionality to read analog signals from a mat file and to add
        them to the neo data structure from the nex file
    """

    def __init__(self, filename=None):
        """init and check whether matfile available"""
        super(NexIOplus, self).__init__(filename)
        self.matname = self.filename[:-3] + 'mat'
        if not os.path.exists(self.matname):
            raise IOError('corresponding mat file not found')

    def read(self):
        """read the nex file and add the analog signals from matfile"""
        seg = super(NexIOplus, self).read_segment()
        mat = sio.loadmat(matname, squeeze_me=True)
        data_length = len(mat['data'])
        n_channels = len(mat['titles'])
        data = np.reshape(mat['data'], (n_channels, data_length/n_channels))
        channel_names = ['stepper', 'mechanic', 'heat']
        for i in range(n_channels):
            ansig = AnalogSignal(signal=np.array(data[i]),
                                 name=channel_names[i],
                                 units=mat['unittext'].item(),
                                 sampling_rate=mat['samplerate'][i, 0],
                                 t_start=0)
            seg.analogsignals.append(ansig)
        return seg


data_folder = '/Users/dedan/projects/fu/data/'

nexlist = glob.glob(data_folder + '*.nex')

# actual data processing
for nexname in nexlist:
    seg = NexIOplus(nexname)

    if len(seg.spiketrains) > 1:
        logger.error('more then one unit was defined in nexfile')
