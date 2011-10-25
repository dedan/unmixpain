#!/usr/bin/env python
# encoding: utf-8

"""
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

This class is used to read in the `.nex` and the `.mat` file and
to convert it into the neo-datastructure:
http://packages.python.org/neo/classes.html

The read function of NexIOplus returns a segment:

* the spiketrains and analog signals can be found in the corresponding
  variables of the segment
* The events in the resulting segment correspond to the times when recording
  in labchart was switched on and off (blue lines)

Created by Stephan Gabler on 2011-10-20.
"""

import os
import numpy as np
import scipy.io as sio
from neo.io.neuroexplorerio import NeuroExplorerIO
from neo.core import AnalogSignal


class NexIOplus(NeuroExplorerIO):
    """Read from .nex file and add the info from .mat file

        this class inherits from the neo nexfile reader and adds
        functionality to read analog signals from a mat file and to add
        them to the neo data structure from the nex file
    """

    def __init__(self, filename=None, downsample=1):
        """init and check whether matfile available

        also provides optional downsampling because the mechanical and temp-
        stimulation have a much higher sampling rate then needed

        Keyword arguments:
        filename -- the file to read from
        downsample -- 1 by default, but if another integer number is given then
                      only each n-th value of the original signal is used
        """
        super(NexIOplus, self).__init__(filename)
        if not isinstance(downsample, int):
            raise ValueError('downsampling value has to be int')
        else:
            self.downsample = downsample
        self.matname = self.filename[:-3] + 'mat'
        if not os.path.exists(self.matname):
            raise IOError('corresponding mat file not found')

    def read(self):
        """read the nex file and add the analog signals from matfile

        all the data from different channels and blocks is found in a big
        data array and has to be extracted by the indexes in other variables
        """
        seg = super(NexIOplus, self).read_segment()

        channel_names = ['stepper', 'mechanic', 'heat']
        mat = sio.loadmat(self.matname, squeeze_me=True)
        n_channels, n_blocks = np.shape(mat['datastart'])

        for i in range(n_channels):
            rate = mat['samplerate'][i, 0] / self.downsample
            tmp = np.array(())
            for j in range(n_blocks):
                start = mat['datastart'][i,j] -1
                end = mat['dataend'][i,j]
                tmp = np.append(tmp, mat['data'][start:end])
            ansig = AnalogSignal(signal=tmp[::self.downsample],
                                 name=channel_names[i],
                                 units=mat['unittext'].item(),
                                 sampling_rate=rate,
                                 t_start=0)
            seg.analogsignals.append(ansig)
        return seg
