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
from neo.core import AnalogSignal, Block, Segment, SpikeTrain, Epoch
import quantities as pq


class NexIOplus(NeuroExplorerIO):
    """Read from .nex file and add the info from .mat file

        this class inherits from the neo nexfile reader and adds
        functionality to read analog signals from a mat file and to add
        them to the neo data structure from the nex file
    """

    def __init__(self,
                 filename=None,
                 downsample=1,
                 boring_threshold=1,
                 stepper_threshold=0.02,
                 win_width=100):
        """init and check whether matfile available

        also provides optional downsampling because the mechanical and temp-
        stimulation have a much higher sampling rate then needed

        Keyword arguments:
        filename -- the file to read from
        downsample -- 1 by default, but if another integer number is given then
                      only each n-th value of the original signal is used
        """
        super(NexIOplus, self).__init__(filename)
        self.f_down = float(downsample)
        self.matname = self.filename[:-3] + 'mat'
        self.boring_thresh = boring_threshold
        self.stepper_thresh = stepper_threshold
        self.win_width = win_width
        if not os.path.exists(self.matname):
            raise IOError('corresponding mat file not found')

    def read(self):
        """read the nex file and add the analog signals from matfile

        all the data from different channels and blocks is found in a big
        data array and has to be extracted by the indexes in other variables
        """

        nex_block = super(NexIOplus, self).read_segment()

        # create a new block to return in the end
        block = Block(name='mechanical and heat stimulation recording',
                      description=nex_block.annotations,
                      file_origin=self.filename)
        train = nex_block.spiketrains[0]

        # load data from matlab file
        mat = sio.loadmat(self.matname, squeeze_me=True)
        n_channels, n_segments = np.shape(mat['datastart'])

        # convert blocktimes to posix format (from stupid matlab convention)
        blockt_pos = (mat['blocktimes'] - 719529) * 86400.0
        blockt_pos = blockt_pos - blockt_pos[0]

        for segment in range(n_segments):

            seg = Segment(name=str(segment))
            for channel in range(n_channels):

                rate = mat['samplerate'][channel, segment] / self.f_down
                start = mat['datastart'][channel, segment] - 1
                end = mat['dataend'][channel, segment]
                tmp_sig = mat['data'][start:end][::self.f_down]

                # check for mechanical stimulation
                if (n_channels, channel) in [(3, 0), (4, 1)]:
                    onoffs = self._extract_onsets(tmp_sig)
                    for i, (x1, x2) in enumerate(onoffs):
                        seg.epochs.append(Epoch(x1, x2-x1, i+1))
                    if onoffs:
                        seg.annotate(mechanical=True)
                if (n_channels, channel) in [(3, 2), (4, 3)]:
                    if np.max(tmp_sig) - np.min(tmp_sig) > self.boring_thresh:
                        seg.annotate(temp=True)

                ansig = AnalogSignal(signal=tmp_sig,
                                     name=mat['titles'][channel],
                                     # TODO use unittextmap properly
                                     units=mat['unittext'].item(),
                                     sampling_rate=rate * (1/pq.s),
                                     t_start=blockt_pos[segment] * pq.s)
                seg.analogsignals.append(ansig)

            # ignore segments without heat or mechanical stimulation
            if not seg.annotations:
                continue

            # last segment has to be treated differently
            if segment + 1 < n_segments:
                t = train[(train > blockt_pos[segment]) &
                          (train < blockt_pos[segment+1])]
                end = blockt_pos[segment+1]
            else:
                t = train[train > blockt_pos[segment]]
                end = blockt_pos[segment] + len(ansig)/rate

            seg.spiketrains.append(SpikeTrain(times=t.magnitude,
                                   units=train.units,
                                   t_start=blockt_pos[segment],
                                   t_stop=end))
            block.segments.append(seg)
        return block

    def _extract_onsets(self, signal):
        '''find onsets in stepper signal (stupid but workin version)

        Simple looks whether a certain threshold is exceeded. Signal is first
        convolved with a rectangular kernel to be less susceptible to small
        spikes.
        '''
        kernel = np.ones(self.win_width) / self.win_width
        s = np.convolve(abs(signal), kernel, 'same')
        x = 0
        onoffs = []
        for i, y in enumerate(s):
            if y > self.stepper_thresh:
                if x == 0:
                    x = i
            elif not x == 0:
                onoffs.append((x, i))
                x = 0
        # remove half the kernel width because a shift is introduced by the convolution
        onoffs = [(o[0]+self.win_width/2, o[1]-self.win_width/2) for o in onoffs]
        # sort them differently because my algo is so stupid
        onoffs = [(onoffs[i][0], onoffs[i+1][1]) for i in range(len(onoffs)-1)]
        # don't use the times in between
        return onoffs[::2]
