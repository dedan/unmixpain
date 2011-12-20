#!/usr/bin/env python
# encoding: utf-8
"""
test_nexioplus.py

Created by Stephan Gabler on 2011-12-13.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import unittest
import os
import numpy as np
from math import floor
from unmixpain.nexioplus import NexIOplus
import pylab as plt


class test_nexioplus(unittest.TestCase):

    def setUp(self):
        basepath = os.path.join(os.path.dirname(__file__), 'test_data')
        self.fname = os.path.join(basepath, '2010-09-10-u3_test.nex')
        self.fname_without = os.path.join(basepath, 'without_mat.nex')
        self.block = NexIOplus(filename=self.fname).read()

    def test_filename_stored(self):
        """original filename should be stored in the block"""
        self.assertEqual(self.block.file_origin, self.fname)

    def test_fail_when_no_matlab_file_exists(self):
        """should raise an exception if corresponding mat file not found"""
        with self.assertRaises(IOError):
            NexIOplus(filename=self.fname_without)

    def test_segment_boring_detected_and_neglected(self):
        """should not read the segments without any stimulation"""
        self.assertEqual(len(self.block.segments), 3)

    def test_segment_tempchange_detected(self):
        """segment with change in temperature should be labeled accordingly"""
        self.assertIn('temp', self.block.segments[-1].annotations)

    def test_segment_mechanical_stim_detected(self):
        """segment with mechanical stimulation should be labeled accordingly"""
        self.assertIn('mechanical', self.block.segments[0].annotations)

    def test_mechanical_stim_markers_detected(self):
        """detect the on- and offsets of mechanical stimulation"""
        self.assertEqual(len(self.block.segments[0].epochs), 1)
        self.assertEqual(len(self.block.segments[1].epochs), 2)

    def test_mechanical_stim_markers_detected_correctly(self):
        """check whether the computed epoch markers are in the right place

        correct values from visual inspection
        """
        epoch = self.block.segments[0].epochs[0]
        self.assert_(39980 < epoch.time < 40000)
        self.assert_(240120 < epoch.time + epoch.duration < 240140)

    def test_spike_alignment(self):
        """test the spike alignment

            check whether the spike time from the .nex file fits the
            spike in the raw analog signal where the .nex file was
            extracted from. The peak of the spike should not differ by more
            than half the spikewidth (1ms), so it should not differ by more
            than 0.5 ms which is 10 indeces at 20 KHz
        """
        spike_signal = self.block.segments[1].analogsignals[0]
        start = spike_signal.t_start
        rate = spike_signal.sampling_rate
        for spike in self.block.segments[1].spiketrains[0]:
            spike_pos = int(floor((spike - start) * rate))
            # extract 2 ms around spike (40 points @ 20 KHz)
            window_around_spike = spike_signal[spike_pos-20:spike_pos+20]
            spike_pos_signal = np.argmax(window_around_spike)
            self.assert_(10 < spike_pos_signal < 30)

    def test_downsample(self):
        """make sure that downsampling changes the size correctly"""
        downsample_factor = 10
        original_length = len(self.block.segments[1].analogsignals[0])
        down_block = NexIOplus(filename=self.fname,
                               downsample=downsample_factor).read()
        down_length = len(down_block.segments[1].analogsignals[0])
        self.assertEqual(original_length / downsample_factor, down_length)

if __name__ == '__main__':
    unittest.main()