#!/usr/bin/env python3

import unittest
from unittest import mock
import os
import sys

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import numpy as np

from arts_tracking_beams import create_tracking_beam
from arts_tracking_beams.constants import NTAB


class TestCreateTrackingBeam(unittest.TestCase):

    def setUp(self):
        # create some test files
        files = [f'ARTS00000_CB00_TAB{tab:02d}.fits' for tab in range(NTAB)]
        # add one file of different CB
        files.append(f'ARTS00000_CB01_TAB00.fits')
        for f in files:
            open(f, 'w').close()
        self.files = files

    def tearDown(self):
        # remove the test files
        for f in self.files:
            os.remove(f)

    def test_get_input_path(self):
        func = create_tracking_beam.get_input_path
        with self.assertRaises(SystemExit):
            # no files in some random folder
            func(input_folder='/')
        with self.assertRaises(SystemExit):
            # empty files were created in current folder. Checking all should find too many
            func(input_folder='.')
        with self.assertRaises(SystemExit):
            # search for wrong taskid and CB, should not find any files
            func(input_folder='.', taskid='12345678', cb=5)
        # specify correct cb, should find correct number of files
        func(input_folder='.', cb=0)

    def test_get_source_coordinates(self):
        func = create_tracking_beam.get_source_coordinates
        # known name
        func(name='PSR B0531+21')
        # ra and dec
        func(ra='01:00:00', dec='05:00:00')
        with self.assertRaises(SystemExit):
            # providing ra but not dec
            func(ra='01:00:00')
        with self.assertRaises(SystemExit):
            # providing unknown name
            func(name='PSR foo')

    def test_required_tb_resolution(self):
        func = create_tracking_beam.required_tb_resolution
        # create a pointing
        ra = 180 * u.deg
        dec = 30 * u.deg
        pointing = SkyCoord(ra, dec)
        t = Time('2020-01-01T00:00:00')
        # two durations: one with a single HA, one with both positive and negative HA
        durations = [1 * u.s, 12 * u.hour]
        for d in durations:
            func(pointing, t, duration=d, fhi=1520 * u.MHz)

    def test_get_ncpu(self):
        func = create_tracking_beam.get_ncpu
        # no ncpu given
        self.assertEqual(func(), os.cpu_count())
        # too many cpus
        self.assertEqual(func(9999), os.cpu_count())
        # valid amount
        self.assertEqual(func(1), 1)

    def test_get_subint_chunks(self):
        func = create_tracking_beam.get_subint_chunks
        # nsubint smaller than chunksize should not change values
        start = 12500
        nsubint = 100
        chunksize = 1000
        values = func(start, nsubint, chunksize)
        self.assertTrue(len(values) == 1)
        self.assertEqual(values[0][0], start)
        self.assertEqual(values[0][1], nsubint)

        # nsubint integer multiple of chunksize
        n = 2
        nsubint = n * chunksize
        values = func(start, nsubint, chunksize)
        self.assertTrue(len(values) == n)

        # non-integer multiple
        nsubint = int(2.5 * chunksize)
        values = func(start, nsubint, chunksize)
        self.assertTrue(len(values) == int(n + 1))

        # the sum of all nsubint values must equal original value
        self.assertTrue(nsubint == np.sum(values[:, 1]))

    # mock call to main
    @mock.patch('arts_tracking_beams.create_tracking_beam.main')
    def test_main_with_args(self, mocked_main):
        func = create_tracking_beam.main_with_args
        with self.assertRaises(SystemExit):
            # without arguments fails
            func()
        # call with only required argument
        folder = 'foo'
        create_tracking_beam.sys.argv = [sys.argv[0], '--input_folder', folder]
        func()
        # check argument to main is set properly
        args, kwargs = mocked_main.call_args
        self.assertTrue(hasattr(args[0], 'input_folder'))
        self.assertTrue(args[0].input_folder == folder)

    # ToDo: test main


if __name__ == '__main__':
    unittest.main()
