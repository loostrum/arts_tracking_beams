#!/usr/bin/env python3

import unittest
import astropy.units as u
from astropy.time import Time
import numpy as np

from arts_tracking_beams import TrackingBeam
from arts_tracking_beams.constants import NTAB


class TestTrackingBeam(unittest.TestCase):

    def test_source_pointing(self):
        ra = 180 * u.deg
        dec = 30 * u.deg
        nsub = 48
        times = Time('2020-01-01T12:00:00') + np.arange(24) * u.hour
        TB = TrackingBeam(ra, dec, ra, dec, nsub=nsub)
        expected_tabs = [0] * nsub
        for t in times:
            tabs = TB.run(t)
            self.assertListEqual(tabs.tolist(), expected_tabs)

    def test_all_tabs(self):
        ra = 180 * u.deg
        dec = 30 * u.deg
        ra_src = ra
        dec_src = dec + 5 * u.arcmin
        nsub = 1
        times = Time('2020-01-01T12:00:00') + np.linspace(0, 24, 100) * u.hour
        TB = TrackingBeam(ra, dec, ra_src, dec_src, nsub=nsub)
        tabs = [TB.run(t)[0] for t in times]
        for tab in range(NTAB):
            self.assertIn(tab, tabs)

    def test_freqs(self):
        ra = 180 * u.deg
        dec = 30 * u.deg
        ra_src = ra + 3 * u.deg
        dec_src = dec + 3 * u.deg
        nsub = 384
        t = Time('2020-01-01T12:00:00')
        # init tracking beam
        TB = TrackingBeam(ra, dec, ra_src, dec_src, nsub=nsub)
        tabs = TB.run(t)
        for tab in range(NTAB):
            self.assertIn(tab, tabs)

    def test_time_array(self):
        ra = 180 * u.deg
        dec = 30 * u.deg
        ra_src = ra + 3 * u.deg
        dec_src = dec + 3 * u.deg
        nsub = 384
        ntime = 10
        t = Time('2020-01-01T12:00:00') + np.arange(ntime) * u.s
        # init tracking beam
        TB = TrackingBeam(ra, dec, ra_src, dec_src, nsub=nsub)
        tabs = TB.run(t)
        self.assertTupleEqual(tabs.shape, (ntime, nsub))
        # again with list instead of array of Times
        t = t.isot.tolist()
        tabs = TB.run(t)
        self.assertTupleEqual(tabs.shape, (ntime, nsub))


if __name__ == '__main__':
    unittest.main()
