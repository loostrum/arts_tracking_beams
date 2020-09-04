#!/usr/bin/env python3

import unittest
import astropy.units as u
import numpy as np

from arts_tracking_beams import tools


class TestTools(unittest.TestCase):

    def test_radec_hadec(self):
        ra0 = 180 * u.deg
        dec0 = 30 * u.deg
        # give time as string to test conversion to Time instance
        t = '2020-01-01T12:00:00'
        # convert RADEC to HADEC and back, check if values are equal
        ha, dec = tools.radec_to_hadec(ra0, dec0, t)
        ra1, dec1 = tools.hadec_to_radec(ha, dec, t)
        self.assertTrue(np.isclose(ra1.to(u.deg).value, ra0.to(u.deg).value), msg=f'{ra1} != {ra0}')
        self.assertTrue(np.isclose(dec1.to(u.deg).value, dec0.to(u.deg).value), msg=f'{dec1} != {dec0}')


if __name__ == '__main__':
    unittest.main()
