#!/usr/bin/env python3

import unittest
from unittest import mock
import os

from arts_tracking_beams import ARTSFITSReader
from arts_tracking_beams.constants import NTAB


class TestFitsReader(unittest.TestCase):

    @mock.patch(f'{__name__}.ARTSFITSReader._get_fits_params', return_value=('info', 'hdr_p', 'hdr_s'))
    @mock.patch(f'{__name__}.ARTSFITSReader._verify_files')
    def setUp(self, *mocks):
        # get path to current folder
        root = os.path.dirname(os.path.abspath(__file__))
        self.reader = ARTSFITSReader(f'{root}/data/ARTS20200101_CB00_TAB{{tab:02d}}.fits')

    def tearDown(self):
        # ensure files are closed
        self.reader.close_files()

    def test_verify_files(self):
        # run on real files
        self.reader._verify_files()
        # run on non-existing file, this only logs a warning
        fnames = self.reader.fnames
        self.reader.fnames = ['foo.fits']
        self.reader._verify_files()
        # put back the original filenames
        self.reader.fnames = fnames

    def test_get_fits_params(self):
        info, hdr_p, hdr_s = self.reader._get_fits_params()
        # test some basic parameters from the test files
        self.assertTrue(info['nsubint'] == 1)
        self.assertTrue(info['nchan'] == 384)
        self.assertTrue(hdr_p['RA'] == '12:00:00.0')
        self.assertTrue(hdr_p['SCANLEN'] == 1.024)
        self.assertTrue(hdr_s['NAXIS2'] == 1)
        self.assertTrue(hdr_s['TFORM17'] == '24000B')

    def test_open_close_files(self):
        # running open or close twice gives a warning, test that too
        self.reader.open_files()
        self.reader.open_files()
        self.reader.close_files()
        self.reader.close_files()

    def test_read_tabs_to_buffer(self):
        self.reader.info, self.reader.hdr_p, self.reader.hdr_s = self.reader._get_fits_params()
        self.reader.open_files()
        # set TAB00 as the output buffer
        out = self.reader.file_handles[0]
        # read one subint from the start
        start = 0
        nsubint = 1

        with self.assertRaises(SystemExit):
            # nchan not divisible by nsub
            tabs = range(7)
            self.reader.read_tabs_to_buffer(start, nsubint, tabs, out['SUBINT'])
        with self.assertRaises(SystemExit):
            # not a multiple of 8 channels per subband
            tabs = range(384)
            self.reader.read_tabs_to_buffer(start, nsubint, tabs, out['SUBINT'])

        # read the same tab to each subband
        tabs = 10
        self.reader.read_tabs_to_buffer(start, nsubint, tabs, out['SUBINT'])
        # input has weights equal to TAB number, check output is the same
        data = out['SUBINT'].data['DAT_WTS'][0]
        self.assertListEqual(
            data.astype(int).tolist(), [tabs] * len(data))

        # read different TAB for each subband
        # first put all weights to zero
        out['SUBINT'].data['DAT_WTS'][:] = 0
        tabs = range(NTAB)
        self.reader.read_tabs_to_buffer(start, nsubint, tabs, out['SUBINT'])
        # check each TAB is used, reversed as df < 0
        data = out['SUBINT'].data['DAT_WTS'][0]
        nchan_per_tab = int(self.reader.info['nchan'] / NTAB)  # 32 for nchan = 384 and ntab = 12
        for i, tab in enumerate(reversed(tabs)):
            weights = data[i * nchan_per_tab:(i + 1) * nchan_per_tab].astype(int).tolist()
            self.assertListEqual(weights, [tab] * len(weights))


if __name__ == '__main__':
    unittest.main()
