#!/usr/bin/env python3

import unittest
from unittest import mock
import os
import sys

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import numpy as np

from arts_tracking_beams import create_synthesised_beam
from arts_tracking_beams.constants import NTAB


class TestCreateSynthesisedBeam(unittest.TestCase):

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
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
        for f in ['ARTS_SOURCE.fits', 'ARTS20200101_CB00_SB35.fits']:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass

    @staticmethod
    def get_args(switched=False):
        root = os.path.abspath(os.path.dirname(__file__))
        args = {'input_folder': f'{root}/data',
                'taskid': None,
                'cb': None,
                'sb': 35,
                'output': 'out.fits',
                'overwrite': False,
                'source': 'SOURCE',
                'chunksize': 1000,
                'no_progress_bar': False,
                'verbose': False}
        args_switched = {'overwrite': True,
                         'output': None,
                         'verbose': True}
        if switched:
            args.update(args_switched)
        return args

    # mock call to main
    @mock.patch('arts_tracking_beams.create_synthesised_beam.main')
    def test_main_with_args(self, mocked_main):
        func = create_synthesised_beam.main_with_args
        with self.assertRaises(SystemExit):
            # without arguments fails
            func()
        # call with only required arguments
        folder = 'foo'
        sb = 35
        create_synthesised_beam.sys.argv = [sys.argv[0], '--input_folder', folder, '--sb', str(sb)]
        func()
        # check argument to main is set properly
        args, kwargs = mocked_main.call_args
        self.assertTrue(hasattr(args[0], 'input_folder'))
        self.assertTrue(args[0].input_folder == folder)
        self.assertTrue(hasattr(args[0], 'sb'))
        self.assertTrue(args[0].sb == sb)

    def test_main(self):
        args = self.get_args()
        create_synthesised_beam.main(args)
        # check if output files exist
        self.assertTrue(os.path.isfile(args['output']))

        with self.assertRaises(SystemExit):
            # output already exists so gives error
            create_synthesised_beam.main(args)

        # remove the output file
        os.remove(args['output'])

        # get another set of args with switched flags and enabled output overwrite
        args = self.get_args(switched=True)
        create_synthesised_beam.main(args)

        # output file name should now be based on source name
        self.assertTrue(os.path.isfile(f'ARTS_{args["source"]}.fits'))

        # remove source name, output is now based on taskid, CB, SB
        args['source'] = None
        create_synthesised_beam.main(args)
        self.assertTrue(os.path.isfile(f'ARTS20200101_CB00_SB{args["sb"]}.fits'))

        args['sb'] = 9999
        with self.assertRaises(SystemExit):
            # SB out of range gives error
            create_synthesised_beam.main(args)


if __name__ == '__main__':
    unittest.main()
