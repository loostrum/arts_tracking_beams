#!usr/bin/env python3

import os
import logging
from shutil import copyfile

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import tqdm

from arts_tracking_beams import ARTSFITSReader

logging.basicConfig(format='%(asctime)s.%(levelname)s.%(name)s: %(message)s', level=logging.DEBUG)


def main():
    # load the TB info
    tb_data = np.loadtxt('TB.txt', dtype=int)
    #tb_res = 60  # subints
    tb_res = 1  # subints
    # first column is time step, remainder is TAB indices
    time_steps = tb_data[:, 0]
    tab_indices = tb_data[:, 1:]
    # get number of subbands
    nsub = len(tab_indices[0])

    src = 'PSR J2215+1538'
    src_coord = SkyCoord.from_name(src)
    ra, dec = src_coord.to_string('hmsdms', sep=':', precision=1).split(' ')

    # input/output files
    # input_file = 'X:/PSRJ2215+1538_202000409_CB14/fits/ARTS200409008_CB14_TAB{tab}.fits'
    #input_file = 'F:/Leon/Desktop/python/tracking_beams_testing/data/fits/ARTS200409008_CB14_TAB{tab}.fits'
    #output_file = 'F:/Leon/Desktop/python/arts_tracking_beams/TB.fits'
    input_file = 'input/ARTS200826172_CB26_TAB{tab}.fits'
    output_file = 'TB.fits'
    try:
        os.remove(output_file)
    except FileNotFoundError:
        pass

    # init FITS reader
    reader = ARTSFITSReader(input_file)
    # open the input files
    reader.open_files()

    # use TAB00 file to generate output
    output_handle = reader.file_handles[0]

    #n = 172
    # n = 10
    # update output header keys
    output_handle['PRIMARY'].header['RA'] = ra
    output_handle['PRIMARY'].header['DEC'] = dec
    output_handle['PRIMARY'].header['SRC_NAME'] = src
    #output_handle['SUBINT'].header['NAXIS2'] = n * tb_res
    output_handle['SUBINT'].header['NAXIS2'] = reader.info['nsubint']
    # update truncated scanlen
    output_handle['PRIMARY'].header['SCANLEN'] = reader.info['duration'].to('second').value

    # process subints
    nsubint_to_read = tb_res
    nsubint = reader.info['nsubint']
    for subint in tqdm.tqdm(range(nsubint)):
    #for sub_subint in tqdm.tqdm(range(n)):
        #subint = sub_subint * nsubint_to_read
        # get index in TB list using floor division
        tb_index = subint // tb_res
        # extract tabs to read at this time step
        try:
            tabs = tab_indices[tb_index]
        except IndexError:
            print("USING LAST INDEX")
            tabs = tab_indices[-1]
        # extract data
        reader.read_data(subint, nsubint_to_read, tabs, output_handle['SUBINT'])

    # write the new file
    output_handle.writeto(output_file)
    output_handle.close()
    # close the input files
    reader.close_files()
