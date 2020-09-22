#!/usr/bin/env python3

import os
import sys
import argparse
import logging

import numpy as np
import tqdm

from arts_tracking_beams import ARTSFITSReader
from arts_tracking_beams.constants import NSB, SB_TABLE
from arts_tracking_beams.create_tracking_beam import get_input_path, get_subint_chunks


def get_tabs_from_sb(sb):
    """
    Get the TAB indices corresponding to the given SB, in low-high frequency order

    :param int sb: Synthesised beam index
    :return: TAB indices (array)
    """
    if sb < 0 or sb >= NSB:
        logging.error(f'Synthesised beam index ({sb:02d}) out of range, must be between 0 and {NSB-1}')
        sys.exit(1)

    # load the SB table
    logging.debug(f'Loading SB table {SB_TABLE}')
    sb_mapping = np.loadtxt(SB_TABLE, dtype=int)

    # get the TABs of the requested sb
    try:
        tabs = sb_mapping[sb]
    except IndexError:  # pragma: no cover (because provided SB table is valid)
        logging.error(f'Could not get SB{sb:02d} from SB table with {len(sb_mapping)} entries; verify SB table.')
        sys.exit(1)

    logging.debug(f'SB{sb:02d} consists of TABs: {", ".join([f"{tab:02d}" for tab in tabs])}')
    return tabs


def main(args):
    """
    The main synthesised beam program.

    :param argparse.Namespace/dict args: arguments
    """
    # convert args to namespace if given as dict, so they can be accessed as attributes
    if isinstance(args, dict):
        args = argparse.Namespace(**args)

    # logger setup
    if args.verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO
    logging.basicConfig(format='%(asctime)s.%(levelname)s: %(message)s', level=level)

    # identify input files
    input_fits_path = get_input_path(args.input_folder, taskid=args.taskid, cb=args.cb)

    # get output file name
    if args.output is not None:
        output_file = args.output
    elif args.source is not None:
        output_file = f'ARTS_{args.source}.fits'
    else:
        # get taskid and CB index
        prefix, cb, _ = os.path.basename(input_fits_path).split('_')
        output_file = f'{prefix}_{cb}_SB{args.sb:02d}.fits'

    # check if output file already exists
    if os.path.isfile(output_file) and not args.overwrite:
        logging.error(f'Output file already exists: {output_file}. Use --overwrite to overwrite.')
        sys.exit(1)

    # get the TAB indices of the given SB
    tabs = get_tabs_from_sb(args.sb)

    # initialize the FITS reader
    fits_reader = ARTSFITSReader(input_fits_path)
    fits_reader.open_files()

    # write the output to TAB00
    # this does not change the TAB00 file on disk, but anything in the FITS data that is not changed
    # will have the TAB00 values
    output_handle = fits_reader.file_handles[0]

    # update output header keys
    if args.source is None:
        src_name = f'SB{args.sb:02d}'
    else:
        src_name = args.source
    output_handle['PRIMARY'].header['SRC_NAME'] = src_name
    # update scanlen, as its original value is truncated to 1 decimal place
    output_handle['PRIMARY'].header['SCANLEN'] = fits_reader.info['duration'].to('second').value
    output_handle['SUBINT'].header['NAXIS2'] = fits_reader.info['nsubint']

    # loop over the subints in steps of the max chunksize
    nsubint_tot = fits_reader.info['nsubint']
    for subint_start, nsubint in tqdm.tqdm(get_subint_chunks(0, nsubint_tot, args.chunksize),
                                           desc='Creating synthesised beam', disable=args.no_progress_bar):
        fits_reader.read_tabs_to_buffer(subint_start, nsubint, tabs, output_handle['SUBINT'])

    # create the output file
    output_handle.writeto(output_file, overwrite=args.overwrite)
    output_handle.close()
    # close the input files
    fits_reader.close_files()
    logging.info('Done')


def main_with_args():
    parser = argparse.ArgumentParser(description='Create a synthesised beam from ARTS tied-array beam data. '
                                                 'Example: arts_create_synthesised_beam --input_folder /data/fits '
                                                 '--sb 36'
                                     )

    parser.add_argument('--input_folder', required=True,
                        help='Folder with input FITS data')
    parser.add_argument('--sb', type=int, required=True,
                        help='Index of SB to create')
    parser.add_argument('--taskid',
                        help='Task ID of observation (required if files from multiple observations are '
                             'present in the input folder')
    parser.add_argument('--cb', type=int,
                        help='CB index of input data (required if multiple CBs of the same observation are '
                             'present in the input folder')
    parser.add_argument('--source',
                        help='Source name to set in FITS header (Default: set SB index as source name)')
    parser.add_argument('--output',
                        help='Path to output FITS file (Default: ARTS_<source name>.fits '
                             'if source name is specified, else ARTS_<taskid>_CB<cb index>_SB<sb index>.fits)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file if it already exists')
    parser.add_argument('--chunksize', type=int, default=1000,
                        help='Maximum number of subints to load at a time (Default: %(default)s')
    parser.add_argument('--no_progress_bar', action='store_true',
                        help='Disable progress bars')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose logging')

    # parse the arguments
    args = parser.parse_args()

    main(args)
