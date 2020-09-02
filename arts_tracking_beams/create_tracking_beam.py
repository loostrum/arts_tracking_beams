#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import logging

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.coordinates.name_resolve import NameResolveError
from astropy.time import Time
import astropy.constants as const
import astropy.units as u

from arts_tracking_beams import ARTSFITSReader, TrackingBeam
from arts_tracking_beams.tools import radec_to_hadec
from arts_tracking_beams.constants import NTAB, CB_HPBW, REF_FREQ, BMAX


def get_input_path(input_folder, taskid=None, cb=None):
    """
    Get input file path

    :param str input_folder: Folder containing FITS files
    :param str taskid: observation taskid (required if more than one observation is present)
    :param int cb: observation compound beam (required if more than one CB is present)
    :return: file path formatted as <path>/ARTS<taskid>_CB<cb>_{tab:02d}.fits
    """
    # construct glob pattern
    if taskid is not None:
        if cb is not None:
            pattern = f'{input_folder}/ARTS{taskid}_CB{cb:02d}*.fits'
        else:
            pattern = f'{input_folder}/ARTS{taskid}*.fits'
    else:
        pattern = f'{input_folder}/*.fits'

    # get list of files
    files = glob.glob(pattern)
    if not files:
        logging.error("No input files found")
        sys.exit(1)

    # there should be one file per TAB
    if not len(files) == NTAB:
        logging.error(f'Expected {NTAB} files but found only {len(files)}')
        sys.exit(1)

    # construct the file path with formatter for TAB index
    # first sort so TAB00 file is first
    pth = sorted(files)[0].replace('TAB00', 'TAB{tab:02d}')
    logging.debug(f'Input FITS path: {pth}')

    return pth


def get_source_coordinates(name=None, ra=None, dec=None):
    """
    Get source coordinates

    :param name: source name (to resolve with SkyCoord.from_name)
    :param ra: Right Ascension (hh:mm:ss.s, ignored if name is given)
    :param dec: Declination (dd:mm:ss.s, ignored if name is given)
    :return: source coordinates (SkyCoord)
    """
    if name is not None:
        try:
            coord = SkyCoord.from_name(name)
        except NameResolveError:
            logging.error(f'Could not find coordinates for source {name}')
            sys.exit(1)
    else:
        if ra is None or dec is None:
            logging.error('RA and Dec must be specified if source name is not specified')
            sys.exit(1)
        coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    logging.debug(f'Source coordinates: {coord.to_string("hmsdms", precision=1)}')

    return coord


def required_tb_resolution(pointing, tstart, duration, fhi):
    """
    Calculate the required time interval for calculating the TB TAB indices.
    This is set by the time required to cross a single TAB at 15"/second Earth rotation,
    at the point in the observation where the TABs are narrowest (i.e. HA closest to zero)
    Worst case is roughly 1.1 second, which nicely fits the typical FITS subint length of 1.024s

    :param SkyCoord pointing: CB pointing
    :param Time tstart: Observation start time
    :param Quantity duration: Observation duration
    :param Quantity fhi: Highest frequency of observation
    :return: TB time resolution (Quantity)
    """
    # get HA, Dec at start and end of observation
    ha_s, dec_s = radec_to_hadec(pointing.ra, pointing.dec, tstart)
    ha_e, dec_e = radec_to_hadec(pointing.ra, pointing.dec, tstart + duration)
    # if one HA is negative and the other positive, the minimum HA is zero
    if (ha_s * ha_e).value < 0:
        ha = 0 * u.deg
    else:
        # find value closest to zero, sign does not matter as projection effect is symmetric around HA=0
        ha = min(np.abs(ha_s), np.abs(ha_e))
    # take the average Dec (these should be the same anyway)
    dec = .5 * (dec_s + dec_e)
    # get TAB projection angle (estimate, assumes only By is non-zero which is only approximately true)
    cos_tab_proj = np.sqrt(np.cos(ha)**2 + np.sin(dec)**2 * np.sin(ha)**2)

    # get the TAB half-power width (approximation), the 0.78 is a numerically determined scaling factor valid
    # for the 8-dish Apertif setup
    tab_hpbw = np.arcsin(.8 * const.c / (fhi * BMAX * cos_tab_proj))

    # calculate crossing time
    tcross = tab_hpbw / (360 * u.deg / u.sday)
    # require 2 steps per TAB
    tb_res = .5 * tcross

    logging.debug(f'Required TB time resolution: {tb_res.to(u.s):.1f}')
    return tb_res.to(u.s)


def main(args):
    """
    The main tracking beam program.

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
    logging.basicConfig(format='%(asctime)s.%(levelname)s.%(name)s: %(message)s', level=level)

    # identify input files
    input_fits_path = get_input_path(args.input_folder, taskid=args.taskid, cb=args.cb)

    # check if output file already exists
    if os.path.isfile(args.output) and not args.overwrite:
        logging.error("Output file already exists. Use --overwrite to overwrite.")
        sys.exit(1)

    # get source coordinates
    source_coord = get_source_coordinates(args.source, ra=args.ra, dec=args.dec)

    # initialize the FITS reader
    fits_reader = ARTSFITSReader(input_fits_path)

    # get some observation parameters
    # pointing of the CB
    pointing = SkyCoord(fits_reader.info['cb_ra'], fits_reader.info['cb_dec'])
    # observation start time and duration
    tstart = fits_reader.info['tstart']
    duration = fits_reader.info['duration']
    # frequency range
    bw = np.abs(fits_reader.info['df'] * fits_reader.info['nchan'])
    freq = fits_reader.info['freq']
    # fmin is the lowest frequency of the lowest channel, _not_ the centre frequency of the lowest channel
    fmin = freq - .5 * bw

    # check if the source is in the compound beam, just give a warning if not
    cb_hpbw = CB_HPBW * REF_FREQ / freq
    sep = pointing.separation(source_coord)
    if sep > cb_hpbw:
        logging.warning(f'Source separation from CB centre ({sep.to(u.arcmin):.1f}) is more than '
                        f'CB half-power width at {freq.to(u.MHz)} ({cb_hpbw.to(u.arcmin):.1f})')

    # determine the required TB time resolution
    tb_res = required_tb_resolution(pointing, tstart, duration, fmin + bw)

    # round tb_res to closest factor of subintegration length
    tsub = fits_reader.info['nsamp_per_subint'] * fits_reader.info['tsamp']
    tb_res = np.round(tb_res / tsub) * tsub
    logging.debug(f'Setting TB time resolution to {tb_res.to(u.s)}')


def main_with_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_folder', required=True,
                        help='Folder with input FITS data')
    parser.add_argument('--taskid',
                        help='Task ID of observation (required if files from multiple observations are '
                             'present in the input folder')
    parser.add_argument('--cb', type=int,
                        help='CB index of input data (required if multiple CBs of the same observation are '
                             'present in the input folder')
    parser.add_argument('--output', required=True,
                        help='Path to output FITS file (Default: ARTS_<RA>_<Dec>.fits,'
                             'where <RA> and <Dec> are the tracking beam coordinates in hh:mm:ss and '
                             'dd:mm:ss format, respectively)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite output file if it already exists')
    parser.add_argument('--source',
                        help='Source to create tracking beam for. Name should be resolvable by '
                             'astropy.coordinates.SkyCoord.from_name')
    parser.add_argument('--ra',
                        help='Right ascension of tracking beam in hh:mm:ss.s format '
                             '(ignored if source argument is specified)')
    parser.add_argument('--dec',
                        help='Declination of tracking beam in dd:mm:ss.s format '
                             '(ignored if source argument is specified)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose logging')

    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    main_with_args()
