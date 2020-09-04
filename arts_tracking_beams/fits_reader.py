#!/usr/bin/env python3

import os
import sys
import logging
import warnings

import numpy as np
from astropy.time import Time
from astropy.io import fits
from astropy.io.fits.fitsrec import FITS_rec
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io.fits.verify import VerifyWarning

from arts_tracking_beams.constants import NTAB


class ARTSFITSReader:

    def __init__(self, fname, ntab=NTAB):
        """
        FITS reader for ARTS data, one file per TAB

        :param str fname: path to FITS files, with {tab:02d} in place of the tab index
        :param int ntab: Number of TABs (Default: NTAB from constants)
        """
        self.logger = logging.getLogger(__class__.__name__)
        self.ntab = ntab

        self.fnames = [fname.format(tab=tab) for tab in range(ntab)]
        self._verify_files()

        self.file_handles = []
        self.files_open = False

        self.info, self.hdr_p, self.hdr_s = self._get_fits_params()

    def _verify_files(self):
        """
        Check if the input TAB files exist, warn if files are missing
        """
        missing_files = []
        for f in self.fnames:
            if not os.path.isfile(f):
                missing_files.append(f)
        # if user did not add {tab} to file name, the same name occurs NTAB times, so remove duplicates
        missing_files = list(set(missing_files))
        if missing_files:
            self.logger.warning(f"Cannot find TAB file(s): {', '.join(missing_files)}")
        else:
            self.logger.debug("Found all TAB files")

    def _get_fits_params(self):
        """
        Reader FITS parameters from the TAB00 file

        :return: FITS info (dict)
        """
        self.logger.debug("Reading parameters from TAB00 file")
        info = {}
        with fits.open(self.fnames[0]) as f:
            hdr_p = f['PRIMARY'].header
            hdr_s = f['SUBINT'].header

        # data from primary header
        # use SkyCoord to convert hh:mm:ss, dd:mm:ss to decimal degrees
        radec = SkyCoord(hdr_p['RA'], hdr_p['DEC'], unit=(u.hourangle, u.deg))
        info['cb_ra'] = radec.ra
        info['cb_dec'] = radec.dec
        info['freq'] = hdr_p['OBSFREQ'] * u.MHz
        # get start MJD from the different keys in the header
        info['tstart'] = Time(hdr_p['STT_IMJD'] + (hdr_p['STT_SMJD'] + hdr_p['STT_OFFS']) / 86400., format='mjd')

        # data from subint header
        info['nchan'] = hdr_s['NCHAN']
        info['tsamp'] = hdr_s['TBIN'] * u.s
        info['df'] = hdr_s['CHAN_BW'] * u.MHz
        info['nsamp_per_subint'] = hdr_s['NSBLK']
        info['table_width'] = hdr_s['NAXIS1']
        info['nsubint'] = hdr_s['NAXIS2']

        # set duration from sample parameters, SCANLEN seems to be truncated
        info['duration'] = info['nsubint'] * info['nsamp_per_subint'] * info['tsamp']
        return info, hdr_p, hdr_s

    def open_files(self):
        """
        Open a handle to each FITS file
        """
        if self.files_open:
            self.logger.warning("Files already opened")
            return
        self.logger.debug("Opening files")
        self.file_handles = [fits.open(f, memmap=True) for f in self.fnames]
        self.files_open = True

    def close_files(self):
        """
        Close the handle to each FITS file
        """
        if not self.files_open:
            self.logger.warning("Files already closed")
            return
        self.logger.debug("Closing files")
        for h in self.file_handles:
            h.close()
        self.files_open = False

    def read_tabs_to_buffer(self, subintstart, nsubint, tabs, out):
        """
        Read a chunk of data into output FITS record. A different TAB can be loaded per frequency subbands
        by setting tabs to an array of TAB indices per subband, lowest frequency first.

        :param int subintstart: First subintegration to load
        :param int nsubint: Number of subintegrations to loadd
        :param int/list tabs: tab to load per subband, from low to high freq (same tab for all freqs if int)
        :param FITS_rec out: output array to store data to
        """
        # if one tab, turn into list
        if isinstance(tabs, int):
            tabs = [tabs]

        # reverse tabs list if frequencies are ordered high to low (which they usually are)
        if self.info['df'] < 0:
            tabs = list(reversed(tabs))

        # get number of subbands
        nsub = len(tabs)
        nchan = self.info['nchan']
        if not nchan % nsub == 0:
            self.logger.error(f"Number of channels ({nchan}) is not divisible by requested number of subbands ({nsub})")
            sys.exit(1)
        # get number of bytes per subband
        if not (nchan // nsub) % 8 == 0:
            self.logger.error(f"Number of channels per subband must be a multiple of 8")
            sys.exit(1)
        nchan_per_subband = (nchan // nsub)
        nbyte_per_subband = nchan_per_subband // 8

        nsamp = self.info['nsamp_per_subint']
        ind = np.zeros((nsamp * nbyte_per_subband), dtype=int)
        # load data of each subband
        # ignore the VerifyWarning that happens because 8 samples are packed into one byte
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', category=VerifyWarning)
            for subband, tab in enumerate(tabs):
                # get the file handle for this tab
                handle = self.file_handles[tab]['SUBINT']
                # get indices of frequency chunks to extract
                freq_start = subband * nchan_per_subband
                freq_end = (subband + 1) * nchan_per_subband
                # the same in bytes
                freq_byte_start = subband * nbyte_per_subband
                freq_byte_end = (subband + 1) * nbyte_per_subband
                # set weights, scales, offsets for this part of the frequency band
                for col in ['DAT_WTS', 'DAT_SCL', 'DAT_OFFS']:
                    # here use normal frequency indices
                    out.data[col][subintstart:subintstart + nsubint, freq_start:freq_end] = \
                        handle.data[col][subintstart:subintstart + nsubint, freq_start:freq_end]
                # we need to extract these indices for each subint, each time sample in that subint
                # construct array of all samples that need to be extract from a subint
                for s in range(nsamp):
                    # nchan // 8 here indicates the total number of bytes of one time step
                    # use byte indices for frequency as this these are packed data
                    samp_start = s * (nchan // 8) + freq_byte_start
                    samp_end = s * (nchan // 8) + freq_byte_end
                    ind[s * nbyte_per_subband:(s + 1) * nbyte_per_subband] = np.arange(samp_start, samp_end)
                # now extract the data
                out.data['DATA'][subintstart:subintstart + nsubint, ind] = \
                    handle.data['DATA'][subintstart:subintstart + nsubint, ind]
