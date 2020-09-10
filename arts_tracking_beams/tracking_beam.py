#!/usr/bin/env python3

import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.time import Time

from arts_tracking_beams.tools import radec_to_hadec
from arts_tracking_beams.constants import WSRT_LOC, ARRAY_ITRF, DISH_ITRF, NTAB, BCQ


class TrackingBeam:

    def __init__(self, ra, dec, ra_src, dec_src, fmin=1220 * u.MHz, bw=300 * u.MHz, nsub=48):
        self.ra0 = ra
        self.dec0 = dec
        self.ra0_src = ra_src
        self.dec0_src = dec_src

        # set subband frequencies
        df_sub = bw / nsub
        self.f_sub = fmin + df_sub / 2 + np.arange(nsub) * df_sub

        # set XYZ baselines
        self.baselines = self._get_baselines()

    def run(self, t):
        """
        Calculate the TAB indices for given time stamp

        :param Time/array t: Time stamp or array of time stamps
        :return: TAB indices, shape (len(t), NTAB)
        """

        # check if t is one value or a list/array
        self.time_is_array = False
        self.ntime = None
        if isinstance(t, np.ndarray) or isinstance(t, list):
            self.time_is_array = True
        elif isinstance(t, Time):
            # if created with Time, the input could be a single Time instances with multiple values, check
            # with Time.value
            if isinstance(t.value, np.ndarray):
                self.time_is_array = True
        if self.time_is_array:
            self.ntime = len(t)

        # get HA, Dec coordinates of phase center and source
        self.ha, self.dec = radec_to_hadec(self.ra0, self.dec0, t)
        self.ha_src, self.dec_src = radec_to_hadec(self.ra0_src, self.dec0_src, t)
        # get UVW coordinates
        self.uvw = self._get_lambda_uvw()
        # get TAB rotation and projection
        self.tab_rot, self.tab_proj = self._get_tab_angles()
        # get TAB index at each subband
        return self._get_tab_indices()

    @staticmethod
    def _get_baselines():
        """
        Convert ITRF dish and array positions to local XYZ baselines
        """
        # first rotate both dish and array position to local coordinates
        lon = WSRT_LOC.lon
        rot_matrix = np.array([[np.cos(-lon), -np.sin(-lon), 0],
                               [np.sin(-lon), np.cos(-lon), 0],
                               [0, 0, 1]])
        dish_xyz = (rot_matrix @ DISH_ITRF.T).T
        array_xyz = rot_matrix @ ARRAY_ITRF
        return dish_xyz - array_xyz

    def _get_lambda_uvw(self):
        """
        Get uvw coordinates from HA, Dec
        Units are meters, i.e. not scaled by wavelength
        """
        ha = self.ha
        dec = self.dec
        rot_matrix = np.array([[np.sin(ha), np.cos(ha), 0],
                               [-np.sin(dec) * np.cos(ha), np.sin(dec) * np.sin(ha), np.cos(dec)],
                               [np.cos(dec) * np.cos(ha), -np.cos(dec) * np.sin(ha), np.sin(dec)]])

        uvw = np.matmul(rot_matrix, self.baselines.T)
        if self.time_is_array:
            # uvw is a (3, ndish) array, but each value is a Quantity with ntime values
            # turn into one big Quantity object
            uvw_old = uvw.copy()
            ndish = len(self.baselines)
            uvw = np.zeros((3, ndish, self.ntime)) * u.dimensionless_unscaled
            for a in range(3):
                for dish in range(ndish):
                    uvw[a, dish] = uvw_old[a, dish].value
            uvw *= uvw_old.unit
        return uvw

    def _get_tab_angles(self):
        uu, vv, ww = self.uvw
        baseline_lengths = np.sqrt(np.sum(self.uvw ** 2, axis=0))

        # rotation: add 90 deg to get East from North
        # tab_rot = np.arctan2(vv, uu)  # this needs some more checking (and a modulo operation?)
        # version assuming only By is nonzero:
        tab_rot = np.arctan2(np.sin(self.dec) * np.sin(self.ha), np.cos(self.ha))
        # projection: only interested in absolute value
        # could also be calculated as arccos(sqrt((u**2 + v**2) / B**2))
        tab_proj = np.abs(np.arcsin(ww / baseline_lengths))

        # take median over dishes (not mean, because values for RT8 can be close to zero divided by zero)
        tab_rot = np.median(tab_rot)
        tab_proj = np.median(tab_proj)

        return tab_rot, tab_proj

    def _get_tab_indices(self):
        # get source offset from CB centre
        ddec_src = self.dec_src - self.dec
        meandec = .5 * (self.dec + self.dec_src)
        dhacosdec_src = (self.ha_src - self.ha) * np.cos(meandec)
        # de-rotate source to get x-axis position in TAB reference frame
        x_src = np.cos(self.tab_rot) * dhacosdec_src - np.sin(self.tab_rot) * ddec_src
        # find position in grating response
        # first get grating distance, use arcsin to automaticaly get units
        # though small angle approximation should hold here
        tab_gr = np.arcsin(const.c / (self.f_sub * BCQ * np.cos(self.tab_proj)))
        # get pointing shift per tab
        pointing_shift_per_tab = tab_gr / NTAB
        # get TAB index in range [0, NTAB)
        if self.time_is_array:
            # add extra axis so output shape is ntime, nsub
            x_src = x_src[:, None]
        return np.round((x_src / pointing_shift_per_tab).to(1).value).astype(int) % NTAB
