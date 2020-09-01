#!/usr/bin/env python3

from astropy.coordinates import SkyCoord
import astropy.units as u

from arts_tracking_beams.constants import WSRT_LOC


def radec_to_hadec(ra, dec, t):
    """
    Convert RA, Dec to apparent WSRT HA, Dec

    :param Quantity ra: Right ascension
    :param Quantity dec: Declination
    :param Time t: Observing time
    :return: HA (Quantity), Dec (Quantity)
    """
    coord = SkyCoord(ra, dec, frame='icrs', obstime=t)
    ha = WSRT_LOC.lon - coord.itrs.spherical.lon
    ha.wrap_at(12 * u.hourangle, inplace=True)
    dec = coord.itrs.spherical.lat

    return ha, dec
