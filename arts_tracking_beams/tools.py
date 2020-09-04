#!/usr/bin/env python3

from astropy.coordinates import SkyCoord, SphericalRepresentation
import astropy.units as u
from astropy.time import Time

from arts_tracking_beams.constants import WSRT_LOC


def radec_to_hadec(ra, dec, t):
    """
    Convert RA, Dec to apparent WSRT HA, Dec

    :param Quantity ra: Right ascension
    :param Quantity dec: Declination
    :param Time/str t: Observing time
    :return: HA (Quantity), Dec (Quantity)
    """

    # Convert time to Time object if given as string
    if isinstance(t, str):
        t = Time(t)

    coord = SkyCoord(ra, dec, frame='icrs', obstime=t)
    ha = WSRT_LOC.lon - coord.itrs.spherical.lon
    ha.wrap_at(12 * u.hourangle, inplace=True)
    dec = coord.itrs.spherical.lat

    return ha, dec


def hadec_to_radec(ha, dec, t):
    """
    Convert apparent HA, Dec to J2000 RA, Dec

    :param ha: hour angle with unit
    :param dec: declination with unit
    :param Time/str t: Observing time
    :return: SkyCoord object of J2000 coordinates
    """

    # Convert time to Time object if given as string
    if isinstance(t, str):
        t = Time(t)

    # create spherical representation of ITRS coordinates of given ha, dec
    itrs_spherical = SphericalRepresentation(WSRT_LOC.lon - ha, dec, 1.)
    # create ITRS object, which requires cartesian input
    coord = SkyCoord(itrs_spherical.to_cartesian(), frame='itrs', obstime=t)
    # convert to J2000
    return coord.icrs.ra, coord.icrs.dec
