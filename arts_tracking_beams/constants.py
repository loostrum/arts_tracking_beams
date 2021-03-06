#!/usr/bin/env python3
#
# WSRT constants

import os
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation

# path to static files
STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'static')

#: File with definition of synthesised beams
SB_TABLE = os.path.join(STATIC_DIR, 'sbtable-sc4-12tabs-71sbs.txt')
#: Number of compound beams
NCB = 40
#: Number of tied-array beams
NTAB = 12
#: Number of synthesised beams
NSB = 71
#: CB half-power width
CB_HPBW = 28.0835088 * u.arcmin
#: Reference frequency for CB half-power width
REF_FREQ = 1770 * u.MHz
#: Common-quotient baseline between the dishes
BCQ = 144 * u.m
#: Maximum baseline in Apertif-8 setup
BMAX = 1008 * u.m

#: ITRF WSRT reference position
ARRAY_ITRF = np.array([3828630.63486200943589211, 443593.39226634375518188, 5064922.99755000043660402]) * u.m
WSRT_LOC = EarthLocation.from_geocentric(*ARRAY_ITRF)

#: ITRF positions of RT2 - RT9
DISH_ITRF = np.array([[3828729.99081358872354031, 442735.17696416645776480, 5064923.00829000025987625],
                      [3828713.43109884625300765, 442878.21189340209821239, 5064923.00435999967157841],
                      [3828696.86994427768513560, 443021.24917263782117516, 5064923.00396999996155500],
                      [3828680.31391932582482696, 443164.28596862131962553, 5064923.00035000033676624],
                      [3828663.75159173039719462, 443307.32138055720133707, 5064923.00203999970108271],
                      [3828647.19342757249251008, 443450.35604637680808082, 5064923.00229999981820583],
                      [3828630.63486200943589211, 443593.39226634375518188, 5064922.99755000043660402],
                      [3828614.07606798363849521, 443736.42941620573401451, 5064923.00000000000000000]]) * u.m
