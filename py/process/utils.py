#!/usr/bin/env python

"""utils.py: module is dedicated to utlity for processing."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import sys

import numpy as np
import pandas as pd

sys.path.extend(["py/", "py/fetch/", "py/process/"])


def detrend_magnetic_field(B, L=120, p=None):
    """
    This method is resposible for detrend
    magnetic field data and taper it to reduce
    spurious frequency components.
    """
    p = p if p is None else p
    fmed = np.median(B[:L])
    Bd = B - fmed
    if p:
        Bd *= get_tapering_function(p)
    return Bd


def calculate_GC_distance(lat1, lon1, lat2, lon2, method="GC", R=6371.0):
    """
    Compute distance between two points over a circle
    """
    from math import acos, asin, cos, radians, sin, sqrt

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    if method == "GC":
        gc = R * (
            acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
        )
    if method == "haversine":
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        gc = 2 * R * asin(sqrt(a))
    return gc


def get_all_stations(_dir_=".data/EMTF/"):
    """
    Fetch all station names from the
    EMTF directory.
    """
    stations = pd.read_csv("config/stations.csv", sep="|")
    return stations


def get_nearest_station(glat, glon):
    """
    Get nearest goegraphic station given by a
    laitiude longitude position
    """
    from get_efields import load_MT_site

    stations = get_all_stations()
    stations["GC_km"] = stations.apply(
        lambda r: calculate_GC_distance(glat, glon, r["Latitude"], r["Longitude"]),
        axis=1,
    )
    idx = stations.GC_km.argmin()
    station = stations.loc[idx]
    site = load_MT_site(station.Station)
    return station, site


def compute_Egeo_from_Bgeo(magX, magY, station=None, glat=None, glon=None):
    """
    Compute Egeo from Bgeo
    given station name or nearest
    geographic position
    """
    from get_efields import load_MT_site

    if station:
        site = load_MT_site(station)
    elif (glat is not None) and (glon is not None):
        station, site = get_nearest_station(glat, glon)
    Ex, Ey = site.convolve_fft(magX, magY, dt=1)
    return Ex, Ey, site
