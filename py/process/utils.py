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


def get_tapering_function(B, p, dT=1):
    """
    This method is resposible for generateing
    tapering function based on time sequence t
    and tapering coefficient p
    """
    t = np.arange(len(B)) * dT
    T = len(B)
    P, P2 = int(T * p), int(T * p / 2)
    w = np.zeros_like(B)
    w[:P2] = 0.5 * (1 - np.cos(2 * np.pi * t[:P2] / P))
    w[P2 : T - P2] = 1.0
    w[T - P2 :] = 0.5 * (1 - np.cos(2 * np.pi * (t[-1] - t[T - P2 :]) / P))
    return w


def detrend_magnetic_field(B, L=120, p=None, dT=1):
    """
    This method is resposible for detrend
    magnetic field data and taper it to reduce
    spurious frequency components.
    """
    p = p if p is None else p
    fmed = np.median(B[:L])
    Bd = B - fmed
    if p:
        Bd *= get_tapering_function(B, p, dT)
    return Bd


def calculate_GC_distance(lat1, lon1, lat2, lon2, method="GC", R=6371.0):
    """
    Compute distance between two points over a circle
    """
    from math import acos, asin, cos, radians, sin, sqrt

    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    if method == "GC":
        arg = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)
        arg = max(-1.0, min(1.0, arg))
        gc = R * acos(arg)
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


def get_nearest_station(glat, glon, MT_data_folder=".data/EMTF/"):
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
    stations = stations.sort_values("GC_km").reset_index(drop=True)

    for _, station in stations.iterrows():
        site = load_MT_site(station.Station, MT_data_folder=MT_data_folder)
        if site is not None:
            return station, site

    return None, None


def compute_Egeo_from_Bgeo(magX, magY, station=None, glat=None, glon=None, dT=1, MT_data_folder=".data/EMTF/"):
    """
    Compute Egeo from Bgeo
    given station name or nearest
    geographic position
    """
    from get_efields import load_MT_site
    if station:
        site = load_MT_site(station, MT_data_folder=MT_data_folder)
        if site is None and (glat is not None) and (glon is not None):
            station, site = get_nearest_station(glat, glon, MT_data_folder=MT_data_folder)
    elif (glat is not None) and (glon is not None):
        station, site = get_nearest_station(glat, glon, MT_data_folder=MT_data_folder)
    else:
        station = None
    if site:
        Ex, Ey = site.convolve_fft(magX, magY, dt=dT)
    else:
        Ex, Ey = [], []
    return Ex, Ey, site, station
