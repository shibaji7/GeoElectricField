#!/usr/bin/env python

"""get_mag_data.py: module is dedicated to get magnetic data from CSV files."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import os

import matplotlib.dates as mdates
import pandas as pd
from loguru import logger
from supermag import *
from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np


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


class SuperMAG(object):
    """
    This class is dedicated to load magnetic dataset from
    SuperMAG, analyze, and plot them. The functionalities
    includes -
    1. Load magnetic station info and data for a given date
    2. Plot stations on cartopy
    3. Plot magnetic data with additional information
    """

    def __init__(
        self,
        base,
        dates,
        uid="shibaji7",
        stations=[],
        flags="mlt,mag,geo,decl,sza",
        baseline=None,
    ):
        """
        Parameters:
        -----------
        base: Base folder to store data
        dates: Start and end dates
        uid: User ID
        stations: List of station codes
        flags: Flags of data tech
        baseline: SuperMAG baseline mode. Use ``"none"`` for raw magnetic
            field values, ``"yearly"`` for explicitly yearly-baselined
            values, or ``None`` to preserve the historical request/cache.
        """
        self.base = base
        os.makedirs(base, exist_ok=True)
        self.dates = dates
        self.logon = uid
        self.flags = flags
        if baseline not in (None, "none", "yearly"):
            raise ValueError("baseline must be None, 'none', or 'yearly'")
        self.baseline = baseline
        self.stations = self.load_stations() if len(stations)==0 else stations
        self.fetch_baseline_substracted_datasets()
        return

    def load_stations(self):
        """
        Load stations from SuperMAG
        """
        extent = int((self.dates[1] - self.dates[0]).total_seconds())
        (status, stations) = SuperMAGGetInventory(self.logon, self.dates[0], extent)
        logger.info(f"SM inventory fetch stats: {status}, {len(stations)}")
        return stations

    def fetch_baseline_substracted_datasets(self):
        """
        Fetch baseline substracted datasets.
        """
        self.sm_data = pd.DataFrame()
        extent = int((self.dates[1] - self.dates[0]).total_seconds())
        if self.baseline == "none":
            fname = self.base + "supermag-raw.csv"
        elif self.baseline == "yearly":
            fname = self.base + "supermag-baseline-yearly.csv"
        else:
            fname = self.base + "supermag.csv"
        request_flags = self.flags
        if self.baseline is not None:
            request_flags = ",".join(
                flag for flag in request_flags.split(",")
                if not flag.strip().lower().startswith("baseline=")
                and flag.strip().lower() != "delta=start"
                and flag.strip()
            )
            request_flags += ",baseline=" + self.baseline
        if not os.path.exists(fname):
            logger.info(f"Total stations to be fetched {len(self.stations)}")
            for i, stn in tqdm(enumerate(self.stations)):
                (status, data) = SuperMAGGetData(
                    self.logon, self.dates[0], extent, request_flags, stn
                )
                if len(data) > 0:
                    required = {"N", "E", "Z", "tval"}
                    missing = required.difference(data.columns)
                    if missing:
                        logger.warning(
                            f"Skipping station {stn}: missing fields {sorted(missing)}"
                        )
                        logger.info(
                            f"SM inventory fetch data stats[{stn}]: {status}, idx:{i}"
                        )
                        continue
                    if len(data) < 2:
                        logger.warning(
                            f"Skipping station {stn}: fewer than two data records"
                        )
                        logger.info(
                            f"SM inventory fetch data stats[{stn}]: {status}, idx:{i}"
                        )
                        continue
                    # if np.mean(data.sza) < 110.0:
                    for comp in ["N", "E", "Z"]:
                        for cord in ["nez", "geo"]:
                            data[comp + "_" + cord] = sm_grabme(data, comp, cord)
                    data.drop(["N", "E", "Z"], axis=1, inplace=True)
                    data.tval = data.tval.apply(
                        lambda x: dt.datetime.utcfromtimestamp(x)
                    )
                    data["B"] = np.sqrt(
                        data["N_geo"]**2 + data["E_geo"]**2 + data["Z_geo"]**2
                    )
                    delt = (data.tval.iloc[1] - data.tval.iloc[0]).total_seconds()
                    data["dBdt"] = np.diff(data.B, prepend=data.B.iloc[0]) / delt
                    self.sm_data = pd.concat([self.sm_data, data])
                logger.info(f"SM inventory fetch data stats[{stn}]: {status}, idx:{i}")
            self.sm_data.to_csv(fname, header=True, index=False, float_format="%g")
        else:
            logger.info(f"SM inventory is local: {0}")
            if os.path.getsize(fname) == 0:
                logger.warning(f"SuperMAG cache is empty: {fname}")
                return
            try:
                self.sm_data = pd.read_csv(fname, parse_dates=["tval"])
            except pd.errors.EmptyDataError:
                logger.warning(f"SuperMAG cache is empty: {fname}")
                self.sm_data = pd.DataFrame()
                return
            self.sm_data["glon"] = np.round(np.mod(self.sm_data["glon"]+180,360)-180, 2)
        return

    def plot_TS_dataset(
        self,
        stn,
        ax=None,
        ylim=[-100, 100],
        coords="geo",
        comps={
            "N": {"color": "r", "ls": "-", "lw": 0.5},
            "E": {"color": "g", "ls": "-", "lw": 0.5},
            "Z": {"color": "b", "ls": "-", "lw": 0.5},
        },
        xlabel="UT",
        ylabel=r"$\delta$, nT",
        loc=2,
        figname=None,
    ):
        """
        Overlay station data into axes
        """
        if ax == None:
            fig = plt.figure(dpi=240, figsize=(5, 3))
            ax = fig.add_subplot(111)
        data = self.sm_data[
            self.sm_data.iaga == stn
        ]
        if ylim:
            ax.set_ylim(ylim)
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        ax.set_xlim(self.dates)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        for comp in comps.keys():
            co = comps[comp]
            ax.plot(
                data.tval,
                data[comp + "_" + coords],
                ls=co["ls"],
                color=co["color"],
                lw=co["lw"],
                label=comp,
            )
        ax.legend(loc=loc, fontsize=6)
        txt = f"{stn}" + "\n" + rf"$\theta,\phi$={data.glat.tolist()[0]}, {data.glon.tolist()[0]}" + \
           "\n"+ fr"$\chi$={'%.2f'%np.mean(data.sza)}"
        ax.text(0.9, 0.9, txt, ha="center", va="center", transform=ax.transAxes)

        if figname:
            fig.savefig(figname, bbox_inches="tight")
        return ax

    @staticmethod
    def FetchSM(base, dates, uid="shibaji7", stations=None):
        """
        Static method to call the supermag
        inevntory functions, pre-process them
        and store to object.
        """
        sm = SuperMAG(base, dates, uid, stations)
        return sm


if __name__ == "__main__":
    SuperMAG.FetchSM(
        "data/2017-09-06-11-56/",
        [dt.datetime(2017, 9, 6, 11, 50), dt.datetime(2017, 9, 6, 12, 20)],
    )
