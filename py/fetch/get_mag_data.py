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
        stations=None,
        flags="mlt,mag,geo,decl,sza",
    ):
        """
        Parameters:
        -----------
        base: Base folder to store data
        dates: Start and end dates
        uid: User ID
        stations: List of station codes
        flags: Flags of data tech
        """
        self.base = base
        os.makedirs(base, exist_ok=True)
        self.dates = dates
        self.logon = uid
        self.flags = flags
        self.stations = self.load_stations() if stations is None else stations
        self.fetch_baseline_substracted_datasets()
        return

    def load_stations(self):
        """
        Load stations from SuperMAG
        """
        extent = int((self.dates[1] - self.dates[0]).total_seconds())
        (status, stations) = SuperMAGGetInventory(self.logon, self.dates[0], extent)
        logger.info(f"SM inventory fetch stats: {status}")
        return stations

    def fetch_baseline_substracted_datasets(self):
        """
        Fetch baseline substracted datasets.
        """
        self.sm_data = pd.DataFrame()
        extent = int((self.dates[1] - self.dates[0]).total_seconds())
        fname = self.base + "supermag.csv"
        if not os.path.exists(fname):
            logger.info(f"Total stations to be fetched {len(self.stations)}")
            for i, stn in tqdm(enumerate(self.stations)):
                (status, data) = SuperMAGGetData(
                    self.logon, self.dates[0], extent, self.flags, stn
                )
                if len(data) > 0:
                    # if np.mean(data.sza) < 110.0:
                    for comp in ["N", "E", "Z"]:
                        for cord in ["nez", "geo"]:
                            data[comp + "_" + cord] = sm_grabme(data, comp, cord)
                    data.drop(["N", "E", "Z"], axis=1, inplace=True)
                    data.tval = data.tval.apply(
                        lambda x: dt.datetime.utcfromtimestamp(x)
                    )
                    self.sm_data = pd.concat([self.sm_data, data])
                logger.info(f"SM inventory fetch data stats[{stn}]: {status}, idx:{i}")
            self.sm_data.to_csv(fname, header=True, index=False, float_format="%g")
        else:
            logger.info(f"SM inventory is local: {0}")
            self.sm_data = pd.read_csv(fname, parse_dates=["tval"])
        return

    def plot_TS_dataset(
        self,
        stn,
        ax,
        ylim=[-50, 50],
        coords="geo",
        comps={
            "N": {"color": "r", "ls": "-", "lw": 0.5},
            "E": {"color": "r", "ls": "-", "lw": 0.5},
            "Z": {"color": "r", "ls": "-", "lw": 0.5},
        },
        xlabel="UT",
        ylabel=r"$\delta$, nT",
        loc=2,
    ):
        """
        Overlay station data into axes
        """
        data = self.sm_data[stn]
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
