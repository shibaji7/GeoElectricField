#!/usr/bin/env python

"""pre_process.py: module is dedicated to fetch, filter, and save data for post_processing."""

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
import sys

import pandas as pd

sys.path.extend(["py/", "py/fetch/", "py/process/"])


class Hopper(object):
    """
    This class is responsible for following
    operations for each flare event.
        i) Fetching SD, SM, GOES, and
        FISM dataset using .fetch module.
        ii) Populate coordinates (mlat, mlon).
        iii) Median filter.
        iv) Filter the data using SZA, tfreq, slant-range.
        v) Summary plots for the dataset.
        vi) Store the data for post-processing.
    """

    def __init__(
        self,
        base,
        dates,
        rads,
        elec_stations=[],
        sat_resolution=2,
        sat_station=15,
        fism_spectrum=[0.01, 40, 1.0],
        uid="shibaji7",
    ):
        """
        Populate all data tables from GOES, FISM2, SuperDARN, and SuperMAG.
        """
        from get_fit_data import FetchData
        from get_flare_data import GOES

        self.base = base
        self.dates = dates
        self.rads = rads
        self.sat_resolution = sat_resolution
        self.sat_station = sat_station
        self.fism_spectrum = fism_spectrum
        self.uid = uid
        self.elec_stations = elec_stations

        if not os.path.exists(base):
            os.makedirs(base)
        self.sds = []
        du = dates[0]
        du = du.replace(hour=0, minute=0)
        for r in rads:
            sd = FetchData.FetchSD(
                base,
                r,
                [du, du + dt.timedelta(1)],
            )
            self.sds.append(sd)
        self.g = GOES.FetchGOES(
            base,
            dates,
            sat_station,
            sat_resolution,
        )
        self.fetch_USGS_database()
        self.summary_plots()
        return

    def fetch_USGS_database(self):
        """
        Load USGS dataset
        """
        from get_efields import (
            download_save_EarthScope_GeoEB,
            load_USArray_one_station_multi_cha,
        )

        self.usgs = {}

        download_save_EarthScope_GeoEB(
            self.elec_stations, self.base + "geomag_electric/USArray/"
        )
        for sta in self.elec_stations:
            self.usgs[sta] = load_USArray_one_station_multi_cha(
                sta,
                self.dates[0],
                self.dates[1],
                data_dir=self.base + "geomag_electric/USArray/",
            )
        return

    def summary_plots(self):
        """
        Create summary plots for analysis
        """
        from get_efields import plot_TS_USGS_dataset
        import matplotlib.pyplot as plt

        fig_folder = self.base + "figures/"
        if not os.path.exists(fig_folder):
            os.makedirs(fig_folder)

        plt.style.use(["science", "ieee"])
        fig = plt.figure(dpi=240, figsize=(5, 3 * (1 + len(self.elec_stations))))
        
        # Plot GOES dataset
        ax = fig.add_subplot(100 * (1 + len(self.elec_stations)) + 11)
        self.g.plot_TS_dataset(ax=ax)
        
        for i, sta in enumerate(self.elec_stations):
            o = self.usgs[sta]
            if len(o) > 0:
                ax = fig.add_subplot(100 * (1 + len(self.elec_stations)) + 11 + i)
                plot_TS_USGS_dataset(
                    o,
                    ax,
                    self.dates,
                    sta=sta,
                )
            
        fig.savefig(fig_folder + "summary.png", bbox_inches="tight")
        return


def fork_event_based_mpi(file="config/events.csv"):
    """
    Load all the events from
    events list files and fork Hopper
    to pre-process and store the
    dataset.
    """
    o = pd.read_csv(file, parse_dates=["event", "s_time", "e_time"])
    for i, row in o.iterrows():
        ev = row["event"]
        base = "tmp/{Y}-{m}-{d}-{H}-{M}/".format(
            Y=ev.year,
            m="%02d" % ev.month,
            d="%02d" % ev.day,
            H="%02d" % ev.hour,
            M="%02d" % ev.minute,
        )
        dates = [row["s_time"], row["e_time"]]
        rads = row["rads"].split("-")
        elec_stations = row["efield_stn"].split("-")
        Hopper(base, dates, rads, elec_stations)
    return


if __name__ == "__main__":
    fork_event_based_mpi()
