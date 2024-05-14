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

import numpy as np
import pandas as pd

sys.path.extend(["py/", "py/fetch/", "py/process/", "py/plots/"])


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
        event,
        base,
        dates,
        supermag=[],
        intermagnet=[],
        elec_stations=[],
        uid="shibaji7",
    ):
        """
        Populate all data tables from GOES, FISM2, SuperMAG and Efield.
        """
        from get_mag_data import SuperMAG
        from get_flare_data import FISM
        from goes import FlareTS

        self.event = event
        self.base = base
        self.dates = dates
        self.uid = uid
        self.elec_stations = elec_stations
        self.supermag = supermag
        self.intermagnet = intermagnet

        if not os.path.exists(base):
            os.makedirs(base)
        du = dates[0]
        du = du.replace(hour=0, minute=0)

        if len(self.elec_stations):
            self.fetch_USGS_database()
        
        self.smag = SuperMAG(
            self.base,
            self.dates,
            uid=self.uid,
            stations=self.supermag
        )
        #################################################
        # TODO: Did not implemented the intermagnet
        #################################################

        # self.fism = FISM.FetchFISM(
        #     self.base,
        #     self.dates
        # )
        self.g = FlareTS(dates)
        self.summary_plots()
        return

    def fetch_USGS_database(self):
        """
        Load USGS dataset
        """
        from get_efields import (
            download_save_EarthScope_GeoEB,
            load_MT_site,
            load_USArray_one_station_multi_cha,
        )
        from utils import compute_Egeo_from_Bgeo, detrend_magnetic_field

        self.usgs = {}
        self.usgs_sites = {}

        download_save_EarthScope_GeoEB(
            self.elec_stations, self.base + "geomag_electric/USArray/"
        )
        for sta in self.elec_stations:
            self.usgs_sites[sta] = load_MT_site(sta)
            o = load_USArray_one_station_multi_cha(
                sta,
                self.dates[0],
                self.dates[1],
                data_dir=self.base + "geomag_electric/USArray/",
            )

            if len(o) > 0:
                magX, magY = (
                    detrend_magnetic_field(np.array(o.BN)),
                    detrend_magnetic_field(np.array(o.BE)),
                )
                o["EN_sim"], o["EE_sim"], _ = compute_Egeo_from_Bgeo(
                    magX, magY, station=sta
                )
            self.usgs[sta] = o
        return

    def summary_plots(self):
        """
        Create summary plots for analysis
        """
        import os
        os.makedirs(f"{self.base}/figures/mag/", exist_ok=True)
        os.makedirs(f"{self.base}/figures/elec/", exist_ok=True)
        os.makedirs(f"{self.base}/figures/map/", exist_ok=True)
        self.g.plot_TS()
        self.g.save("goes.png", f"{self.base}/figures/")
        self.g.close()
        # self.fism.plot_TS()
        # self.fism.save("fism.png", f"{self.base}/figures/")
        # self.fism.close()

        from get_efields import plot_TS_USGS_dataset
        for i, sta in enumerate(self.elec_stations):
            o = self.usgs[sta]
            if len(o) > 0:
                plot_TS_USGS_dataset(
                    o,
                    None,
                    self.dates,
                    sta=sta,
                    figname=f"{self.base}/figures/elec/{sta}.png"
                )
        
        for stn in self.smag.stations:
            o = self.smag.sm_data[
                (self.smag.sm_data.tval == self.event)
                & (self.smag.sm_data.iaga == stn)
                & (self.smag.sm_data.sza <= 90)
            ]
            if len(o) > 0:
                self.smag.plot_TS_dataset(
                    stn,
                    figname=f"{self.base}/figures/mag/{stn}.png"
                )

        from maps import Map
        dates = [
            self.dates[0] + dt.timedelta(minutes=i) 
            for i in range(int((self.dates[1]-self.dates[0]).total_seconds()/60))
        ]
        for d in dates:
            map = Map(d)
            for i, sta in enumerate(self.smag.stations):
                o = self.smag.sm_data[
                    (self.smag.sm_data.tval == d)
                    & (self.smag.sm_data.iaga == sta)
                ]
                if len(o) > 0:
                    inst = dict(
                        lon=np.mod(o.glon.tolist()[0]+180,360)-180,
                        lat=o.glat.tolist()[0],
                        code=o.iaga.tolist()[0],
                        E_geo=o.E_geo.tolist()[0],
                        N_geo=o.N_geo.tolist()[0],
                        Z_geo=o.Z_geo.tolist()[0],
                    )
                    map.plot_instrument(inst)
            map.save(f"{self.base}/figures/map/map_{d.strftime('%H%M')}.png")
        return


def fork_event_based_mpi(file="config/eventlist.csv", ix=0):
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
        if not os.path.exists(base):
            dates = [row["s_time"], row["e_time"]]
            elec_stations = [] if "N" == row["efield_stn"] else row["efield_stn"].split("-")
            supermag = [] if "N" == row["supermag"] else row["supermag"].split("-")
            intermagnet = [] if "N" == row["intermagnet"] else row["intermagnet"].split("-")
            Hopper(ev, base, dates, supermag, intermagnet, elec_stations)
    return


if __name__ == "__main__":
    fork_event_based_mpi()
