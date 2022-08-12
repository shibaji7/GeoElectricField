#!/usr/bin/env python

"""simulation.py: module is dedicated to Ex,Ey estimates from simulation study."""

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

sys.path.extend(["py/", "py/fetch/", "py/process/"])


class Simulate(object):
    """
    This class is dedicated to simulate Ex, Ey from magnetic data.
    """

    def __init__(
        self,
        event,
        dstart=20,
        dend=30,
        mt_stations=[],
        sm_stations=[],
        sat_resolution=3,
        sat_station=10,
        base="tmp/{Y}-{m}-{d}-{H}-{M}/",
    ):
        """
        Initialize all parameters
        """
        from get_flare_data import GOES

        self.event = event
        self.dstart = self.event - dt.timedelta(minutes=dstart)
        self.dend = self.event + dt.timedelta(minutes=dend)
        self.sm_stations = sm_stations
        self.mt_station = mt_stations
        self.base = base.format(
            Y=event.year,
            m="%02d" % event.month,
            d="%02d" % event.day,
            H="%02d" % event.hour,
            M="%02d" % event.minute,
        )
        if not os.path.exists(self.base):
            os.makedirs(self.base)
        self.g = GOES.FetchGOES(
            self.base,
            [self.dstart, self.dend],
            sat_station,
            sat_resolution,
        )
        if len(self.mt_station) > 0:
            self.fetch_USGS_database()
        if len(self.sm_stations) > 0:
            self.fetch_SM_database()
            
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

        self.usgs = {}
        self.usgs_sites = {}

        download_save_EarthScope_GeoEB(
            self.mt_stations, self.base + "geomag_electric/USArray/"
        )
        for sta in self.mt_stations:
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
                o["EN_sim"], o["EE_sim"] = compute_Egeo_from_Bgeo(
                    magX, magY, station=sta
                )
            self.usgs[sta] = o
        return

    def fetch_SM_database(self):
        """
        Fetch SuperMAG data and
        compute Ex, Ey from simulation
        """
        from get_mag_data import SuperMAG
        from utils import compute_Egeo_from_Bgeo, detrend_magnetic_field

        self.sm_site = SuperMAG(
            self.base,
            [self.dstart, self.dend],
            stations=self.sm_stations,
        )
        self.sm_efield_site = {}
        for sta in self.sm_stations:
            o = self.sm_site.sm_data[self.sm_site.sm_data.iaga == sta]
            if len(o) > 0:
                df = pd.DataFrame()
                magX, magY = (
                    detrend_magnetic_field(np.array(o.N_geo)),
                    detrend_magnetic_field(np.array(o.E_geo)),
                )
                glat, glon = o.glat.tolist()[0], o.glon.tolist()[0]
                o["EN_sim"], o["EE_sim"], site = compute_Egeo_from_Bgeo(
                    magX, magY, glat=glat, glon=glon
                )
                self.sm_efield_site[sta] = o
        return

    def summary_plots(self):
        """
        Create summary plots for analysis
        """
        import matplotlib.pyplot as plt

        fig_folder = self.base + "figures/"
        if not os.path.exists(fig_folder):
            os.makedirs(fig_folder)

        L = len(self.sm_stations) + 1 if len(self.sm_stations) < 8 else 9
        plt.style.use(["science", "ieee"])
        fig = plt.figure(dpi=240, figsize=(5, 3 * L))

        # Plot GOES dataset
        ax = fig.add_subplot(100 * L + 11)
        self.g.plot_TS_dataset(ax=ax)

        # Plot E_geo and B_geo data

        fig.savefig(fig_folder + "summary.png", bbox_inches="tight")
        return


if __name__ == "__main__":
    sim = Simulate(
        event=dt.datetime(2005, 9, 9, 19, 30),
        dstart=30,
        dend=120,
        sm_stations=[],
        sat_resolution=3,
        sat_station=12,
    )
