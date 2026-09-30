#!/usr/bin/env python

"""get_flare_data.py: module is dedicated to get flare data and store to local directory."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import calendar
import datetime as dt
import os

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
import requests
import swifter
from loguru import logger

def setup(science=True):
    if science:
        import mplstyle
        plt.rcParams.update(
            {
                "figure.figsize": np.array([8, 6]),
                "text.usetex": True,
                "font.family": "sans-serif",
                "font.sans-serif": [
                    "Tahoma",
                    "DejaVu Sans",
                    "Lucida Grande",
                    "Verdana",
                ],
                "font.size": 12,
            }
        )
    return


class FISM(object):
    """
    This class is dedicated to load magnetic dataset from
    FISM2, analyze, and plot them. The functionalities
    includes -
    1. Load FISM2 data for a given date
    2. Plot data with additional information
    """

    def __init__(
        self,
        base,
        dates,
        spectrum=[0.05, 5, 1.0],
    ):
        """
        Parameters:
        -----------
        base: Base folder to store data
        dates: Start and end dates
        spectrum: Irradiance spectrum in nm
        """
        self.base = base
        os.makedirs(base, exist_ok=True)
        self.dates = dates
        self.wavelengths = np.arange(spectrum[0], spectrum[1], spectrum[2])
        self.spectrum = spectrum
        self.df = self.fetch_1m_datasets()
        return

    def fetch_1m_datasets(self):
        """
        Fetch 1m averaged FISM2-Flare data.
        """
        URL = """https://lasp.colorado.edu/lisird/latis/dap/fism_flare_hr.csv?\
                &time>={t0}.000Z&time<={t1}.000Z&wavelength~{w}"""
        spectrum = np.arange(self.spectrum[0], self.spectrum[1], self.spectrum[2])
        tmpfname, fism_fname = (
            self.base + "fism_flare_hr.csv",
            self.base + "FISM2.high.csv",
        )
        if not os.path.exists(fism_fname):
            o = pd.DataFrame()
            for wv in spectrum:
                uri = URL.format(
                    w=wv,
                    t0=self.dates[0].strftime("%Y-%m-%dT%H:%M:%S"),
                    t1=self.dates[1].strftime("%Y-%m-%dT%H:%M:%S"),
                ).replace(" ", "")
                logger.info(f"{uri}")
                r = requests.get(url=uri)
                if r.status_code == 200:
                    with open(tmpfname, "w") as h:
                        h.write(r.text)
                    d = pd.read_csv(tmpfname)
                    d = d.rename(
                        columns={
                            "time (Julian Date)": "time",
                            "wavelength (nm)": "wavelength",
                            "irradiance (W/m^2/nm)": "irradiance",
                            "uncertainty (unitless)": "uncertainty",
                        }
                    )
                    d.time = d.time.swifter.apply(
                        lambda x: pd.to_datetime(x, origin='julian', unit='D')
                    )
                    o = pd.concat([o, d])
            os.remove(tmpfname)
            if len(o) > 0:
                o.to_csv(fism_fname, index=False, header=True, float_format="%g")
        else:
            logger.info(f"File {fism_fname} exists!")
            o = pd.read_csv(fism_fname, parse_dates=["time"])
        return o

    @staticmethod
    def FetchFISM(
        base,
        dates,
        spectrum=[0.05, 5, 1.0],
    ):
        """
        Static method to call the FISM
        inevntory functions, pre-process them
        and store to object.
        """
        f = FISM(base, dates, spectrum)
        return f

    def plot_TS(self, dates=None):
        """
        Plot time series data in-memory
        """
        dates = dates if dates else self.dates
        setup()
        self.fig = plt.figure(figsize=(6, 2.5), dpi=150)
        ax = self.fig.add_subplot(111)
        self.plot_TS_from_axes(ax, dates)
        return

    def plot_TS_from_axes(self, ax, dates=None, xlabel="Time (UT)"):
        """
        Plot time series data in-memory
        """
        dates = dates if dates else self.dates
        ax.set_xlabel(xlabel, fontdict={"size": 15, "fontweight": "bold"})
        ax.set_ylabel(
            r"Irradiance ($W/m^2$)", fontdict={"size": 15, "fontweight": "bold"}
        )
        ax.text(
            1.05,
            0.99,
            "FISM2",
            ha="center",
            va="top",
            fontdict={"size": 12, "fontweight": "bold"},
            transform=ax.transAxes,
            rotation=90,
        )
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        hours = mdates.HourLocator(byhour=range(0, 24, 1))
        ax.xaxis.set_major_locator(hours)
        dtime = (dates[1] - dates[0]).total_seconds() / 3600.0
        if dtime < 4.0:
            minutes = mdates.MinuteLocator(byminute=range(0, 60, 10))
            ax.xaxis.set_minor_locator(minutes)
            ax.xaxis.set_minor_formatter(mdates.DateFormatter(r"%H^{%M}"))
        for wv in self.wavelengths:
            o = self.df[self.df.wavelength==wv]
            if len(o) > 0:
                ax.semilogy(
                    o.time,
                    o.irradiance,
                    marker="o",
                    ls="None",
                    ms=0.4,
                    alpha=0.9,
                    label=rf"$\lambda_0\sim ({wv})$ nm",
                )
        ax.legend(loc=1, fontsize=7)
        ax.set_xlim(dates)
        ax.set_ylim(1e-8, 1e-2)
        if xlabel == "":
            ax.set_xticklabels([])
        return

    def save(self, figname=None, folder="tmp/figures/goes/"):
        os.makedirs(folder, exist_ok=True)
        figname = folder + figname if figname else folder + f"{self.dates[0].strftime('%Y%m%d')}.png"
        self.fig.savefig(figname, bbox_inches="tight")
        return

    def close(self):
        plt.close()
        return


class GOES(object):
    """
    This class is dedicated to load magnetic dataset from
    GOES, analyze, and plot them. The functionalities
    includes -
    1. Load GOES station info and data for a given date
    2. Plot X-ray data with additional information
    """

    def __init__(
        self,
        base,
        dates,
        station=15,
        sec_resolution=2,
    ):
        """
        Parameters:
        -----------
        base: Base folder to store data
        dates: Start and end dates
        station: GOES station ID
        sec_resolution: High resolution seconds
        """
        self.base = base
        os.makedirs(base, exist_ok=True)
        self.dates = dates
        self.station = station
        self.sec_resolution = sec_resolution
        self.low_res_data = self.fetch_1m_datasets()
        self.high_res_data = self.fetch_highres_datasets()
        return

    def __key_map__(self, res="low"):
        """
        Get keys for netcdf files
        """
        _map_ = {
            "low": {
                "15": ["A_AVG", "B_AVG"],
                "12": ["xs", "xl"],
            },
            "high": {
                "15": ["A_FLUX", "B_FLUX"],
                "12": ["xs", "xl"],
            },
        }
        return _map_[res][str(self.station)]

    def fetch_1m_datasets(self):
        """
        Fetch 1m averaged GOES-Flare data
        """
        fn = self.base + "GOES%02d.low.csv" % self.station
        if not os.path.exists(fn):
            o = pd.DataFrame()
            date = self.dates[0]
            _, d0 = calendar.monthrange(date.year, date.month)
            fname = "g{s}_xrs_1m_{y}{m}{d}_{y}{m}{d0}.nc".format(
                s=self.station,
                y=date.year,
                m="%02d" % date.month,
                d="%02d" % 1,
                d0=d0,
            )
            uri = "https://satdat.ngdc.noaa.gov/sem/goes/data/avg/{y}/{m}/goes{s}/netcdf/".format(
                s=self.station,
                y=date.year,
                m="%02d" % date.month,
            )
            url = uri + fname
            logger.info(f"File {url}")
            fname = self.base + fname
            os.system(f"wget -O {fname} {url}")
            ds = nc.Dataset(fname)
            keys = self.__key_map__()
            o["hxr"], o["sxr"] = ds.variables[keys[0]][:], ds.variables[keys[1]][:]
            tunit, tcal = (
                ds.variables["time_tag"].units,
                ds.variables["time_tag"].calendar,
            )
            o["tval"] = nc.num2date(
                ds.variables["time_tag"][:], units=tunit, calendar=tcal
            )
            o = o[(o.tval >= self.dates[0]) & (o.tval <= self.dates[1])]
            o.to_csv(fn, index=False, header=True, float_format="%g")
            os.remove(fname)
        else:
            logger.info(f"Local file {fn}")
            o = pd.read_csv(fn, parse_dates=["tval"])
        return o

    def fetch_highres_datasets(self):
        """
        Fetch high-resolution GOES-Flare data
        """
        fn = self.base + "GOES%02d.high.csv" % self.station
        if not os.path.exists(fn):
            o = pd.DataFrame()
            date = self.dates[0]
            fname = "g{s}_xrs_{r}s_{y}{m}{d}_{y}{m}{d}.nc".format(
                s=self.station,
                r=self.sec_resolution,
                y=date.year,
                m="%02d" % date.month,
                d="%02d" % date.day,
            )
            uri = "https://satdat.ngdc.noaa.gov/sem/goes/data/full/{y}/{m}/goes{s}/netcdf/".format(
                s=self.station,
                y=date.year,
                m="%02d" % date.month,
            )
            url = uri + fname
            logger.info(f"File {url}")
            fname = self.base + fname
            os.system(f"wget -O {fname} {url}")
            ds = nc.Dataset(fname)
            keys = self.__key_map__()
            o["hxr"], o["sxr"] = ds.variables[keys[0]][:], ds.variables[keys[1]][:]
            tunit, tcal = (
                ds.variables["time_tag"].units,
                ds.variables["time_tag"].calendar,
            )
            o["tval"] = nc.num2date(
                ds.variables["time_tag"][:], units=tunit, calendar=tcal
            )
            o = o[(o.tval >= self.dates[0]) & (o.tval <= self.dates[1])]
            o.to_csv(fn, index=False, header=True, float_format="%g")
            os.remove(fname)
        else:
            logger.info(f"Local file {fn}")
            o = pd.read_csv(fn, parse_dates=["tval"])
        return o

    def plot_TS_dataset(
        self,
        ax=None,
        ylim=[1e-8, 1e-3],
        comps={
            "hxr": {"color": "b", "ls": "--", "lw": 0.5},
            "sxr": {"color": "r", "ls": "--", "lw": 0.5},
        },
        high_res=False,
        xlabel="UT",
        ylabel=r"$I_{\infty}^{GOES}$, $Wm^{-2}$",
        loc=2,
        fname=None,
    ):
        """
        Overlay station data into axes
        """
        #plt.style.use(["science", "ieee"])
        import mplstyle
        if ax is None:
            fig = plt.figure(dpi=180, figsize=(5, 3))
            ax = fig.add_subplot(111)
        else:
            fig = plt.gcf()
        data = self.high_res_data if high_res else self.low_res_data
        if ylim:
            ax.set_ylim(ylim)
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        ax.set_xlim(self.dates)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        for comp in comps.keys():
            co = comps[comp]
            ax.semilogy(
                data.tval,
                data[comp],
                ls=co["ls"],
                color=co["color"],
                lw=co["lw"],
                alpha=0.9,
                label=comp.upper(),
            )
        ax.set_ylim(1e-8, 1e-2)
        ax.legend(loc=loc, fontsize=6)
        if fname:
            fig.savefig(fname, bbox_inches="tight")
        return fig, ax

    @staticmethod
    def FetchGOES(
        base,
        dates,
        station=15,
        sec_resolution=2,
    ):
        """
        Static method to call the GOES
        inevntory functions, pre-process them
        and store to object.
        """
        g = GOES(base, dates, station, sec_resolution)
        return g


# Test codes
if __name__ == "__main__":
    GOES.FetchGOES(
        "data/2015-03-11-16-20/",
        [dt.datetime(2015, 3, 11, 16), dt.datetime(2015, 3, 11, 16, 30)],
    )
    FISM.FetchFISM(
        "data/2015-03-11-16-20/",
        [dt.datetime(2015, 3, 11, 16), dt.datetime(2015, 3, 11, 16, 30)],
    )
