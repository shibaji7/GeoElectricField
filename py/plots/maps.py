#!/usr/bin/env python

"""carto.py: module is dedicated to for carto plots."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys
sys.path.extend(["py/", "py/fetch/", "py/process/", "py/process/plots/"])
import mplstyle
import matplotlib.pyplot as plt
import cartopy
import numpy as np
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from carto import SDCarto

sys.path.extend(["py/", "py/fetch/", "py/process/", "py/plots/"])

class Map(object):

    def __init__(
        self, date,
        fig_title=None,
        nrows=1,
        ncols=1,
        coord="geo",
    ):
        self.date = date
        self.nrows, self.ncols = nrows, ncols
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(3 * ncols, 3 * nrows), dpi=300)
        self.coord = coord
        self.ax = self.add_axes()
        # self.ax.text(
        #     0.1, 1.05,
        #     f"{fig_title}: {self.date_string()}"
        #     if fig_title
        #     else f"{self.date_string()}",
        #     ha="left",
        #     fontweight="bold",
        #     fontsize=10,
        #     transform=self.ax.transAxes,
        # )
        return

    def add_axes(self):
        """
        Instatitate figure and axes labels
        """
        from carto import SDCarto
        self._num_subplots_created += 1
        #self.proj = cartopy.crs.NorthPolarStereo(central_longitude=-90.0)
        self.proj = cartopy.crs.Stereographic(central_longitude=-90.0, central_latitude=45.0)
        #proj = cartopy.crs.PlateCarree(central_longitude=-90.0)
        ax = self.fig.add_subplot(
            100 * self.nrows + 10 * self.ncols + self._num_subplots_created,
            projection="SDCarto",
            map_projection=self.proj,
            coords=self.coord,
            plot_date=self.date,
        )
        ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        self.lon_lat_bb = [-130, -70, 30, 80]
        ax.set_extent(self.lon_lat_bb, crs=cartopy.crs.PlateCarree())
        plt_lons = np.arange(-180, 181, 15)
        mark_lons = np.arange(-180, 181, 30)
        plt_lats = np.arange(30, 80, 10)
        gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=0.5)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        #ax.mark_latitudes(plt_lats, fontsize="xx-small", color="k")
        #ax.mark_longitudes(mark_lons, fontsize="xx-small", color="k")
        self.geo = cartopy.crs.PlateCarree()
        # ax.text(
        #     -0.02,
        #     0.99,
        #     "Coord: Geo",
        #     ha="center",
        #     va="top",
        #     transform=ax.transAxes,
        #     fontsize=10,
        #     rotation=90,
        # )
        from cartopy.feature.nightshade import Nightshade
        ax.add_feature(Nightshade(self.date, alpha=0.2))
        return ax

    def date_string(self, label_style="web"):
        # Set the date and time formats
        dfmt = "%d %b %Y" if label_style == "web" else "%d %b %Y,"
        tfmt = "%H:%M"
        stime = self.date
        date_str = "{:{dd} {tt}} UT".format(stime, dd=dfmt, tt=tfmt)
        return date_str

    def plot_instrument(
        self, 
        instrument,
        station_point=dict(
            ms=2, marker="D",
            color="r", fontSize=4, 
            xOff=-1.5, yOff=-1.5
        ),
        quiver_style=dict(
            scale=300,
            headaxislength=0,
            width=0.1,
            scale_units="inches",
        )
    ):
        if (self.lon_lat_bb[2] <= instrument["lat"] <= self.lon_lat_bb[3]) and\
         (self.lon_lat_bb[0] <= instrument["lon"] <= self.lon_lat_bb[1]):
            self.ax.scatter(
                [instrument["lon"]], [instrument["lat"]], transform=self.geo, 
                s=station_point["ms"], marker=station_point["marker"], 
                color=station_point["color"]
            )
            x, y = self.proj.transform_point(
                instrument["lon"]+station_point["xOff"], 
                instrument["lat"]+station_point["yOff"], 
                self.geo
            )
            self.ax.text(
                x, y, instrument["code"], 
                ha="center", va="center", 
                transform=self.proj,
                fontdict={"color":station_point["color"], "size":station_point["fontSize"]}, 
                alpha=0.7
            )
            xy = self.proj.transform_points(
                self.geo, 
                np.array([instrument["lon"]]), 
                np.array([instrument["lat"]]),
            )
            x, y = xy[:, 0], xy[:, 1]
            ql = self.ax.quiver(
                x,
                y,
                instrument["E_geo"],
                instrument["N_geo"],
                scale=quiver_style["scale"],
                headaxislength=quiver_style["headaxislength"],
                scale_units=quiver_style["scale_units"],
                color="blue",
                linewidth=0.4,
            )
            self.ax.quiverkey(
                ql,
                0.8,
                0.9,
                50,
                r"$\vec{B}$ [50 nT]",
                labelpos="S",
                transform=self.ax.transAxes,
                color="k",
                fontproperties={"size": 8},
            )
        return

    def save(self, filepath):
        self.fig.savefig(filepath, facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return