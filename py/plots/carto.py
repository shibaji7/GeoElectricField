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

import cartopy
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
from cartopy.mpl.geoaxes import GeoAxes
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.projections import register_projection
from shapely.geometry import LineString, MultiLineString

sys.path.extend(["py/", "py/fetch/", "py/process/"])


class SDCarto(GeoAxes):
    name = "sdcarto"

    def __init__(self, *args, **kwargs):
        if "map_projection" in kwargs:
            map_projection = kwargs.pop("map_projection")
        else:
            map_projection = cartopy.crs.NorthPolarStereo()
        super().__init__(map_projection=map_projection, *args, **kwargs)
        return

    def overaly_coast_lakes(self, resolution="50m", color="black", **kwargs):
        """
        Overlay AACGM coastlines and lakes
        """
        kwargs["edgecolor"] = color
        kwargs["facecolor"] = "none"
        # overaly coastlines
        feature = cartopy.feature.NaturalEarthFeature(
            "physical", "coastline", resolution, **kwargs
        )
        self.add_feature(cartopy.feature.COASTLINE, **kwargs)
        self.add_feature(cartopy.feature.LAKES, **kwargs)

    def coastlines(self, resolution="50m", color="black", **kwargs):
        # details!
        kwargs["edgecolor"] = color
        kwargs["facecolor"] = "none"
        feature = cartopy.feature.NaturalEarthFeature(
            "physical", "coastline", resolution, **kwargs
        )
        return self.add_feature(feature, **kwargs)

    def mark_latitudes(self, lat_arr, lon_location=90, **kwargs):
        """
        mark the latitudes
        Write down the latitudes on the map for labeling!
        we are using this because cartopy doesn"t have a
        label by default for non-rectangular projections!
        """
        if isinstance(lat_arr, list):
            lat_arr = np.array(lat_arr)
        else:
            if not isinstance(lat_arr, np.ndarray):
                raise TypeError("lat_arr must either be a list or numpy array")
        # make an array of lon_location
        lon_location_arr = np.full(lat_arr.shape, lon_location)
        proj_xyz = self.projection.transform_points(
            cartopy.crs.Geodetic(), lon_location_arr, lat_arr
        )
        # plot the lats now!
        out_extent_lats = False
        for _np, _pro in enumerate(proj_xyz[..., :2].tolist()):
            # check if lats are out of extent! if so ignore them
            lat_lim = self.get_extent(crs=cartopy.crs.Geodetic())[2::]
            if (lat_arr[_np] >= min(lat_lim)) and (lat_arr[_np] <= max(lat_lim)):
                self.text(
                    _pro[0], _pro[1], r"$%s^{\circ}$" % str(lat_arr[_np]), **kwargs
                )
            else:
                out_extent_lats = True
        if out_extent_lats:
            print("some lats were out of extent ignored them")

    def mark_longitudes(self, lon_arr=np.arange(-180, 180, 60), **kwargs):
        """
        mark the longitudes
        Write down the longitudes on the map for labeling!
        we are using this because cartopy doesn"t have a
        label by default for non-rectangular projections!
        This is also trickier compared to latitudes!
        """
        if isinstance(lon_arr, list):
            lon_arr = np.array(lon_arr)
        else:
            if not isinstance(lon_arr, np.ndarray):
                raise TypeError("lat_arr must either be a list or numpy array")
        # get the boundaries
        [x1, y1], [x2, y2] = self.viewLim.get_points()
        bound_lim_arr = []
        right_bound = LineString(([-x1, y1], [x2, y2]))
        top_bound = LineString(([x1, -y1], [x2, y2]))
        bottom_bound = LineString(([x1, y1], [x2, -y2]))
        left_bound = LineString(([x1, y1], [-x2, y2]))
        plot_outline = MultiLineString(
            [right_bound, top_bound, bottom_bound, left_bound]
        )
        # get the plot extent, we"ll get an intersection
        # to locate the ticks!
        plot_extent = self.get_extent(cartopy.crs.Geodetic())
        line_constructor = lambda t, n, b: np.vstack(
            (np.zeros(n) + t, np.linspace(b[2], b[3], n))
        ).T
        for t in lon_arr[:-1]:
            xy = line_constructor(t, 30, plot_extent)
            # print(xy)
            proj_xyz = self.projection.transform_points(
                cartopy.crs.Geodetic(), xy[:, 0], xy[:, 1]
            )
            xyt = proj_xyz[..., :2]
            ls = LineString(xyt.tolist())
            locs = plot_outline.intersection(ls)
            if not locs:
                continue
            # we need to get the alignment right
            # so get the boundary closest to the label
            # and plot it!
            closest_bound = min(
                [
                    right_bound.distance(locs),
                    top_bound.distance(locs),
                    bottom_bound.distance(locs),
                    left_bound.distance(locs),
                ]
            )
            if closest_bound == right_bound.distance(locs):
                ha = "left"
                va = "top"
            elif closest_bound == top_bound.distance(locs):
                ha = "left"
                va = "bottom"
            elif closest_bound == bottom_bound.distance(locs):
                ha = "left"
                va = "top"
            else:
                ha = "right"
                va = "top"
            #             if self.coords == "aacgmv2_mlt":
            #                 marker_text = str(int(t/15.))
            #             else:
            #                 marker_text = str(t)
            marker_text = r"$%s^{\circ}$" % str(t)
            self.text(
                locs.bounds[0] + 0.02 * locs.bounds[0],
                locs.bounds[1] + 0.02 * locs.bounds[1],
                marker_text,
                ha=ha,
                va=va,
                **kwargs
            )
        return


register_projection(SDCarto)


class MapPlot(object):
    """
    Plot stations / data
    """

    def __init__(self, hemi="north", coords="", dtime=None):
        self.hemi = hemi
        self.coords = coords
        self.dtime = dtime
        self.ini_figure()
        return

    def ini_figure(self):
        """
        Instatitate figure and axes labels
        """
        proj = (
            cartopy.crs.NorthPolarStereo(-120)
            if self.hemi == "north"
            else cartopy.crs.SouthPolarStereo(-120)
        )

        self.fig = plt.figure(dpi=300, figsize=(3, 3))
        self.ax = self.fig.add_subplot(111, projection="sdcarto", map_projection=proj)
        self.ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        if self.hemi == "north":
            self.ax.set_extent([-180, 180, 30, 90], crs=cartopy.crs.PlateCarree())
        else:
            self.ax.set_extent([-180, 180, -90, -30], crs=cartopy.crs.PlateCarree())
        plt_lons = np.arange(-180, 181, 30)
        mark_lons = np.arange(-180, 181, 30)
        plt_lats = (
            np.arange(30, 90, 10) if self.hemi == "north" else np.arange(-90, -30, 10)
        )
        gl = self.ax.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=0.5)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        self.ax.mark_latitudes(plt_lats, fontsize="xx-small", color="darkblue")
        self.ax.mark_longitudes(plt_lons, fontsize="xx-small", color="darkblue")
        #         #self.ax.text(0.5, 0.95, self.date_string(), ha="center", va="center",
        #         #             transform=self.ax.transAxes, fontsize="medium")
        self.proj = proj
        self.geo = cartopy.crs.Geodetic()
        self.ax.text(
            -0.02,
            0.99,
            "Coords: Geo",
            ha="center",
            va="top",
            transform=self.ax.transAxes,
            fontsize="xx-small",
            rotation=90,
        )
        return

    def put_stations(
        self,
        lat,
        lon,
        name=None,
        marker="D",
        zorder=2,
        markerColor="m",
        markerSize=2,
        fontSize="xx-small",
        font_color="m",
        xOffset=-5,
        yOffset=-1.5,
    ):
        self.ax.scatter(
            [lon],
            [lat],
            s=markerSize,
            marker=marker,
            color=markerColor,
            zorder=zorder,
            transform=self.geo,
            lw=0.8,
            alpha=0.4,
        )

        if name:
            x, y = self.projection.transform_point(
                lon + xOffset, lat + yOffset, src_crs=self.proj
            )
            self.ax.text(
                x,
                y,
                name.upper(),
                ha="center",
                va="center",
                transform=self.proj,
                fontdict={"color": font_color, "size": fontSize},
                alpha=0.4,
            )
        return

    def draw_GC(self, lat1, lon1, lat2, lon2, color="k"):
        self.ax.plot(
            [lon1, lon2], [lat1, lat2], color=color, transform=self.geo, ls="-", lw=0.3
        )
        return

    def save(self, fname):
        self.fig.savefig(fname, bbox_inches="tight")
        return


if __name__ == "__main__":
    mp = MapPlot()
    mp.save("out.png")
