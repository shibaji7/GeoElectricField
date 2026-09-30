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

sys.path.extend(["py/", "py/fetch/", "py/process/", "py/plots/"])


class SDCarto(GeoAxes):
    name = "SDCarto"

    def __init__(self, *args, **kwargs):
        if "map_projection" in kwargs:
            map_projection = kwargs.pop("map_projection")
        else:
            map_projection = cartopy.crs.NorthPolarStereo()
        # first check if datetime keyword is given!
        # it should be since we need it for aacgm
        if "plot_date" in kwargs:
            self.plot_date = kwargs.pop("plot_date")
        else:
            raise TypeError(
                "need to provide a date using 'plot_date' keyword for aacgmv2 plotting"
            )
        # Now work with the coords!
        supported_coords = ["geo", "aacgmv2", "aacgmv2_mlt"]
        if "coords" in kwargs:
            self.coords = kwargs.pop("coords")
            if self.coords not in supported_coords:
                err_str = "coordinates not supported, choose from : "
                for _n, _sc in enumerate(supported_coords):
                    if _n + 1 != len(supported_coords):
                        err_str += _sc + ", "
                    else:
                        err_str += _sc
                raise TypeError(err_str)
        else:
            self.coords = "geo"
            print("coords keyword not set, setting it to aacgmv2")
        # finally, initialize te GeoAxes object
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
