#!/usr/bin/env python3
"""Plot smoothed 1-degree binned EMTF electric-field magnitudes.

This diagnostic is intentionally separate from ``event_study.py``.  It shows
the spatial field used for threshold-contour inspection: station values are
reduced to the maximum |E| in each 1-degree CONUS bin and then smoothed.

Example
-------
    python smooth_efield_bins.py \
        --config config/event_study/2005-09-07-17-40.toml
"""

import argparse
import datetime as dt
import os
import sys

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import toml
from matplotlib.colors import LogNorm
from scipy.ndimage import (
    binary_closing,
    binary_dilation,
    binary_fill_holes,
    gaussian_filter,
    label,
)


PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_EMTF = os.path.join(PROJECT_ROOT, ".data", "EMTF")
CONUS_BOUNDS = (-125.0, -66.9, 24.4, 49.4)

sys.path.extend(
    [
        os.path.join(PROJECT_ROOT, "py"),
        os.path.join(PROJECT_ROOT, "py", "fetch"),
        os.path.join(PROJECT_ROOT, "py", "process"),
        os.path.join(PROJECT_ROOT, "py", "plots"),
    ]
)

from efield import _emtf_site_fields, _is_conus_location  # noqa: E402
from get_mag_data import SuperMAG  # noqa: E402


def parse_datetime(value):
    return dt.datetime.fromisoformat(value.replace("Z", "+00:00")).replace(
        tzinfo=None
    )


def parse_case_time(value, event):
    """Parse an ISO datetime or an HH:MM value on the event day."""
    value = str(value)
    if len(value) <= 5 and ":" in value:
        hour, minute = (int(part) for part in value.split(":")[:2])
        return event.replace(hour=hour, minute=minute, second=0, microsecond=0)
    return parse_datetime(value)


def event_base(event):
    return os.path.join(PROJECT_ROOT, "tmp", event.strftime("%Y-%m-%d-%H-%M"))


def cache_covers_window(filename, dates):
    if not os.path.exists(filename):
        return False
    try:
        values = pd.read_csv(filename, usecols=["tval"], parse_dates=["tval"])["tval"]
    except (OSError, ValueError, pd.errors.EmptyDataError):
        return False
    values = values.dropna()
    return not values.empty and values.min() <= dates[0] and values.max() >= dates[1]


def bin_efield_statistics(fields, timestamp, bin_size, excluded,
                          statistic="percentile", percentile=75.0):
    """Bin site magnitudes and calculate one robust statistic per cell."""
    west, east, south, north = CONUS_BOUNDS
    nlon = int(np.ceil((east - west) / bin_size))
    nlat = int(np.ceil((north - south) / bin_size))
    values_by_bin = {}
    station_points = []

    for site, _, data in fields:
        name = str(getattr(site, "name", "")).upper()
        if name in excluded or data.empty:
            continue
        delta = (data["tval"] - timestamp).abs()
        row = data.loc[delta.idxmin()]
        if abs((row["tval"] - timestamp).total_seconds()) > 90:
            continue
        east_value = float(row["induced_E_east"])
        north_value = float(row["induced_E_north"])
        lon = (float(site.longitude) + 180.0) % 360.0 - 180.0
        lat = float(site.latitude)
        if not np.isfinite(east_value + north_value):
            continue
        if not _is_conus_location(lon, lat):
            continue
        magnitude = float(np.hypot(east_value, north_value))
        if not np.isfinite(magnitude) or magnitude <= 0.0:
            continue
        lon_index = min(max(int((lon - west) / bin_size), 0), nlon - 1)
        lat_index = min(max(int((lat - south) / bin_size), 0), nlat - 1)
        key = (lat_index, lon_index)
        values_by_bin.setdefault(key, []).append(magnitude)
        station_points.append((lon, lat, magnitude))

    grid = np.full((nlat, nlon), np.nan, dtype=float)
    counts = np.zeros((nlat, nlon), dtype=int)
    for (lat_index, lon_index), bin_values in values_by_bin.items():
        counts[lat_index, lon_index] = len(bin_values)
        if statistic == "maximum":
            value = np.max(bin_values)
        elif statistic == "median":
            value = np.median(bin_values)
        elif statistic == "percentile":
            value = np.percentile(bin_values, percentile)
        else:
            raise ValueError(f"Unsupported bin statistic: {statistic}")
        grid[lat_index, lon_index] = value
    lon_edges = west + np.arange(nlon + 1) * bin_size
    lat_edges = south + np.arange(nlat + 1) * bin_size
    return lon_edges, lat_edges, grid, counts, station_points


def smooth_grid(grid, sigma):
    """Smooth a sparse positive grid without treating missing bins as zero."""
    valid = np.isfinite(grid)
    values = np.where(valid, grid, 0.0)
    weights = valid.astype(float)
    smoothed_values = gaussian_filter(values, sigma=sigma)
    smoothed_weights = gaussian_filter(weights, sigma=sigma)
    with np.errstate(divide="ignore", invalid="ignore"):
        smoothed = smoothed_values / smoothed_weights
    smoothed[smoothed_weights < 0.15] = np.nan
    return smoothed


def contour_support(grid, counts, smoothed, threshold, minimum_stations,
                    minimum_neighbor_bins):
    """Return a binary mask for connected above-threshold regions.

    The region mask is thresholded after smoothing, but only where the
    original bins are sufficiently sampled or immediately supported by a
    sampled bin.  Boundaries are extracted from this binary mask rather than
    from an interpolated scalar-field isoline.  This prevents open, broken
    contour fragments in sparse parts of the grid.
    """
    sampled = counts >= int(minimum_stations)
    sampled_support = binary_dilation(
        sampled, structure=np.ones((3, 3), dtype=bool), iterations=1
    )
    qualifying = (
        np.isfinite(smoothed)
        & (smoothed >= float(threshold))
        & sampled_support
    )

    # Close one-cell gaps and fill holes so each retained region has a clean
    # binary footprint.  Components are filtered after closing so tiny noise
    # cannot become a contour on its own.
    qualifying = binary_closing(
        qualifying, structure=np.ones((3, 3), dtype=bool), iterations=1
    )
    qualifying = binary_fill_holes(qualifying)
    labels, number = label(qualifying, structure=np.ones((3, 3), dtype=bool))
    retained = np.zeros_like(qualifying, dtype=bool)
    for component in range(1, number + 1):
        component_mask = labels == component
        if int(component_mask.sum()) >= int(minimum_neighbor_bins):
            retained |= component_mask
    if not retained.any():
        return np.zeros_like(qualifying, dtype=bool), 0
    return retained, int(retained.sum())


def main():
    parser = argparse.ArgumentParser(
        description="Plot smoothed 1-degree binned EMTF |E| fields."
    )
    parser.add_argument("--config", required=True, help="Event-study TOML file.")
    parser.add_argument("--bin-size", type=float, default=None)
    parser.add_argument("--sigma", type=float, default=1.5)
    parser.add_argument("--statistic", choices=["maximum", "median", "percentile"],
                        default=None)
    parser.add_argument("--percentile", type=float, default=None)
    parser.add_argument("--minimum-stations", type=int, default=None)
    parser.add_argument("--minimum-neighbor-bins", type=int, default=None)
    parser.add_argument("--threshold", type=float, default=None,
                        help="Contour threshold in mV/km.")
    parser.add_argument("--output", default=None,
                        help="Output PNG path; defaults to the event pub_fig folder.")
    args = parser.parse_args()

    config = toml.load(args.config)
    event_config = config["event"]
    map_config = config.get("map", {})
    plot_config = config.get("plot", {})
    event = parse_datetime(event_config["time"])
    dates = [
        parse_datetime(event_config["start"]),
        parse_datetime(event_config["end"]),
    ]
    baseline_hours = float(plot_config.get("case_baseline_hours", 2.0))
    case_xlim = plot_config.get("case_xlim", [event_config["start"], event_config["end"]])
    reference_time = parse_case_time(case_xlim[0], event)
    magnetic_dates = [
        min(dates[0], reference_time - dt.timedelta(hours=baseline_hours)),
        dates[1],
    ]
    base = event_base(event)
    raw_baseline = "none" if event_config.get("mag_baseline") == "none" else event_config.get("mag_baseline")
    cache_name = "supermag-raw.csv" if raw_baseline == "none" else "supermag.csv"
    magnetic_base = base + os.sep
    if not cache_covers_window(os.path.join(base, cache_name), magnetic_dates):
        magnetic_base = os.path.join(base, "case-study-extended") + os.sep

    smag = SuperMAG(
        magnetic_base,
        magnetic_dates,
        uid=event_config.get("supermag_user", "shibaji7"),
        baseline=raw_baseline,
    )
    fields = _emtf_site_fields(
        smag,
        event_config.get("emtf_data_dir", DEFAULT_EMTF),
        dT=60,
        bounds=CONUS_BOUNDS,
    )
    time_strings = map_config.get("case_map_times", ["17:20", "17:40"])
    timestamps = [
        event.replace(
            hour=int(value.split(":")[0]),
            minute=int(value.split(":")[1]),
            second=0,
            microsecond=0,
        )
        for value in time_strings[:2]
    ]
    bin_size = float(args.bin_size or plot_config.get("efield_contour_bin_size_deg", 1.0))
    threshold = float(args.threshold or plot_config.get("efield_contour_threshold_mV_per_km", 250.0))
    excluded = {str(name).upper() for name in plot_config.get("case_efield_exclude_sites", [])}
    statistic = args.statistic or plot_config.get("efield_bin_statistic", "percentile")
    percentile = float(args.percentile or plot_config.get("efield_bin_percentile", 75.0))
    minimum_stations = int(args.minimum_stations or plot_config.get("efield_min_stations_per_bin", 3))
    minimum_neighbor_bins = int(args.minimum_neighbor_bins or plot_config.get("efield_min_neighbor_bins", 2))

    binned = [
        bin_efield_statistics(fields, timestamp, bin_size, excluded, statistic, percentile)
        for timestamp in timestamps
    ]
    smoothed = [smooth_grid(item[2], args.sigma) for item in binned]
    contour_fields = [
        contour_support(item[2], item[3], values, threshold,
                        minimum_stations, minimum_neighbor_bins)
        for item, values in zip(binned, smoothed)
    ]
    finite_values = np.concatenate([
        values[np.isfinite(values)] for values in smoothed
        if np.isfinite(values).any()
    ])
    if finite_values.size == 0:
        raise RuntimeError("No finite smoothed E-field values were produced.")
    lower = max(1.0, float(np.nanpercentile(finite_values, 1)))
    upper = max(lower * 10.0, float(np.nanpercentile(finite_values, 99)))

    projection = ccrs.LambertConformal(
        central_longitude=-96.0,
        central_latitude=39.0,
        standard_parallels=(29.0, 45.0),
    )
    fig, axes = plt.subplots(
        1, len(timestamps), figsize=(9.0, 4.0),
        subplot_kw={"projection": projection},
        squeeze=False,
    )
    axes = axes[0]
    mesh = None
    for index, (axis, timestamp, item, values, contour_item) in enumerate(
        zip(axes, timestamps, binned, smoothed, contour_fields)
    ):
        lon_edges, lat_edges, _, _, station_points = item
        contour_mask, active_bins = contour_item
        mesh = axis.pcolormesh(
            lon_edges,
            lat_edges,
            values,
            transform=ccrs.PlateCarree(),
            cmap="viridis",
            norm=LogNorm(vmin=lower, vmax=upper),
            shading="auto",
        )
        if contour_mask.any():
            # Pad with a below-threshold border so even regions touching the
            # sampled grid edge produce closed boundaries.
            padded_mask = np.pad(contour_mask.astype(float), 1, constant_values=0.0)
            lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
            lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
            bin_size_x = float(np.median(np.diff(lon_centers)))
            bin_size_y = float(np.median(np.diff(lat_centers)))
            padded_lon = np.concatenate(
                ([lon_centers[0] - bin_size_x], lon_centers,
                 [lon_centers[-1] + bin_size_x])
            )
            padded_lat = np.concatenate(
                ([lat_centers[0] - bin_size_y], lat_centers,
                 [lat_centers[-1] + bin_size_y])
            )
            axis.contour(
                padded_lon,
                padded_lat,
                padded_mask,
                levels=[0.5],
                colors="#0066FF",
                linewidths=0.8,
                transform=ccrs.PlateCarree(),
            )
        axis.set_extent(CONUS_BOUNDS, crs=ccrs.PlateCarree())
        axis.coastlines(resolution="50m", color="0.25", linewidth=0.6)
        axis.add_feature(cfeature.BORDERS, edgecolor="0.25", linewidth=0.4)
        axis.gridlines(draw_labels=False, linewidth=0.25, color="0.75", linestyle=":")
        axis.text(0.02, 0.95, f"({chr(97 + index)}) {timestamp:%H:%M} UT",
                  transform=axis.transAxes, va="top", fontsize=11)
        axis.text(0.02, 0.04, f"N={len(station_points)} stations; bins={active_bins}",
                  transform=axis.transAxes, fontsize=8)

    colorbar = fig.colorbar(mesh, ax=axes, orientation="vertical", fraction=0.035, pad=0.02)
    colorbar.set_label(r"Smoothed binned $|\mathcal{E}|$ (mV/km)")
    fig.suptitle(
        f"1° binned and smoothed E-field; binary contour threshold = "
        f"{threshold:g} mV/km",
        fontsize=11,
    )
    fig.subplots_adjust(left=0.02, right=0.91, bottom=0.04, top=0.88, wspace=0.03)
    output = args.output or os.path.join(
        PROJECT_ROOT, "pub_fig", event.strftime("%Y-%m-%d-%H-%M"),
        "smoothed-efield-bins.png",
    )
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
