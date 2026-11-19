"""Induced electric-field plots for event studies."""

import glob
import datetime as dt
import csv
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
import matplotlib.dates as mdates
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.interpolate import griddata
from scipy.ndimage import (
    binary_closing,
    binary_dilation,
    binary_fill_holes,
    gaussian_filter,
    label,
)
from scipy.spatial import cKDTree

import bezpy
from utils import calculate_GC_distance, detrend_magnetic_field
from utils import compute_Egeo_from_Bgeo
from maps import Map


EMTF_DATA_FOLDER = ".data/EMTF"
CONUS_BOUNDS = (-125.0, -66.9, 24.4, 49.4)
# A compact geographic mask used in addition to the rectangular map extent.
# The rectangular extent alone includes southern Canada and northern Mexico.
CONUS_POLYGON = np.asarray([
    (-124.8, 49.0), (-95.0, 49.0), (-89.5, 48.0), (-84.5, 46.0),
    (-82.0, 43.0), (-79.0, 43.0), (-75.0, 45.0), (-67.0, 45.0),
    (-67.0, 25.0), (-80.0, 25.0), (-82.0, 27.0), (-81.0, 30.0),
    (-89.0, 29.0), (-94.0, 29.0), (-97.0, 26.0), (-106.0, 31.0),
    (-114.0, 32.0), (-117.0, 32.0), (-120.0, 34.0), (-124.0, 40.0),
    (-125.0, 46.0),
])
_CONUS_PATH = Path(CONUS_POLYGON)


def _is_conus_location(longitude, latitude):
    """Return whether a longitude/latitude point lies in the CONUS mask."""
    return bool(_CONUS_PATH.contains_point((float(longitude), float(latitude))))


def _clip_contour_to_conus(contour_set, ax):
    """Clip a geographic contour collection to the CONUS polygon."""
    clip_path = PathPatch(
        _CONUS_PATH,
        transform=ccrs.PlateCarree()._as_mpl_transform(ax),
    )
    for collection in contour_set.collections:
        collection.set_clip_path(clip_path)


def _baseline_correct_component(data, component, baseline_start, reference_time):
    """Remove a linear pre-event trend and reference a component at event start."""
    times = data["tval"]
    values = np.asarray(data[f"{component}_geo"], dtype=float)
    seconds = (times - baseline_start).dt.total_seconds().to_numpy(dtype=float)
    baseline_mask = (
        (times >= baseline_start)
        & (times < reference_time)
        & np.isfinite(values)
    ).to_numpy()
    if baseline_mask.sum() < 2:
        raise ValueError(
            f"At least two finite {component}_geo samples are required between "
            f"{baseline_start:%H:%M} and {reference_time:%H:%M} UT."
        )

    coefficients = np.polyfit(seconds[baseline_mask], values[baseline_mask], 1)
    corrected = values - np.polyval(coefficients, seconds)
    reference_mask = (
        (times >= reference_time) & np.isfinite(corrected)
    ).to_numpy()
    if not reference_mask.any():
        raise ValueError(
            f"No finite {component}_geo sample was found at or after "
            f"{reference_time:%H:%M} UT."
        )
    reference_index = np.flatnonzero(reference_mask)[0]
    corrected -= corrected[reference_index]
    return corrected


def _case_baseline_corrected_components(
    data, baseline_start, reference_time
):
    """Return N/E components detrended before event start and zeroed at it."""
    return (
        _baseline_correct_component(data, "N", baseline_start, reference_time),
        _baseline_correct_component(data, "E", baseline_start, reference_time),
    )


def _reference_series_at_time(values, times, reference_time):
    """Subtract the first finite series value at or after ``reference_time``."""
    referenced = np.asarray(values, dtype=float).copy()
    reference_mask = (
        (times >= reference_time) & np.isfinite(referenced)
    ).to_numpy()
    if not reference_mask.any():
        raise ValueError(
            f"No finite series sample was found at or after "
            f"{reference_time:%H:%M} UT."
        )
    referenced -= referenced[np.flatnonzero(reference_mask)[0]]
    return referenced


def _bin_efield_statistics(
    longitudes,
    latitudes,
    magnitudes_mV_per_km,
    bounds=CONUS_BOUNDS,
    bin_size_deg=1.0,
    statistic="percentile",
    percentile=75.0,
):
    """Bin station magnitudes and calculate one robust value per cell."""
    west, east, south, north = bounds
    bin_size = float(bin_size_deg)
    nlon = int(np.ceil((east - west) / bin_size))
    nlat = int(np.ceil((north - south) / bin_size))
    longitude = np.asarray(longitudes, dtype=float)
    latitude = np.asarray(latitudes, dtype=float)
    magnitude = np.asarray(magnitudes_mV_per_km, dtype=float)
    valid = np.isfinite(longitude) & np.isfinite(latitude) & np.isfinite(magnitude)
    grid = np.full((nlat, nlon), np.nan, dtype=float)
    counts = np.zeros((nlat, nlon), dtype=int)
    values_by_bin = {}
    if not valid.any():
        lon_edges = west + np.arange(nlon + 1) * bin_size
        lat_edges = south + np.arange(nlat + 1) * bin_size
        return lon_edges, lat_edges, grid, counts

    for site_lon, site_lat, site_value in zip(
        longitude[valid], latitude[valid], magnitude[valid]
    ):
        if not _is_conus_location(site_lon, site_lat):
            continue
        lon_index = min(max(int((site_lon - west) / bin_size), 0), nlon - 1)
        lat_index = min(max(int((site_lat - south) / bin_size), 0), nlat - 1)
        values_by_bin.setdefault((lat_index, lon_index), []).append(float(site_value))

    for (lat_index, lon_index), values in values_by_bin.items():
        counts[lat_index, lon_index] = len(values)
        if statistic == "maximum":
            cell_value = np.max(values)
        elif statistic == "median":
            cell_value = np.median(values)
        elif statistic == "percentile":
            cell_value = np.percentile(values, float(percentile))
        else:
            raise ValueError(f"Unsupported E-field bin statistic: {statistic}")
        grid[lat_index, lon_index] = float(cell_value)

    lon_edges = west + np.arange(nlon + 1) * bin_size
    lat_edges = south + np.arange(nlat + 1) * bin_size
    return lon_edges, lat_edges, grid, counts


def _smooth_sparse_efield_grid(grid, sigma=1.5, minimum_weight=0.15):
    """Gaussian-smooth sparse bins without treating missing cells as zero."""
    values = np.asarray(grid, dtype=float)
    valid = np.isfinite(values)
    filled = np.where(valid, values, 0.0)
    weights = valid.astype(float)
    smoothed_values = gaussian_filter(filled, sigma=float(sigma))
    smoothed_weights = gaussian_filter(weights, sigma=float(sigma))
    with np.errstate(divide="ignore", invalid="ignore"):
        smoothed = smoothed_values / smoothed_weights
    smoothed[smoothed_weights < float(minimum_weight)] = np.nan
    return smoothed


def _binary_efield_contour_mask(
    grid,
    counts,
    smoothed,
    threshold,
    minimum_stations=3,
    minimum_neighbor_bins=2,
):
    """Build connected binary regions for closed E-field boundaries."""
    sampled = np.asarray(counts) >= int(minimum_stations)
    sampled_support = binary_dilation(
        sampled, structure=np.ones((3, 3), dtype=bool), iterations=1
    )
    qualifying = (
        np.isfinite(smoothed)
        & (smoothed >= float(threshold))
        & sampled_support
    )
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
    return retained


def plot_efield_threshold_contour(
    ax,
    longitudes,
    latitudes,
    magnitudes_mV_per_km,
    threshold_mV_per_km,
    bounds=CONUS_BOUNDS,
    bin_size_deg=1.0,
    smoothing_sigma=1.5,
    statistic="percentile",
    percentile=75.0,
    minimum_stations=3,
    minimum_neighbor_bins=2,
    color="#0066FF",
):
    """Draw closed boundaries of connected, above-threshold E-field regions."""
    try:
        threshold = float(threshold_mV_per_km)
        bin_size = float(bin_size_deg)
    except (TypeError, ValueError):
        return None
    if (
        not np.isfinite(threshold)
        or threshold <= 0.0
        or not np.isfinite(bin_size)
        or bin_size <= 0.0
    ):
        return None

    longitude = np.asarray(longitudes, dtype=float)
    latitude = np.asarray(latitudes, dtype=float)
    magnitude = np.asarray(magnitudes_mV_per_km, dtype=float)
    valid = np.isfinite(longitude) & np.isfinite(latitude) & np.isfinite(magnitude)
    if valid.sum() < 4:
        return None

    west, east, south, north = bounds
    lon_edges, lat_edges, binned, counts = _bin_efield_statistics(
        longitude[valid], latitude[valid], magnitude[valid],
        bounds=bounds, bin_size_deg=bin_size,
        statistic=statistic, percentile=percentile,
    )
    if not np.isfinite(binned).any():
        return None
    smoothed = _smooth_sparse_efield_grid(binned, sigma=smoothing_sigma)
    if not np.isfinite(smoothed).any():
        return None
    contour_mask = _binary_efield_contour_mask(
        binned,
        counts,
        smoothed,
        threshold,
        minimum_stations=minimum_stations,
        minimum_neighbor_bins=minimum_neighbor_bins,
    )
    if not contour_mask.any():
        return None

    lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    padded_mask = np.pad(contour_mask.astype(float), 1, constant_values=0.0)
    padded_lon = np.concatenate(
        ([lon_centers[0] - bin_size], lon_centers,
         [lon_centers[-1] + bin_size])
    )
    padded_lat = np.concatenate(
        ([lat_centers[0] - bin_size], lat_centers,
         [lat_centers[-1] + bin_size])
    )
    try:
        contour_set = ax.contour(
            padded_lon,
            padded_lat,
            padded_mask,
            levels=[0.5],
            colors=[color],
            linewidths=0.8,
            transform=ccrs.PlateCarree(),
            zorder=3,
        )
        _clip_contour_to_conus(contour_set, ax)
    except (ValueError, RuntimeError):
        return None
    # Keep the complete threshold contour set, including paths that intersect
    # the CONUS map boundary. This matches the standalone smoothed-bin figure.
    return contour_set


def _valid_vertical_component(data):
    """Return True only when the complete station Z series is usable."""
    if "Z_geo" not in data.columns:
        return False
    vertical = np.asarray(data["Z_geo"], dtype=float)
    return vertical.size > 1 and np.isfinite(vertical).all()


def _load_emtf_sites(mt_data_folder):
    """Load one EMTF response site for every available EMTF station."""
    sites = {}
    pattern = os.path.join(mt_data_folder, "MT_TF_*", "*.xml")
    for filename in glob.glob(pattern):
        try:
            site = bezpy.mt.read_xml(filename)
        except Exception:
            continue
        name = getattr(site, "name", None)
        latitude = getattr(site, "latitude", None)
        longitude = getattr(site, "longitude", None)
        if name and latitude is not None and longitude is not None:
            sites.setdefault(name, site)
    return list(sites.values())


def _emtf_site_fields(smag, mt_data_folder, dT=60, bounds=None):
    """Calculate E fields at EMTF sites using nearest SuperMAG stations."""
    station_data = {}
    for station_name in smag.sm_data.iaga.dropna().unique():
        data = smag.sm_data[smag.sm_data.iaga == station_name].copy()
        if data.empty or not _valid_vertical_component(data):
            continue
        row = data.iloc[0]
        station_data[station_name] = (data, float(row.glat), float(row.glon))

    sites = []
    for site in _load_emtf_sites(mt_data_folder):
        site_longitude = float(site.longitude)
        site_latitude = float(site.latitude)
        if bounds is not None:
            west, east, south, north = bounds
            if not (west <= site_longitude <= east and south <= site_latitude <= north):
                continue
            if bounds == CONUS_BOUNDS and not _is_conus_location(
                site_longitude, site_latitude
            ):
                continue
        nearest = min(
            station_data,
            key=lambda station: calculate_GC_distance(
                float(site.latitude),
                float(site.longitude),
                station_data[station][1],
                station_data[station][2],
            ),
            default=None,
        )
        if nearest is None:
            continue
        data, _, _ = station_data[nearest]
        try:
            # EMTF expects the horizontal magnetic components in N/E order.
            magnetic_north = detrend_magnetic_field(data["N_geo"].to_numpy(), dT=dT)
            magnetic_east = detrend_magnetic_field(data["E_geo"].to_numpy(), dT=dT)
            e_north, e_east = site.convolve_fft(
                magnetic_north,
                magnetic_east,
                dt=dT,
            )
        except Exception:
            continue
        if len(e_north) != len(data) or len(e_east) != len(data):
            continue
        data = data.copy()
        data["induced_E_north"] = np.asarray(e_north)
        data["induced_E_east"] = np.asarray(e_east)
        sites.append((site, nearest, data))
    return sites


def plot_efield(smag, event, dates, stations, output, dT=60,
                mt_data_folder=EMTF_DATA_FOLDER):
    """Plot induced geographic electric fields for the requested stations.

    The station order is preserved.  Stations without a matching EMTF
    magnetotelluric response are skipped because an electric field cannot be
    calculated for them.
    """
    datasets = []
    for station_name in stations:
        data = smag.sm_data[smag.sm_data.iaga == station_name].copy()
        if data.empty or not _valid_vertical_component(data):
            continue

        row = data.iloc[0]
        magnetic_north = detrend_magnetic_field(data["N_geo"].to_numpy(), dT=dT)
        magnetic_east = detrend_magnetic_field(data["E_geo"].to_numpy(), dT=dT)
        north, east, site, _ = compute_Egeo_from_Bgeo(
            magnetic_north,
            magnetic_east,
            glat=row.glat,
            glon=row.glon,
            dT=dT,
            MT_data_folder=mt_data_folder,
        )
        if site is None or len(north) != len(data) or len(east) != len(data):
            continue

        data["E_north"] = north
        data["E_east"] = east
        datasets.append((station_name, data, site))

    if not datasets:
        raise RuntimeError("No requested stations have usable induced E-field data.")

    fig = plt.figure(figsize=(5, 2.2 * len(datasets)), dpi=240)
    for index, (station_name, data, site) in enumerate(datasets):
        ax = fig.add_subplot(len(datasets), 1, index + 1)
        ax.set_ylabel(r"$E_{geo}$ (mV/km)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24)))
        ax.set_xlim(dates)
        ax.plot(data.tval, data["E_north"], color="r", lw=0.6, label=r"$E_N$")
        ax.plot(data.tval, data["E_east"], color="b", lw=0.6, label=r"$E_E$")
        ax.axvline(event, color="k", lw=0.5)
        ax.set_ylim(-50, 100)
        ax.text(
            0.03,
            0.82,
            f"({chr(97 + index)}) {station_name} / {getattr(site, 'name', station_name)}",
            transform=ax.transAxes,
            fontsize=8,
        )
        if index == 0:
            ax.legend(loc=1, fontsize=8)
        if index == len(datasets) - 1:
            ax.set_xlabel("Time (UT)")
        else:
            ax.set_xticklabels([])

    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_efield_map(smag, map_times, region, output, dT=60,
                    mt_data_folder=EMTF_DATA_FOLDER, quiver_scale=None,
                    quiver_limit=None):
    """Plot induced E-field vectors at the same times as the magnetic map."""
    if not map_times:
        raise ValueError("At least one map time is required for maps-efield.png.")
    if len(map_times) > 6:
        raise ValueError("maps-efield.png supports at most six map times.")

    fields = _emtf_site_fields(
        smag,
        mt_data_folder,
        dT=dT,
        bounds=CONUS_BOUNDS,
    )

    if not fields:
        raise RuntimeError("No stations have usable induced E-field map data.")

    magnitudes = []
    for _, _, data in fields:
        selected = data[data.tval.isin(map_times)]
        magnitudes.extend(
            np.hypot(
                selected["induced_E_east"].to_numpy(),
                selected["induced_E_north"].to_numpy(),
            )
        )
    magnitudes = np.asarray(magnitudes, dtype=float)
    magnitudes = magnitudes[np.isfinite(magnitudes)]
    if magnitudes.size == 0:
        raise RuntimeError("The EMTF E-field values are all non-finite.")
    # Limit only the arrow display. The unmodified values remain in efield.png.
    arrow_scale = (
        float(quiver_scale)
        if quiver_scale is not None
        else max(3.0 * float(np.percentile(magnitudes, 98)), 1.0)
    )

    nrows, ncols = (2, 2) if len(map_times) <= 4 else (3, 2)
    figure_size = (8.0, 4.4) if nrows == 2 else (8.0, 7.2)
    map_figure = Map(
        map_times[0],
        nrows=nrows,
        ncols=ncols,
        region=region,
        figsize=figure_size,
    )
    for index, timestamp in enumerate(map_times):
        map_figure.ax.text(
            0.1,
            0.9,
            f"({chr(97 + index)}) / {timestamp.strftime('%H:%M UT')}",
            ha="left",
            va="center",
            fontsize=12,
            transform=map_figure.ax.transAxes,
        )
        lon, lat, ex, ey = [], [], [], []
        for site, station_name, data in fields:
            rows = data[data.tval == timestamp]
            if rows.empty:
                continue
            row = rows.iloc[0]
            if np.hypot(row.induced_E_east, row.induced_E_north) < 1e-6:
                continue
            lon.append((float(site.longitude) + 180) % 360 - 180)
            lat.append(float(site.latitude))
            ex.append(float(row.induced_E_east) / 500.0)
            ey.append(float(row.induced_E_north) / 500.0)
        if lon:
            quiver = map_figure.ax.quiver(
                np.asarray(lon),
                np.asarray(lat),
                np.asarray(ex),
                np.asarray(ey),
                transform=ccrs.PlateCarree(),
                color="#d62728",
                pivot="tail",
                scale=8,
                scale_units="width",
                width=0.002,
                headwidth=4,
                headlength=5,
                headaxislength=4,
                zorder=4,
            )
            if index == 0:
                map_figure.ax.quiverkey(
                    quiver,
                    0.50,
                    0.035,
                    1,
                    r"$\vec{E}$ [500 mV/km]",
                    labelpos="N",
                    coordinates="figure",
                    color="k",
                    fontproperties={"size": 8},
                )
        if index < len(map_times) - 1:
            map_figure.ax = map_figure.add_axes()

    map_figure.fig.subplots_adjust(
        left=0.015,
        right=0.985,
        bottom=0.14,
        top=0.985,
        wspace=0.02,
        hspace=0.0,
    )
    map_figure.save(output)
    map_figure.close()
    return output


def plot_case_study(goes, smag, event, dates, stations, case_station, output,
                    mt_data_folder=EMTF_DATA_FOLDER, dT=60,
                    map_times=None, case_mt_site=None, case_xlim=None,
                    case_top_stations=None, case_efield_ylim=None,
                    case_efield_quiver_limit=5000.0,
                    case_efield_panel_c_threshold_mV_per_km=500.0,
                    case_efield_exclude_sites=None,
                    case_baseline_hours=2.0,
                    efield_contour_bin_size_deg=1.0,
                    efield_contour_threshold_mV_per_km=0.25,
                    efield_contour_smoothing_sigma=1.5,
                    efield_bin_statistic="percentile",
                    efield_bin_percentile=75.0,
                    efield_min_stations_per_bin=3,
                    efield_min_neighbor_bins=2):
    """Create the case-study time series plus two E-field map panels."""
    station_colors = {
        "FRD": "#4c78a8",
        "FRN": "#f58518",
        "TUC": "#54a24b",
    }
    data_by_station = {
        station: smag.sm_data[smag.sm_data.iaga == station].copy()
        for station in stations
    }
    data_by_station = {
        station: data for station, data in data_by_station.items() if not data.empty
    }
    if case_station not in data_by_station:
        raise ValueError(f"Case-study station '{case_station}' has no SuperMAG data.")

    # The reference case figure compares these three stations.  Keep the
    # larger station list for stack-mag, but avoid making this publication
    # panel visually crowded when the configuration contains extra stations.
    if case_top_stations:
        top_stations = [
            name for name in case_top_stations if name in data_by_station
        ]
    else:
        reference_stations = [name for name in ("FRD", "FRN", "TUC")
                              if name in data_by_station]
        top_stations = reference_stations or list(data_by_station)

    if case_xlim is None:
        case_xlim = (
            event.replace(hour=17, minute=0, second=0, microsecond=0),
            event.replace(hour=20, minute=0, second=0, microsecond=0),
        )
    reference_time = case_xlim[0]
    baseline_start = reference_time - dt.timedelta(hours=float(case_baseline_hours))

    case_data = data_by_station[case_station]
    magnetic_north, magnetic_east = _case_baseline_corrected_components(
        case_data, baseline_start, reference_time
    )
            
    site = None
    if case_mt_site:
        requested = case_mt_site.upper()
        candidates = _load_emtf_sites(mt_data_folder)
        site = next(
            (candidate for candidate in candidates
             if str(getattr(candidate, "name", "")).upper() == requested),
            None,
        )
        if site is None:
            raise RuntimeError(
                f"Requested EMTF site '{case_mt_site}' was not found under "
                f"{mt_data_folder}."
            )
        e_north, e_east = site.convolve_fft(
            magnetic_north, magnetic_east, dt=dT
        )
    else:
        e_north, e_east, site, _ = compute_Egeo_from_Bgeo(
            magnetic_north,
            magnetic_east,
            glat=float(case_data.iloc[0].glat),
            glon=float(case_data.iloc[0].glon),
            dT=dT,
            MT_data_folder=mt_data_folder,
        )
    e_north = _reference_series_at_time(e_north, case_data["tval"], reference_time)
    e_east = _reference_series_at_time(e_east, case_data["tval"], reference_time)
    if site is None:
        raise RuntimeError(f"No EMTF response found for case station '{case_station}'.")

    # # Calculate a second E-field response at VAQ55 using the same magnetic
    # # station and magnetic perturbation time series.
    # vaq55_site = None
    # e_north_55 = e_east_55 = None
    # site_name = str(getattr(site, "name", "")).upper()
    # if site_name != "VAQ55":
    #     available_sites = candidates if case_mt_site else _load_emtf_sites(mt_data_folder)
    #     vaq55_site = next(
    #         (candidate for candidate in available_sites
    #          if str(getattr(candidate, "name", "")).upper() == "VAQ55"),
    #         None,
    #     )
    #     if vaq55_site is None:
    #         raise RuntimeError(f"EMTF site 'VAQ55' was not found under {mt_data_folder}.")
    #     e_north_55, e_east_55 = vaq55_site.convolve_fft(
    #         magnetic_north, magnetic_east, dt=dT
    #     )
    #     e_north_55 = _reference_series_at_time(
    #         e_north_55, case_data["tval"], reference_time
    #     )
    #     e_east_55 = _reference_series_at_time(
    #         e_east_55, case_data["tval"], reference_time
    #     )

    apply_style = {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.linewidth": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 10,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "savefig.dpi": 300,
        "figure.dpi": 150,
    }
    with plt.rc_context(apply_style):
        fig = plt.figure(figsize=(7.4*0.8, 7.4*0.8), dpi=300)
        grid = fig.add_gridspec(3, 2, height_ratios=(1.0, 1.2, 0.95),
                                hspace=0.28, wspace=0.06)
        ax_top = fig.add_subplot(grid[0, 0])
        fig.delaxes(ax_top)
        ax_top = fig.add_subplot(grid[0, :])
        ax_bottom = fig.add_subplot(grid[1, 0], sharex=ax_top)
        fig.delaxes(ax_bottom)
        ax_bottom = fig.add_subplot(grid[1, :], sharex=ax_top)

        goes_data = goes.dfs["goes"]
        ax_top.semilogy(
            goes_data.time,
            goes_data.xrsb,
            color="#c51b7d",
            lw=1.2,
            label="GOES 0.1-0.8 nm",
        )
        ax_top.set_ylabel("GOES 0.1-0.8 nm\n(W m$^{-2}$)", color="#c51b7d")
        ax_top.tick_params(axis="y", colors="#c51b7d")
        ax_top.set_ylim(1e-7, 1e-2)
        ax_top.grid(axis="y", which="major", linestyle=":", color="0.55", linewidth=0.45)
        ax_top.axvline(event, color="0.15", lw=0.7)
        ax_top.set_xlim(case_xlim)
        mag_ax = ax_top.twinx()
        for station in top_stations:
            data = data_by_station[station]
            north, east = _case_baseline_corrected_components(
                data, baseline_start, reference_time
            )
            field = np.sqrt(north ** 2 + east ** 2)
            field = _reference_series_at_time(field, data["tval"], reference_time)
            field *= -1.0
            mag_ax.plot(
                data.tval,
                field,
                lw=0.9,
                ls="-",
                color=station_colors.get(station, "#333333"),
                label=station,
            )
        mag_ax.set_ylabel(r"$\Delta\vec{B}$ (nT)")
        mag_ax.set_ylim(-150, 50)
        mag_ax.axvline(event, color="0.15", lw=0.7)
        mag_ax.legend(loc="lower right", ncol=3, frameon=False,
                      handlelength=1.6, columnspacing=1.0)
        ax_top.text(0.015, 0.94, "(a)", transform=ax_top.transAxes, fontweight="bold")

        ax_bottom.grid(axis="y", linestyle=":", color="0.55", linewidth=0.45)
        e_magnitude = np.hypot(e_east, e_north)
        primary_site_name = str(getattr(site, "name", None) or case_mt_site or "MT").upper()
        primary_plot_magnitude = np.where(e_magnitude > 0.0, e_magnitude, np.nan)
        primary_line, = ax_bottom.plot(
            case_data.tval,
            primary_plot_magnitude,
            color="#1f77b4",
            lw=1.0,
            ls="-",
            label=rf"$|\vec{{\mathcal{{E}}}}|$ ({primary_site_name})",
        )
        magnitude_series = [e_magnitude]
        legend_handles = [primary_line]
        # if e_north_55 is not None and e_east_55 is not None:
        #     e_magnitude_55 = np.hypot(e_east_55, e_north_55)
        #     secondary_plot_magnitude = np.where(
        #         e_magnitude_55 > 0.0, e_magnitude_55, np.nan
        #     )
        #     secondary_line, = ax_bottom.plot(
        #         case_data.tval,
        #         secondary_plot_magnitude,
        #         color="#d62728",
        #         lw=1.0,
        #         ls="-",
        #         label=r"$|\vec{\mathcal{E}}|$ (VAQ55)",
        #     )
        #     magnitude_series.append(e_magnitude_55)
        #     legend_handles.append(secondary_line)
        ax_bottom.set_ylabel(r"$|\vec{\mathcal{E}}|$ (mV/km)")
        finite_magnitude = np.concatenate([
            np.asarray(series)[np.isfinite(series)] for series in magnitude_series
        ])
        # ax_bottom.set_yscale("log")
        if case_efield_ylim is None:
            ax_bottom.set_ylim(1e-1, 2000)
        else:
            lower, upper = case_efield_ylim
            ax_bottom.set_ylim(max(float(lower), 1e-1), upper)
        # ax_bottom.set_yscale("log")
        # if positive_magnitude.size:
        #     magnitude_min = float(np.min(positive_magnitude))
        #     magnitude_max = float(np.max(positive_magnitude))
        # else:
        #     magnitude_min, magnitude_max = 1.0, 10.0
        # if case_efield_ylim is None:
        #     # ax_bottom.set_ylim(1e-1, 2e3)
        #     # ax_bottom.set_ylim(magnitude_min * 0.8, magnitude_max * 1.2)
        # else:
        #     lower, upper = case_efield_ylim
        #     lower = lower if lower > 0.0 else magnitude_min * 0.8
        #     upper = upper if upper > lower else magnitude_max * 1.2
        #     ax_bottom.set_ylim(1e-1, 2e3)
        ax_bottom.axvline(event, color="0.15", lw=0.7)
        ax_bottom.set_xlabel(f"UT on {event:%Y-%m-%d}")
        ax_bottom.legend(
            handles=legend_handles,
            loc="upper right",
            frameon=False,
            handlelength=1.6,
        )
        ax_bottom.text(0.015, 0.94, "(b)", transform=ax_bottom.transAxes, fontweight="bold")
        mt_site_label = str(
            getattr(site, "name", None) or case_mt_site or "MT"
        ).upper()
        equation = (
            rf"$\vec{{\mathcal{{E}}}}(\mathrm{{MT}}) = "
            rf"\mathbf{{Z}}(\omega,\mathrm{{MT}})"
            rf"\otimes \Delta\vec{{B}}(\mathrm{{{case_station.upper()}}})$"
        )
        ax_bottom.text(
            1.01,
            1.08,
            equation,
            transform=ax_bottom.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            clip_on=False,
        )

        # if hasattr(site, "periods") and hasattr(site, "resistivity"):
        #     inset = ax_bottom.inset_axes([0.70, 0.49, 0.23, 0.40])
        #     periods = np.asarray(site.periods)
        #     resistivity = np.asarray(site.resistivity)
        #     phase = np.asarray(site.phase)
        #     order = np.argsort(periods)
        #     periods = periods[order]
        #     resistivity = np.atleast_2d(resistivity)[:, order]
        #     phase = np.atleast_2d(phase)[:, order]
        #     period_mask = (periods >= 1.0) & (periods <= 1.0e4)
        #     periods = periods[period_mask]
        #     resistivity = resistivity[:, period_mask]
        #     phase = phase[:, period_mask]
        #     for row in np.atleast_2d(resistivity)[:2]:
        #         inset.loglog(periods, row, color="#3182bd", lw=0,
        #                      ls="None", marker=".", ms=4.0)
        #     inset.set_xlabel("Period (s)", fontsize=7)
        #     inset.set_ylabel(r"$\rho_a$ ($\Omega$m)", fontsize=7, color="#3182bd")
        #     inset.set_xlim(1.0, 1.0e4)
        #     inset.set_ylim(1e-1, 1e2)
        #     inset.tick_params(labelsize=6)
        #     inset.tick_params(which="both", direction="in")
        #     phase_ax = inset.twinx()
        #     for row in np.atleast_2d(phase)[:2]:
        #         phase_ax.semilogx(periods, row, color="#8c6bb1", lw=0,
        #                           ls="None", marker=".", ms=4.0)
        #     phase_ax.set_ylabel("Phase", fontsize=7, color="#8c6bb1")
        #     phase_ax.set_ylim(-180, 180)
        #     phase_ax.tick_params(labelsize=6)
        #     phase_ax.tick_params(which="both", direction="in")

        # Add two compact CONUS map panels below the time-series panels.
        if map_times is None or len(map_times) == 0:
            map_times = [event - dt.timedelta(minutes=50), event + dt.timedelta(minutes=20)]
        map_times = list(map_times)[:2]
        map_fields = _emtf_site_fields(
            smag, mt_data_folder, dT=dT, bounds=CONUS_BOUNDS
        )
        map_axes = []
        panel_c_excluded_sites = set()
        panel_c_exclusion_rows = []
        map_projection = ccrs.LambertConformal(
            central_longitude=-96.0,
            central_latitude=39.0,
            standard_parallels=(29.0, 45.0),
        )
        for index, timestamp in enumerate(map_times):
            map_ax = fig.add_subplot(grid[2, index], projection=map_projection)
            map_axes.append(map_ax)
            map_ax.set_extent(CONUS_BOUNDS, crs=ccrs.PlateCarree())
            map_ax.add_feature(
                Nightshade(timestamp, alpha=0.18, color="#4d4d4d"),
                zorder=0.5,
            )
            map_ax.coastlines(resolution="50m", color="0.25", linewidth=0.5)
            map_ax.add_feature(cfeature.BORDERS, edgecolor="0.25", linewidth=0.35)
            gridliner = map_ax.gridlines(
                draw_labels=False, linewidth=0.25, color="0.75",
                alpha=0.6, linestyle=":"
            )
            map_ax.text(
                0.05, 0.1, timestamp.strftime("%H:%M UT"),
                transform=map_ax.transAxes, ha="center", va="bottom",
                fontsize=9, rotation=90
            )
            map_ax.text(0.02, 0.94, f"({chr(99 + index)})",
                        transform=map_ax.transAxes, fontsize=10,
                        fontweight="bold", va="top")

            lon, lat, ex, ey = [], [], [], []
            contour_lon, contour_lat, contour_magnitude = [], [], []
            reference = 500.0  # mV/km
            excluded_sites = {
                str(name).upper() for name in (case_efield_exclude_sites or [])
            }
            for site_item, station_name, data in map_fields:
                if data.empty:
                    continue
                site_name = str(getattr(site_item, "name", "")).upper()
                if site_name in excluded_sites:
                    continue
                delta = (data.tval - timestamp).abs()
                row = data.loc[delta.idxmin()]
                if abs((row.tval - timestamp).total_seconds()) > 90:
                    continue
                east_value = float(row.induced_E_east)
                north_value = float(row.induced_E_north)
                if not (np.isfinite(east_value) and np.isfinite(north_value)):
                    continue
                magnitude = np.hypot(east_value, north_value)
                site_lon = (float(site_item.longitude) + 180.0) % 360.0 - 180.0
                site_lat = float(site_item.latitude)
                if index == 0 and (
                    case_efield_panel_c_threshold_mV_per_km is not None
                    and magnitude > float(case_efield_panel_c_threshold_mV_per_km)
                ):
                    panel_c_excluded_sites.add(site_name)
                    try:
                        rating = float(getattr(site_item, "rating", np.nan))
                    except (TypeError, ValueError):
                        rating = np.nan
                    panel_c_exclusion_rows.append({
                        "emtf_site": site_name,
                        "nearest_supermag": station_name,
                        "latitude": site_lat,
                        "longitude": site_lon,
                        "mt_rating": rating,
                        "map_time_utc": timestamp.strftime("%Y-%m-%dT%H:%M:%S"),
                        "E_north_mV_per_km": north_value,
                        "E_east_mV_per_km": east_value,
                        "E_magnitude_mV_per_km": magnitude,
                        "threshold_mV_per_km": float(
                            case_efield_panel_c_threshold_mV_per_km
                        ),
                    })
                    continue
                if index == 1 and site_name in panel_c_excluded_sites:
                    continue
                contour_lon.append(site_lon)
                contour_lat.append(site_lat)
                contour_magnitude.append(magnitude)
                # Keep quivers restricted to stations with the requested
                # quality rating, while contours use every valid response
                # site so regional interpolation is not undersampled.
                if site_item.rating < 3:
                    continue
                if (case_efield_quiver_limit is not None and
                        magnitude > float(case_efield_quiver_limit)):
                    continue
                # Compare in physical units before normalizing by the
                # 500 mV/km reference used by the quiver scale.
                lon.append(site_lon)
                lat.append(site_lat)
                ex.append(east_value / reference)
                ey.append(north_value / reference)

            if index == 1:
                plot_efield_threshold_contour(
                    map_ax,
                    contour_lon,
                    contour_lat,
                    contour_magnitude,
                    efield_contour_threshold_mV_per_km,
                    bounds=CONUS_BOUNDS,
                    bin_size_deg=efield_contour_bin_size_deg,
                    smoothing_sigma=efield_contour_smoothing_sigma,
                    statistic=efield_bin_statistic,
                    percentile=efield_bin_percentile,
                    minimum_stations=efield_min_stations_per_bin,
                    minimum_neighbor_bins=efield_min_neighbor_bins,
                )
                
            quiver = map_ax.quiver(
                np.asarray(lon), np.asarray(lat),
                np.asarray(ex), np.asarray(ey),
                transform=ccrs.PlateCarree(), color="#d62728",
                pivot="tail", scale=8, scale_units="width", width=0.002,
                headwidth=4, headlength=5, headaxislength=4, zorder=4,
            )
            
            if (quiver is not None) and (index==1):
                qk = map_ax.quiverkey(
                    quiver, 1.05, 0.9, 1.0,
                    label=r"$\vec{\mathcal{E}}=$ 500 mV/km", labelpos="E",
                    coordinates="axes", color="#d62728",
                    fontproperties={"size": 7},
                    angle=90
                )
                qk.text.set_rotation(90)

        exclusion_file = os.path.splitext(output)[0] + "-panel-c-excluded-sites.csv"
        with open(exclusion_file, "w", newline="", encoding="utf-8") as stream:
            fieldnames = [
                "emtf_site", "nearest_supermag", "latitude", "longitude",
                "mt_rating", "map_time_utc", "E_north_mV_per_km",
                "E_east_mV_per_km", "E_magnitude_mV_per_km",
                "threshold_mV_per_km",
            ]
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(panel_c_exclusion_rows)

        for ax in (ax_top, ax_bottom):
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
            ax.tick_params(which="both", top=True, right=True)
        plt.setp(ax_top.get_xticklabels(), visible=False)
        fig.subplots_adjust(left=0.10, right=0.92, bottom=0.08, top=0.96,
                            hspace=0.28, wspace=0.06)
        fig.savefig(output, dpi=300)
        plt.close(fig)
    return output
