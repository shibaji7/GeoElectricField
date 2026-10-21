#!/usr/bin/env python3
"""Generate the standard GeoElectricField figures for one event.

Examples
--------
    python event_study.py --event 2005-09-07T17:40 \
        --start 2005-09-07T17:15 --end 2005-09-07T19:00 \
        --stations OTT FRD PIN

If ``--start`` and ``--end`` are omitted, the study window is 20 minutes
before through 30 minutes after the event.
"""

import argparse
import datetime as dt
import glob
import os
import sys

import numpy as np
import pandas as pd
import toml


PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_GOES_DATA_DIR = "/home/chakras4/sunpy/data"
DEFAULT_OMNI_DATA_DIR = "/home/chakras4/OMNI/"
DEFAULT_EMTF_DATA_DIR = os.path.join(PROJECT_ROOT, ".data", "EMTF")

# The legacy project assumes a different home directory.  Keep SunPy's
# configuration writable before importing any module that imports SunPy.
os.environ.setdefault("SUNPY_CONFIGDIR", "/tmp/gef-sunpy")
os.environ["OMNIDATA_PATH"] = DEFAULT_OMNI_DATA_DIR
sys.path.extend(
    [
        os.path.join(PROJECT_ROOT, "py"),
        os.path.join(PROJECT_ROOT, "py", "fetch"),
        os.path.join(PROJECT_ROOT, "py", "process"),
        os.path.join(PROJECT_ROOT, "py", "plots"),
    ]
)

from get_mag_data import SuperMAG
from efield import (
    CONUS_BOUNDS,
    _emtf_site_fields,
    plot_case_study,
    plot_efield,
    plot_efield_map,
)
from maps import Map
# Import this before simulation.py.  pyomnidata captures OMNIDATA_PATH at
# import time, while simulation.py contains a legacy path assignment.
from pyomnidata import GetOMNI as _GetOMNI
from simulation import getOMNI, setup, stack_mag_mod_plots, stack_plots


def parse_datetime(value):
    """Parse an ISO-8601 datetime accepted by the command-line interface."""
    return dt.datetime.fromisoformat(value.replace("Z", "+00:00")).replace(
        tzinfo=None
    )


def event_base(event):
    return os.path.join(
        PROJECT_ROOT,
        "tmp",
        event.strftime("%Y-%m-%d-%H-%M"),
    )


def cache_covers_window(filename, dates):
    """Return True when a cached SuperMAG CSV spans the requested window."""
    if not os.path.exists(filename):
        return False
    try:
        times = pd.read_csv(filename, usecols=["tval"], parse_dates=["tval"])["tval"]
    except (OSError, ValueError, pd.errors.EmptyDataError):
        return False
    times = times.dropna()
    return not times.empty and times.min() <= dates[0] and times.max() >= dates[1]


def ensure_dbdt(smag):
    """Add the magnetic-field derivative expected by stack_plots."""
    if "dBdt" in smag.sm_data.columns:
        return

    data = smag.sm_data.copy()
    data["B"] = np.sqrt(
        data["N_geo"] ** 2 + data["E_geo"] ** 2 + data["Z_geo"] ** 2
    )
    data["dBdt"] = 0.0
    for _, indices in data.groupby("iaga").groups.items():
        station = data.loc[indices].sort_values("tval")
        delta_t = station["tval"].diff().dt.total_seconds()
        delta_t = delta_t.replace(0, np.nan).fillna(1.0)
        derivative = station["B"].diff().fillna(0.0) / delta_t
        data.loc[station.index, "dBdt"] = derivative
    smag.sm_data = data


def add_baseline_subtracted_horizontal(smag, dates, baseline_minutes=20):
    """Add station-wise baseline-subtracted horizontal magnetic vectors.

    The E/N baselines are the station-wise medians during the first
    ``baseline_minutes`` of the study window.  The map uses the component
    perturbations ``dE`` and ``dN`` directly, preserving the original map
    vector style while removing the station baseline.
    """
    data = smag.sm_data.copy()
    data["H_geo"] = np.hypot(data["N_geo"], data["E_geo"])
    baseline_end = dates[0] + dt.timedelta(minutes=baseline_minutes)
    baseline_rows = data[data["tval"].between(dates[0], baseline_end)]
    baselines = baseline_rows.groupby("iaga")[["E_geo", "N_geo"]].median()
    data["E_baseline"] = data["iaga"].map(baselines["E_geo"])
    data["N_baseline"] = data["iaga"].map(baselines["N_geo"])
    data["dE_geo"] = data["E_geo"] - data["E_baseline"]
    data["dN_geo"] = data["N_geo"] - data["N_baseline"]
    data["dH"] = np.hypot(data["dE_geo"], data["dN_geo"])
    smag.sm_data = data


def add_baseline_subtracted_H(smag, dates, baseline_minutes=20):
    """Add baseline-subtracted magnetic components and ``H - H_baseline``."""
    data = smag.sm_data.copy()
    data["H_geo"] = np.hypot(data["N_geo"], data["E_geo"])
    baseline_end = dates[0] + dt.timedelta(minutes=baseline_minutes)
    baseline_rows = data[data["tval"].between(dates[0], baseline_end)]
    baselines = baseline_rows.groupby("iaga")[["N_geo", "E_geo", "Z_geo"]].median()
    for component in ("N", "E", "Z"):
        data[f"{component}_baseline"] = data["iaga"].map(
            baselines[f"{component}_geo"]
        )
        data[f"dB_{component}_geo"] = (
            data[f"{component}_geo"] - data[f"{component}_baseline"]
        )
    data["dB_H_geo"] = np.hypot(data["dB_N_geo"], data["dB_E_geo"])
    data["H_baseline"] = np.hypot(
        data["N_baseline"], data["E_baseline"]
    )
    data["dH_geo"] = data["H_geo"] - data["H_baseline"]
    smag.sm_data = data


def export_station_db_data(smag, stations, output_root, event, dates):
    """Export selected station magnetic data to the shared publication folder."""
    if not stations:
        return []
    share_dir = os.path.join(output_root, "share")
    os.makedirs(share_dir, exist_ok=True)
    outputs = []
    columns = [
        "tval",
        "N_geo", "E_geo", "Z_geo",
        "dB_N_geo", "dB_E_geo", "dB_Z_geo", "dB_H_geo",
        "H_geo", "H_baseline", "dH_geo", "dBdt",
    ]
    for station in stations:
        data = smag.sm_data[smag.sm_data["iaga"] == station].copy()
        data = data[data["tval"].between(dates[0], dates[1])]
        if data.empty:
            continue
        available = [column for column in columns if column in data.columns]
        filename = os.path.join(share_dir, f"{station}_dB.csv")
        with open(filename, "w", encoding="utf-8", newline="") as stream:
            stream.write("# GeoElectricField SuperMAG magnetic data\n")
            stream.write(f"# Station: {station}\n")
            stream.write(f"# Event: {event:%Y-%m-%dT%H:%M} UT\n")
            stream.write(
                f"# Time range: {dates[0]:%Y-%m-%d %H:%M} to "
                f"{dates[1]:%H:%M} UT\n"
            )
            stream.write("# Magnetic field units: nT; dBdt units: nT/s\n")
            stream.write(
                "# Baseline: station-wise median during the first 20 minutes "
                "of the study window\n"
            )
            stream.write("# dB_N_geo = N_geo - N_baseline\n")
            stream.write("# dB_E_geo = E_geo - E_baseline\n")
            stream.write("# dB_Z_geo = Z_geo - Z_baseline\n")
            stream.write("# dB_H_geo = sqrt(dB_N_geo^2 + dB_E_geo^2)\n")
            stream.write("# dH_geo = H_geo - H_baseline\n")
            data.sort_values("tval")[available].to_csv(stream, index=False)
        outputs.append(filename)
    return outputs


def export_emtf_efield_data(smag, mt_data_folder, output_root, event, dates):
    """Export E-field time series for every CONUS EMTF site."""
    share_dir = os.path.join(output_root, "share")
    os.makedirs(share_dir, exist_ok=True)
    fields = _emtf_site_fields(
        smag,
        mt_data_folder,
        dT=60,
        bounds=CONUS_BOUNDS,
    )
    outputs = []
    for site, magnetometer_station, data in fields:
        site_name = str(getattr(site, "name", "MT")).upper()
        selected = data[data["tval"].between(dates[0], dates[1])].copy()
        if selected.empty:
            continue
        selected["induced_E_magnitude"] = np.hypot(
            selected["induced_E_east"], selected["induced_E_north"]
        )
        filename = os.path.join(share_dir, f"{site_name}_Efield.csv")
        with open(filename, "w", encoding="utf-8", newline="") as stream:
            stream.write("# GeoElectricField induced electric-field data\n")
            stream.write(f"# EMTF site: {site_name}\n")
            stream.write(f"# Nearest SuperMAG station: {magnetometer_station}\n")
            stream.write(f"# Event: {event:%Y-%m-%dT%H:%M} UT\n")
            stream.write(
                f"# Time range: {dates[0]:%Y-%m-%d %H:%M} to "
                f"{dates[1]:%H:%M} UT\n"
            )
            stream.write("# Electric-field units: mV/km\n")
            stream.write(
                "# E_magnitude = sqrt(induced_E_east^2 + "
                "induced_E_north^2)\n"
            )
            selected[
                [
                    "tval",
                    "induced_E_east",
                    "induced_E_north",
                    "induced_E_magnitude",
                ]
            ].to_csv(stream, index=False)
        outputs.append(filename)
    return outputs


class CachedGOES:
    """Load cached GOES NetCDF4 files without SunPy's file parser."""

    def __init__(self, dates, data_dir=DEFAULT_GOES_DATA_DIR):
        self.data_dir = data_dir
        self.dfs = {"goes": self._load(dates)}

    def _files_for_dates(self, dates):
        files = []
        day = dates[0].date()
        while day <= dates[1].date():
            files.extend(
                glob.glob(
                    os.path.join(
                        self.data_dir,
                        f"sci_gxrs*{day.strftime('%Y%m%d')}*.nc",
                    )
                )
            )
            day += dt.timedelta(days=1)
        return sorted(set(files))

    def _download(self, dates):
        from sunpy.net import Fido
        from sunpy.net import attrs as a

        os.makedirs(self.data_dir, exist_ok=True)
        result = Fido.search(
            a.Time(
                dates[0].strftime("%Y-%m-%d %H:%M"),
                dates[1].strftime("%Y-%m-%d %H:%M"),
            ),
            a.Instrument("XRS"),
        )
        if len(result) == 0:
            raise FileNotFoundError("SunPy found no GOES XRS data for this window.")

        Fido.fetch(
            result,
            path=os.path.join(self.data_dir, "{file}"),
            progress=True,
        )

    def _load(self, dates):
        files = self._files_for_dates(dates)
        if not files:
            self._download(dates)
            files = self._files_for_dates(dates)

        if not files:
            raise FileNotFoundError(
                f"No GOES NetCDF files were found under {self.data_dir} after download."
            )

        from netCDF4 import Dataset, num2date

        frames = []
        for filename in sorted(set(files)):
            with Dataset(filename) as dataset:
                variables = dataset.variables
                if "xrsa_flux" in variables and "xrsb_flux" in variables:
                    flux_a, flux_b = "xrsa_flux", "xrsb_flux"
                elif "a_flux" in variables and "b_flux" in variables:
                    flux_a, flux_b = "a_flux", "b_flux"
                elif "A_AVG" in variables and "B_AVG" in variables:
                    flux_a, flux_b = "A_AVG", "B_AVG"
                elif "A_FLUX" in variables and "B_FLUX" in variables:
                    flux_a, flux_b = "A_FLUX", "B_FLUX"
                elif "xs" in variables and "xl" in variables:
                    flux_a, flux_b = "xs", "xl"
                else:
                    continue

                time_name = "time" if "time" in variables else "time_tag"
                time = variables[time_name]
                frames.append(
                    pd.DataFrame(
                        {
                            "time": pd.to_datetime(
                                num2date(
                                    time[:],
                                    units=time.units,
                                    calendar=getattr(time, "calendar", "standard"),
                                    only_use_cftime_datetimes=False,
                                )
                            ),
                            "xrsa": variables[flux_a][:],
                            "xrsb": variables[flux_b][:],
                        }
                    )
                )

        data = pd.concat(frames, ignore_index=True)
        # Fido can return overlapping files in filename order. Sort before
        # slicing so traces cannot jump backward at a file boundary.
        data = data.sort_values("time").drop_duplicates("time")
        return data[(data.time >= dates[0]) & (data.time <= dates[1])].reset_index(
            drop=True
        )


def make_map(smag, dates, output_file, map_times=None, map_region="north_america"):
    """Create the multi-panel magnetic-vector map figure."""
    if map_times:
        map_dates = map_times
    else:
        map_dates = [
            dates[0] + dt.timedelta(minutes=minutes)
            for minutes in range(
                int((dates[1] - dates[0]).total_seconds() / 60)
            )
        ][:6]

    if not map_dates:
        raise ValueError("The study window must be longer than one minute.")
    if len(map_dates) > 6:
        raise ValueError("map.png supports at most six map times.")

    map_figure = Map(
        map_dates[0], nrows=2, ncols=3, region=map_region
    )
    map_figure.ax.text(
        0.1,
        1.05,
        "Coords: Geo",
        ha="left",
        fontweight="bold",
        fontsize=12,
        transform=map_figure.ax.transAxes,
    )

    for index, timestamp in enumerate(map_dates[:6]):
        map_figure.ax.text(
            0.1,
            0.9,
            f"({chr(97 + index)}) {timestamp.strftime('%H:%M UT')}",
            ha="left",
            va="center",
            fontsize=12,
            transform=map_figure.ax.transAxes,
        )

        for idx, station in enumerate(smag.stations):
            rows = smag.sm_data[
                (smag.sm_data.tval == timestamp)
                & (smag.sm_data.iaga == station)
            ]
            if len(rows) == 0:
                continue
            instrument = {
                "lon": (rows.glon.iloc[0] + 180) % 360 - 180,
                "lat": rows.glat.iloc[0],
                "code": rows.iaga.iloc[0],
                "E_geo": rows.E_geo.iloc[0],
                "N_geo": rows.N_geo.iloc[0],
                "Z_geo": rows.Z_geo.iloc[0],
            }
            map_figure.plot_instrument(
                instrument,
                show_label=False,
                quiver_key=index == 0,
                quiver_key_value=50,
                quiver_key_label=r"$\Delta \vec{B}$ [50 nT]",
            )

        if index < min(5, len(map_dates) - 1):
            map_figure.ax = map_figure.add_axes()

    map_figure.fig.subplots_adjust(wspace=0.05, hspace=0.05)
    map_figure.save(output_file)
    map_figure.close()


def generate_event_figures(
    event,
    dates,
    mag_stations=None,
    fetch_stations=None,
    map_times=None,
    efield_map_times=None,
    map_region="north_america",
    output_root=None,
    supermag_user="shibaji7",
    mag_baseline=None,
    goes_data_dir=DEFAULT_GOES_DATA_DIR,
    omni_data_dir=DEFAULT_OMNI_DATA_DIR,
    emtf_data_dir=DEFAULT_EMTF_DATA_DIR,
    efield_quiver_scale=None,
    efield_quiver_limit=None,
    stack_mag_ylim=None,
    stack_ae_ylim=None,
    case_xlim=None,
    case_top_stations=None,
    case_efield_ylim=None,
    case_efield_quiver_limit=5000.0,
    case_efield_exclude_sites=None,
    case_baseline_hours=2.0,
    efield_contour_bin_size_deg=1.0,
    efield_contour_threshold_mV_per_km=0.25,
    efield_contour_smoothing_sigma=1.5,
    efield_bin_statistic="percentile",
    efield_bin_percentile=75.0,
    efield_min_stations_per_bin=3,
    efield_min_neighbor_bins=2,
    case_station=None,
    case_map_times=None,
    case_mt_site=None,
    share_stations=None,
    share_efield_sites=False,
):
    """Fetch event data and generate stack, magnetic, and map figures."""
    base = event_base(event)
    output_root = output_root or os.path.join(PROJECT_ROOT, "pub_fig")
    if not os.path.isabs(output_root):
        output_root = os.path.join(PROJECT_ROOT, output_root)
    figures = os.path.join(output_root, event.strftime("%Y-%m-%d-%H-%M"))
    os.makedirs(figures, exist_ok=True)
    setup(base)
    goes = CachedGOES(dates, data_dir=goes_data_dir)
    case_reference_time = (
        case_xlim[0]
        if case_xlim is not None
        else event.replace(hour=17, minute=0, second=0, microsecond=0)
    )
    case_baseline_start = case_reference_time - dt.timedelta(
        hours=float(case_baseline_hours)
    )
    magnetic_dates = [
        min(dates[0], case_baseline_start),
        dates[1],
    ]
    raw_baseline = "none" if mag_baseline == "none" else mag_baseline
    magnetic_base = base + os.sep
    cache_name = (
        "supermag-raw.csv" if raw_baseline == "none" else "supermag.csv"
    )
    if not cache_covers_window(os.path.join(base, cache_name), magnetic_dates):
        magnetic_base = os.path.join(base, "case-study-extended") + os.sep
    # Fetch the full inventory by default.  The magnetic plot can use a
    # separately ordered subset without restricting the map data.
    smag = SuperMAG(
        magnetic_base,
        magnetic_dates,
        uid=supermag_user,
        stations=fetch_stations if fetch_stations else [],
        baseline=raw_baseline,
    )
    if mag_baseline == "none":
        # Also cache the API's yearly-baseline product separately.  The raw
        # product remains the source for map.png and H-H_baseline below.
        yearly_stations = sorted(
            smag.sm_data["iaga"].dropna().astype(str).unique().tolist()
        )
        smag_dB = SuperMAG(
            magnetic_base,
            magnetic_dates,
            uid=supermag_user,
            stations=yearly_stations,
            baseline=None,
        )
        # ensure_dbdt(smag_dB)
    ensure_dbdt(smag)
    
    if mag_baseline == "none":
        add_baseline_subtracted_H(smag, dates)
    share_stations = [str(station).upper() for station in (share_stations or [])]
    shared_data_files = export_station_db_data(
        smag, share_stations, output_root, event, dates
    )
    shared_efield_files = []
    if share_efield_sites:
        efield_dates = case_xlim or (dates[0], dates[1])
        shared_efield_files = export_emtf_efield_data(
            smag,
            emtf_data_dir,
            output_root,
            event,
            efield_dates,
        )
    # Some SuperMAG responses provide magnetic latitude but omit magnetic
    # longitude.  The legacy stack annotation requires both fields; use the
    # geographic coordinate as a plotting fallback when necessary.
    if "mlat" not in smag.sm_data.columns and "glat" in smag.sm_data.columns:
        smag.sm_data["mlat"] = smag.sm_data["glat"]
    if "mlon" not in smag.sm_data.columns and "glon" in smag.sm_data.columns:
        smag.sm_data["mlon"] = smag.sm_data["glon"]
    
    if "mlat" not in smag_dB.sm_data.columns and "glat" in smag_dB.sm_data.columns:
        smag_dB.sm_data["mlat"] = smag_dB.sm_data["glat"]
    if "mlon" not in smag_dB.sm_data.columns and "glon" in smag_dB.sm_data.columns:
        smag_dB.sm_data["mlon"] = smag_dB.sm_data["glon"]
    candidate_stations = (
        list(fetch_stations)
        if fetch_stations
        else sorted(smag.sm_data["iaga"].dropna().unique().tolist())
    )
    usable_stations = []
    for station in candidate_stations:
        rows = smag.sm_data[smag.sm_data.iaga == station]
        if len(rows) == 0:
            continue
        required = ["mlat", "mlon", "glat", "glon", "E_geo", "N_geo", "Z_geo"]
        if all(column in rows.columns and rows[column].notna().any() for column in required):
            usable_stations.append(station)
    if not usable_stations:
        raise RuntimeError("No SuperMAG stations returned usable data for this event.")

    # simulation.py contains a legacy hard-coded OMNIDATA_PATH.  Override it
    # here after importing that module and before GetOMNI is called.
    os.environ["OMNIDATA_PATH"] = omni_data_dir
    omni = getOMNI(
        event.year,
        [int(dates[0].strftime("%Y%m%d")), int(dates[1].strftime("%Y%m%d"))],
    )

    stack_file = os.path.join(figures, "stack.png")
    requested_mag = mag_stations or usable_stations
    ordered_mag = [station for station in requested_mag if station in usable_stations]
    if not ordered_mag:
        raise RuntimeError("None of the requested magnetic stations have usable data.")
    smag.stations = ordered_mag
    stack_plots(
        event,
        goes,
        smag,
        omni,
        dates,
        stack_file,
        stack_mag_ylim=stack_mag_ylim,
        stack_ae_ylim=stack_ae_ylim,
    )
    stack_mag_mod_file = stack_file.replace(".png", "-mag-mod.png")
    stack_mag_mod_plots(
        event,
        goes,
        smag,
        omni,
        dates,
        stack_file,
        stack_mag_ylim=stack_mag_ylim,
        stack_ae_ylim=stack_ae_ylim,
    )
    efield_file = os.path.join(figures, "efield.png")
    # plot_efield(
    #     smag, event, dates, ordered_mag, efield_file,
    #     mt_data_folder=emtf_data_dir,
    # )
    case_study_file = os.path.join(figures, "case-study.png")
    plot_case_study(
        goes,
        smag,
        event,
        dates,
        ordered_mag,
        case_station or ordered_mag[0],
        case_study_file,
        mt_data_folder=emtf_data_dir,
        map_times=case_map_times or efield_map_times or map_times,
        case_mt_site=case_mt_site,
        case_xlim=case_xlim,
        case_top_stations=case_top_stations,
        case_efield_ylim=case_efield_ylim,
        case_efield_quiver_limit=case_efield_quiver_limit,
        case_efield_exclude_sites=case_efield_exclude_sites,
        case_baseline_hours=case_baseline_hours,
        efield_contour_bin_size_deg=efield_contour_bin_size_deg,
        efield_contour_threshold_mV_per_km=efield_contour_threshold_mV_per_km,
        efield_contour_smoothing_sigma=efield_contour_smoothing_sigma,
        efield_bin_statistic=efield_bin_statistic,
        efield_bin_percentile=efield_bin_percentile,
        efield_min_stations_per_bin=efield_min_stations_per_bin,
        efield_min_neighbor_bins=efield_min_neighbor_bins,
    )
    print(smag)
    smag.stations = usable_stations
    efield_map_file = os.path.join(figures, "maps-efield.png")
    # plot_efield_map(
    #     smag,
    #     efield_map_times or map_times or [
    #         dates[0] + dt.timedelta(minutes=index)
    #         for index in range(min(6, int((dates[1] - dates[0]).total_seconds() / 60)))
    #     ],
    #     "conus",
    #     efield_map_file,
    #     mt_data_folder=emtf_data_dir,
    #     quiver_scale=efield_quiver_scale,
    #     quiver_limit=efield_quiver_limit,
    # )
    make_map(
        smag_dB,
        dates,
        os.path.join(figures, "map.png"),
        map_times=map_times,
        map_region=map_region,
    )

    return {
        "stack": stack_file,
        "stack_mag": stack_file.replace(".png", "-mag.png"),
        "stack_mag_mod": stack_mag_mod_file,
        "efield": efield_file,
        "case_study": case_study_file,
        "efield_map": efield_map_file,
        "map": os.path.join(figures, "map.png"),
        "shared_data": shared_data_files,
        "shared_efield_data": shared_efield_files,
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate stack.png, stack-mag.png, and map.png for an event."
    )
    parser.add_argument("--config", help="Event TOML configuration file")
    parser.add_argument("--event", help="Event time, e.g. 2005-09-07T17:40")
    parser.add_argument("--start", help="Study-window start in ISO format")
    parser.add_argument("--end", help="Study-window end in ISO format")
    parser.add_argument(
        "--stations",
        nargs="*",
        default=[],
        help="SuperMAG station codes separated by spaces; omit to use the inventory",
    )
    parser.add_argument(
        "--mag",
        help=(
            "Ordered station codes for stack-mag.png, comma-separated; "
            "for example RES,RAL,PIN,FRD,FRN,TUC"
        ),
    )
    parser.add_argument(
        "--map-times",
        help=(
            "Map panel times as comma-separated HH:MM values, for example "
            "17:20,17:40,18:00,18:20"
        ),
    )
    parser.add_argument(
        "--map-region",
        choices=["north_america", "europe", "north_polar"],
    )
    parser.add_argument(
        "--efield-map-times",
        help=(
            "E-field map panel times as comma-separated HH:MM values; "
            "overrides map.efield_map_times"
        ),
    )
    parser.add_argument(
        "--efield-quiver-scale",
        type=float,
        help="Matplotlib quiver scale for maps-efield.png.",
    )
    parser.add_argument(
        "--efield-quiver-limit",
        type=float,
        help="Maximum displayed E-field magnitude in mV/km.",
    )
    parser.add_argument(
        "--case-station",
        help="Station used for the lower panel of case-study.png.",
    )
    parser.add_argument(
        "--mag-baseline",
        choices=["none", "yearly"],
        help=(
            "SuperMAG baseline mode. 'none' downloads raw magnetic fields "
            "to supermag-raw.csv; omitted preserves the existing cache."
        ),
    )
    return parser


def main():
    args = build_parser().parse_args()
    config = {}
    config_path = args.config
    if config_path is None and args.event:
        event_id = parse_datetime(args.event).strftime("%Y-%m-%d-%H-%M")
        candidate = os.path.join(PROJECT_ROOT, "config", "event_study", f"{event_id}.toml")
        if os.path.exists(candidate):
            config_path = candidate
    if config_path:
        with open(config_path, "r") as config_file:
            config = toml.load(config_file)

    event_config = config.get("event", {})
    map_config = config.get("map", {})
    plot_config = config.get("plot", {})
    output_config = config.get("output", {})
    event = parse_datetime(args.event or event_config["time"])
    config_start = event_config.get("start")
    config_end = event_config.get("end")
    dates = [
        parse_datetime(args.start)
        if args.start
        else parse_datetime(config_start)
        if config_start
        else event - dt.timedelta(minutes=20),
        parse_datetime(args.end)
        if args.end
        else parse_datetime(config_end)
        if config_end
        else event + dt.timedelta(minutes=30),
    ]
    if dates[1] <= dates[0]:
        raise ValueError("--end must be later than --start.")

    config_mag = event_config.get("mag", [])
    mag_stations = (
        [station.strip() for station in args.mag.split(",") if station.strip()]
        if args.mag
        else config_mag
        if config_mag
        else args.stations
    )
    def parse_map_times(values, option_name):
        if not values:
            return None
        if isinstance(values, str):
            values = values.split(",")
        parsed_times = []
        for value in values:
            value = value.strip()
            try:
                parsed_time = dt.datetime.strptime(
                    value, "%H:%M" if ":" in value else "%H"
                ).time()
            except ValueError as exc:
                raise ValueError(
                    f"Invalid {option_name} value '{value}'. Use HH:MM or HH."
                ) from exc
            parsed_times.append(
                event.replace(
                    hour=parsed_time.hour,
                    minute=parsed_time.minute,
                    second=0,
                    microsecond=0,
                )
            )
        return parsed_times

    map_times = parse_map_times(
        args.map_times or event_config.get("map_times", []),
        "map time",
    )
    efield_map_times = parse_map_times(
        args.efield_map_times or map_config.get("efield_map_times", []),
        "E-field map time",
    )
    case_map_times = parse_map_times(
        map_config.get("case_map_times", []),
        "case-study map time",
    )
    case_xlim_values = parse_map_times(
        plot_config.get("case_xlim", []),
        "case-study x-limit",
    )
    if case_xlim_values and len(case_xlim_values) != 2:
        raise ValueError("plot.case_xlim must contain exactly two HH:MM values.")
    outputs = generate_event_figures(
        event,
        dates,
        mag_stations=mag_stations,
        fetch_stations=args.stations,
        map_times=map_times,
        efield_map_times=efield_map_times,
        case_map_times=case_map_times,
        case_mt_site=event_config.get("case_mt_site"),
        map_region=args.map_region or map_config.get("region", "north_america"),
        output_root=output_config.get("root"),
        supermag_user=event_config.get("supermag_user", "shibaji7"),
        mag_baseline=args.mag_baseline or event_config.get("mag_baseline"),
        stack_mag_ylim=plot_config.get("stack_mag_ylim"),
        stack_ae_ylim=plot_config.get("stack_ae_ylim"),
        case_xlim=tuple(case_xlim_values) if case_xlim_values else None,
        case_top_stations=plot_config.get("case_top_stations"),
        case_efield_ylim=plot_config.get("case_efield_ylim"),
        case_efield_quiver_limit=plot_config.get(
            "case_efield_quiver_limit", 5000.0
        ),
        case_efield_exclude_sites=plot_config.get(
            "case_efield_exclude_sites", []
        ),
        case_baseline_hours=plot_config.get("case_baseline_hours", 2.0),
        efield_contour_bin_size_deg=plot_config.get(
            "efield_contour_bin_size_deg", 1.0
        ),
        efield_contour_threshold_mV_per_km=plot_config.get(
            "efield_contour_threshold_mV_per_km", 0.25
        ),
        efield_contour_smoothing_sigma=plot_config.get(
            "efield_contour_smoothing_sigma", 1.5
        ),
        efield_bin_statistic=plot_config.get(
            "efield_bin_statistic", "percentile"
        ),
        efield_bin_percentile=plot_config.get(
            "efield_bin_percentile", 75.0
        ),
        efield_min_stations_per_bin=plot_config.get(
            "efield_min_stations_per_bin", 3
        ),
        efield_min_neighbor_bins=plot_config.get(
            "efield_min_neighbor_bins", 2
        ),
        share_stations=output_config.get("share_stations", []),
        share_efield_sites=output_config.get("share_efield_sites", False),
        goes_data_dir=event_config.get("goes_data_dir", DEFAULT_GOES_DATA_DIR),
        omni_data_dir=event_config.get("omni_data_dir", DEFAULT_OMNI_DATA_DIR),
        emtf_data_dir=event_config.get("emtf_data_dir", DEFAULT_EMTF_DATA_DIR),
        efield_quiver_scale=(
            args.efield_quiver_scale
            if args.efield_quiver_scale is not None
            else map_config.get("efield_quiver_scale")
        ),
        efield_quiver_limit=(
            args.efield_quiver_limit
            if args.efield_quiver_limit is not None
            else map_config.get("efield_quiver_limit")
        ),
        case_station=args.case_station or event_config.get("case_station"),
    )
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
