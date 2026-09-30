#!/usr/bin/env python3
"""Create a 2x2 spatial MT apparent-resistivity/E-field response figure.

The figure contains:

    (a) determinant apparent resistivity at the selected resistivity period;
    (b) phase-referenced E-field response at the first requested period;
    (c) phase-referenced E-field response at the second requested period;
    (d) phase-referenced E-field response at the third requested period.

The E-field response is calculated in the frequency domain as

    E(f) = Z(f) B(f),

where ``Z`` is the EMTF impedance tensor and ``B`` is the Fourier component
of the nearest valid SuperMAG station during the event window.  The plotted
E vectors are the real part after rotating the complex phasor so that the
larger magnetic input component is real.  This gives a reproducible spatial
orientation while retaining the event-specific amplitude.
"""

from __future__ import annotations

import argparse
import datetime as dt
import glob
import os
import sys
from pathlib import Path

import bezpy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import toml


PROJECT_ROOT = Path(__file__).resolve().parent
for subdir in ("py", "py/fetch", "py/process", "py/plots"):
    sys.path.insert(0, str(PROJECT_ROOT / subdir))

from efield import CONUS_BOUNDS, _load_emtf_sites  # noqa: E402
from utils import calculate_GC_distance, detrend_magnetic_field  # noqa: E402

MU0 = 4.0e-7 * np.pi
DEFAULT_PERIODS = (60.0, 1000.0, 3600.0)
REGION_BOUNDS = {
    "conus": CONUS_BOUNDS,
    "north_america": (-135.0, -55.0, 20.0, 60.0),
    "europe": (-15.0, 40.0, 30.0, 72.0),
    "global": (-180.0, 180.0, -60.0, 80.0),
}


def parse_datetime(value):
    """Parse an ISO timestamp into a naive UTC datetime."""
    return dt.datetime.fromisoformat(str(value).replace("Z", "+00:00")).replace(
        tzinfo=None
    )


def parse_periods(value):
    """Parse comma-separated periods in seconds."""
    if isinstance(value, str):
        values = value.split(",")
    else:
        values = value
    periods = tuple(float(item) for item in values)
    if len(periods) != 3 or any(period <= 0.0 for period in periods):
        raise ValueError("Exactly three positive periods are required.")
    return periods


def normalize_longitude(value):
    return (float(value) + 180.0) % 360.0 - 180.0


def load_magnetic_data(csv_path, start, end):
    """Load and validate the cached SuperMAG CSV for the event."""
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Cached SuperMAG data were not found: {csv_path}\n"
            "Run event_study.py for this event first."
        )
    data = pd.read_csv(csv_path, parse_dates=["tval"])
    required = {"tval", "iaga", "glat", "glon", "N_geo", "E_geo", "Z_geo"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"{csv_path} is missing columns: {sorted(missing)}")
    data = data[(data.tval >= start) & (data.tval <= end)].copy()
    data = data.sort_values(["iaga", "tval"])
    return data


def station_catalog(data):
    """Return usable magnetic station records for nearest-site matching."""
    stations = {}
    for name, rows in data.groupby(data.iaga.astype(str), sort=False):
        rows = rows.copy()
        rows[["glat", "glon", "N_geo", "E_geo", "Z_geo"]] = rows[
            ["glat", "glon", "N_geo", "E_geo", "Z_geo"]
        ].apply(pd.to_numeric, errors="coerce")
        rows = rows.dropna(subset=["glat", "glon", "N_geo", "E_geo", "Z_geo"])
        if len(rows) < 8:
            continue
        first = rows.iloc[0]
        stations[name] = {
            "data": rows,
            "lat": float(first.glat),
            "lon": normalize_longitude(first.glon),
        }
    return stations


def nearest_station(site, stations):
    if not stations:
        return None
    return min(
        stations,
        key=lambda name: calculate_GC_distance(
            float(site.latitude),
            normalize_longitude(site.longitude),
            stations[name]["lat"],
            stations[name]["lon"],
        ),
    )


def magnetic_phasor(values, dT, period):
    """Return the complex Fourier amplitude nearest to ``1 / period``."""
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if good.sum() < 8:
        return np.nan + 1j * np.nan
    if not good.all():
        index = np.arange(values.size)
        values = np.interp(index, index[good], values[good])
    values = detrend_magnetic_field(values, L=min(120, values.size), dT=dT)
    window = np.hanning(values.size)
    if np.all(window == 0.0):
        return np.nan + 1j * np.nan
    nfft = max(4 * values.size, 1024)
    spectrum = np.fft.rfft(values * window, n=nfft)
    frequencies = np.fft.rfftfreq(nfft, d=dT)
    target = 1.0 / float(period)
    index = int(np.argmin(np.abs(frequencies - target)))
    # Convert the positive-frequency coefficient to an approximate sinusoidal
    # amplitude.  The factor is not used to normalize the plotted vectors.
    return 2.0 * spectrum[index] / window.sum()


def impedance_at_period(site, period):
    """Interpolate the EMTF impedance tensor at one period."""
    site_periods = np.asarray(site.periods, dtype=float)
    if period < np.nanmin(site_periods) or period > np.nanmax(site_periods):
        return None
    impedance = np.asarray(site.calcZ(np.array([1.0 / period])))
    if impedance.shape == (4,):
        impedance = impedance[:, None]
    if impedance.shape[0] != 4:
        return None
    return impedance[:, 0]


def determinant_resistivity(site, period):
    """Return determinant apparent resistivity in ohm-m."""
    impedance = impedance_at_period(site, period)
    if impedance is None:
        return np.nan
    # EMTF Z is [mV/km]/nT; multiplying by 1000 converts it to ohm.
    z_si = 1000.0 * impedance
    determinant = z_si[0] * z_si[3] - z_si[1] * z_si[2]
    omega = 2.0 * np.pi / float(period)
    rho = np.abs(determinant) / (MU0 * omega)
    return float(rho) if np.isfinite(rho) and rho > 0.0 else np.nan


def efield_at_period(site, magnetic_rows, period, dT, timestamp):
    """Calculate the period-specific E-field at ``timestamp`` in mV/km."""
    impedance = impedance_at_period(site, period)
    if impedance is None:
        return np.nan, np.nan
    rows = magnetic_rows.sort_values("tval")
    north = magnetic_phasor(rows.N_geo.to_numpy(), dT, period)
    east = magnetic_phasor(rows.E_geo.to_numpy(), dT, period)
    if not (np.isfinite(north) and np.isfinite(east)):
        return np.nan, np.nan
    e_north = impedance[0] * north + impedance[1] * east
    e_east = impedance[2] * north + impedance[3] * east
    elapsed = (timestamp - rows.tval.iloc[0]).total_seconds()
    time_phase = np.exp(1j * 2.0 * np.pi * elapsed / float(period))
    e_north = float(np.real(e_north * time_phase))
    e_east = float(np.real(e_east * time_phase))
    if not (np.isfinite(e_north) and np.isfinite(e_east)):
        return np.nan, np.nan
    return e_north, e_east


def efield_at_time(site, magnetic_rows, timestamp, dT):
    """Calculate the time-domain E-field nearest ``timestamp`` in mV/km."""
    rows = magnetic_rows.sort_values("tval").copy()
    if rows.empty:
        return np.nan, np.nan
    magnetic_north = detrend_magnetic_field(
        rows["N_geo"].to_numpy(dtype=float), dT=dT
    )
    magnetic_east = detrend_magnetic_field(
        rows["E_geo"].to_numpy(dtype=float), dT=dT
    )
    try:
        e_north, e_east = site.convolve_fft(
            magnetic_north, magnetic_east, dt=dT
        )
    except Exception:
        return np.nan, np.nan
    index = int(np.argmin(np.abs(rows["tval"] - timestamp)))
    if abs((rows.iloc[index]["tval"] - timestamp).total_seconds()) > 90:
        return np.nan, np.nan
    e_north = float(e_north[index])
    e_east = float(e_east[index])
    if not (np.isfinite(e_north) and np.isfinite(e_east)):
        return np.nan, np.nan
    return e_north, e_east


def load_sites(mt_data_folder, min_rating=3):
    sites = []
    for site in _load_emtf_sites(str(mt_data_folder)):
        try:
            rating = float(getattr(site, "rating", 0))
        except (TypeError, ValueError):
            rating = 0.0
        if rating < min_rating:
            continue
        sites.append(site)
    return sites


def build_response_records(
    data, sites, periods, rho_period, bounds, dT, peak_time
):
    """Build period-specific MT responses inside the map bounds."""
    del rho_period  # retained in the API for compatibility with the CLI.
    stations = station_catalog(data)
    records = []
    west, east, south, north = bounds
    for site in sites:
        lon = normalize_longitude(site.longitude)
        lat = float(site.latitude)
        if not (west <= lon <= east and south <= lat <= north):
            continue
        station_name = nearest_station(site, stations)
        if station_name is None:
            continue
        magnetic_rows = stations[station_name]["data"]
        rho = {
            period: determinant_resistivity(site, period)
            for period in periods
        }
        period_fields = {}
        for period in periods:
            period_fields[period] = efield_at_period(
                site, magnetic_rows, period, dT, peak_time
            )
        record = {
            "site": str(getattr(site, "name", "MT")),
            "site_object": site,
            "station": station_name,
            "lon": lon,
            "lat": lat,
            "rho": rho,
            "peak_E": efield_at_time(site, magnetic_rows, peak_time, dT),
            "E": period_fields,
        }
        if any(np.isfinite(value) for value in rho.values()) or any(
            np.isfinite(value)
            for pair in period_fields.values()
            for value in pair
        ):
            records.append(record)
    return records


def save_peak_estimates(records, periods, peak_time, output):
    """Save broadband and period-specific E-field estimates at the peak."""
    rows = []
    for record in records:
        e_north, e_east = record["peak_E"]
        row = {
            "event_time_utc": peak_time.strftime("%Y-%m-%dT%H:%M:%S"),
            "emtf_site": record["site"],
            "nearest_supermag": record["station"],
            "latitude": record["lat"],
            "longitude": record["lon"],
            "E_north_mV_per_km": e_north,
            "E_east_mV_per_km": e_east,
            "E_magnitude_mV_per_km": (
                float(np.hypot(e_north, e_east))
                if np.isfinite(e_north) and np.isfinite(e_east)
                else np.nan
            ),
        }
        for period in periods:
            p_north, p_east = record["E"][period]
            suffix = f"{period:g}s"
            row[f"E_north_{suffix}_mV_per_km"] = p_north
            row[f"E_east_{suffix}_mV_per_km"] = p_east
            row[f"E_magnitude_{suffix}_mV_per_km"] = (
                float(np.hypot(p_north, p_east))
                if np.isfinite(p_north) and np.isfinite(p_east)
                else np.nan
            )
        rows.append(row)
    peak_output = output.with_name(
        f"mt-efield-peak-{peak_time:%H%M}.csv"
    )
    peak_output.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame = frame.sort_values(
            "E_magnitude_mV_per_km", ascending=False
        )
    with peak_output.open("w", encoding="utf-8") as stream:
        stream.write(
            f"# EMTF E-field estimates at {peak_time:%Y-%m-%d %H:%M:%S} UTC\n"
        )
        stream.write(
            "# Broadband columns use time-domain convolution; period columns "
            "use E(f)=Z(f)B(f) evaluated at the same timestamp.\n"
        )
        stream.write("# All E-field values are in mV/km.\n")
        frame.to_csv(stream, index=False)
    return peak_output, frame


def map_bounds(region):
    key = str(region or "conus").lower()
    if key not in REGION_BOUNDS:
        raise ValueError(f"Unknown region '{region}'. Choose from {sorted(REGION_BOUNDS)}.")
    return REGION_BOUNDS[key]


def decorate_map(ax, bounds):
    ax.set_extent(bounds, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.97", zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
    ax.coastlines(resolution="50m", color="0.25", linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, edgecolor="0.3", linewidth=0.4)
    if hasattr(cfeature, "STATES"):
        ax.add_feature(cfeature.STATES, edgecolor="0.5", linewidth=0.25)
    ax.gridlines(draw_labels=False, linewidth=0.25, color="0.75", linestyle=":")


def plot_spatial_response(records, periods, rho_period, bounds, output, quiver_reference):
    """Render a 2x3 period-aligned resistivity/E-field figure."""
    del rho_period  # all requested periods are now shown in the top row.
    projection = ccrs.LambertConformal(
        central_longitude=(bounds[0] + bounds[1]) / 2.0,
        central_latitude=(bounds[2] + bounds[3]) / 2.0,
        standard_parallels=(bounds[2] + 10.0, bounds[3] - 10.0),
    )
    with plt.rc_context({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "savefig.dpi": 300,
    }):
        fig, axes = plt.subplots(
            2, 3, figsize=(12.0, 5.0), dpi=220,
            subplot_kw={"projection": projection},
        )
        axes = axes.ravel()
        letters = "abcdef"

        rho_values = np.asarray(
            [item["rho"][period] for item in records for period in periods],
            dtype=float,
        )
        rho_values = rho_values[np.isfinite(rho_values) & (rho_values > 0.0)]
        if rho_values.size:
            lower = 10.0 ** np.floor(np.log10(np.percentile(rho_values, 2)))
            upper = 10.0 ** np.ceil(np.log10(np.percentile(rho_values, 98)))
            if lower == upper:
                lower, upper = lower / 10.0, upper * 10.0
            rho_norm = LogNorm(vmin=lower, vmax=upper)
        else:
            rho_norm = LogNorm(vmin=1.0, vmax=1.0e5)

        rho_mappable = None
        for column, period in enumerate(periods):
            ax = axes[column]
            decorate_map(ax, bounds)
            rho_records = [
                item for item in records
                if np.isfinite(item["rho"][period])
            ]
            if rho_records:
                rho_mappable = ax.scatter(
                    [item["lon"] for item in rho_records],
                    [item["lat"] for item in rho_records],
                    c=[item["rho"][period] for item in rho_records],
                    norm=rho_norm, cmap="viridis", s=10,
                    edgecolors="none", transform=ccrs.PlateCarree(), zorder=3,
                )
            exponent = int(np.floor(np.log10(1000.0 / period)))
            mantissa = (1000.0 / period) / (10.0 ** exponent)
            frequency_label = rf"$f={mantissa:.3g}\times10^{{{exponent}}}$ mHz"
            ax.text(
                0.025, 0.965, f"({letters[column]})",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=11, fontweight="bold", zorder=10,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 1.5},
            )
            ax.text(
                0.05, 1.05, f"T={period:g} s\n{frequency_label}",
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=9, zorder=10,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 1.5},
            )

        if rho_mappable is not None:
            colorbar_axis = fig.add_axes([0.905, 0.565, 0.014, 0.285])
            colorbar = fig.colorbar(rho_mappable, cax=colorbar_axis)
            colorbar.set_label(r"$\rho_a$ ($\Omega$ m)", labelpad=8)
            colorbar_axis.tick_params(labelsize=8, direction="in")

        for column, period in enumerate(periods):
            ax = axes[3 + column]
            decorate_map(ax, bounds)
            vectors = [
                item for item in records
                if np.isfinite(item["E"][period][0])
                and np.isfinite(item["E"][period][1])
            ]
            quiver = None
            if vectors:
                quiver = ax.quiver(
                    np.asarray([item["lon"] for item in vectors]),
                    np.asarray([item["lat"] for item in vectors]),
                    np.asarray([item["E"][period][1] for item in vectors])
                    / quiver_reference,
                    np.asarray([item["E"][period][0] for item in vectors])
                    / quiver_reference,
                    transform=ccrs.PlateCarree(),
                    color="#d62728", pivot="tail", scale=8,
                    scale_units="width", width=0.0022,
                    headwidth=3.5, headlength=4.5, headaxislength=3.5,
                    zorder=4,
                )
            exponent = int(np.floor(np.log10(1000.0 / period)))
            mantissa = (1000.0 / period) / (10.0 ** exponent)
            frequency_label = rf"$f={mantissa:.3g}\times10^{{{exponent}}}$ mHz"
            ax.text(
                0.025, 0.965, f"({letters[3 + column]})",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=11, fontweight="bold", zorder=10,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.7, "pad": 1.5},
            )
            if column == 0 and quiver is not None:
                key = ax.quiverkey(
                    quiver, 0.84, 0.76, 1.0,
                    label=rf"$\vec{{\mathcal{{E}}}}$ = "
                    rf"{quiver_reference:g} mV/km",
                    labelpos="S", coordinates="axes", color="#d62728",
                    fontproperties={"size": 9},
                )
                key.text.set_bbox({
                    "facecolor": "white", "edgecolor": "none",
                    "alpha": 0.85, "pad": 1.5,
                })
        ax = axes[5]
        ax.text(
            1.02, 0.9, 
            "7 Sept. 2017, 17:40 UT", ha="left", va="top", 
            transform=ax.transAxes, 
            fontsize=11, fontweight="bold", zorder=10,
            rotation=90
        )
        fig.subplots_adjust(left=0.025, right=0.88, bottom=0.035, top=0.965,
                            wspace=0.018, hspace=0.0)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
        plt.close(fig)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", required=True,
        help="Event-study TOML configuration file.",
    )
    parser.add_argument(
        "--periods", default=None,
        help="Three periods in seconds, e.g. 60,1000,3600.",
    )
    parser.add_argument(
        "--rho-period", type=float, default=None,
        help="Period in seconds for the apparent-resistivity panel.",
    )
    parser.add_argument("--region", default=None)
    parser.add_argument("--min-rating", type=float, default=3.0)
    parser.add_argument("--dT", type=float, default=60.0)
    parser.add_argument("--quiver-reference", type=float, default=250.0)
    parser.add_argument(
        "--peak-time", default=None,
        help="Peak timestamp as ISO datetime or HH:MM; defaults to event time.",
    )
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    config_path = Path(args.config).expanduser().resolve()
    config = toml.load(config_path)
    event_config = config.get("event", {})
    map_config = config.get("map", {})
    output_config = config.get("output", {})
    plot_config = config.get("plot", {})

    event = parse_datetime(event_config["time"])
    start = parse_datetime(event_config.get("start", event_config["time"]))
    end = parse_datetime(event_config.get("end", event_config["time"]))
    periods = parse_periods(
        args.periods or plot_config.get("mt_response_periods", DEFAULT_PERIODS)
    )
    if args.peak_time:
        if "T" in args.peak_time or "-" in args.peak_time:
            peak_time = parse_datetime(args.peak_time)
        else:
            peak_clock = dt.datetime.strptime(args.peak_time, "%H:%M").time()
            peak_time = event.replace(
                hour=peak_clock.hour, minute=peak_clock.minute,
                second=0, microsecond=0,
            )
    else:
        peak_time = event
    rho_period = args.rho_period
    if rho_period is None:
        rho_period = float(plot_config.get("mt_rho_period", periods[1]))
    if rho_period <= 0.0:
        raise ValueError("--rho-period must be positive.")
    region = args.region or map_config.get("region", "conus")
    bounds = map_bounds(region)

    data_dir = PROJECT_ROOT / "tmp" / event.strftime("%Y-%m-%d-%H-%M")
    magnetic_csv = data_dir / "supermag.csv"
    magnetic_data = load_magnetic_data(magnetic_csv, start, end)
    sites = load_sites(PROJECT_ROOT / ".data" / "EMTF", args.min_rating)
    records = build_response_records(
        magnetic_data, sites, periods, rho_period, bounds, args.dT, peak_time
    )
    if not records:
        raise RuntimeError("No usable EMTF sites produced a spatial response.")

    if args.output:
        output = Path(args.output).expanduser()
        if not output.is_absolute():
            output = PROJECT_ROOT / output
    else:
        output_root = Path(output_config.get("root", "pub_fig"))
        if not output_root.is_absolute():
            output_root = PROJECT_ROOT / output_root
        output = output_root / event.strftime("%Y-%m-%d-%H-%M") / "mt-efield-spatial.png"

    plot_spatial_response(
        records, periods, rho_period, bounds, output, args.quiver_reference
    )
    peak_output, peak_frame = save_peak_estimates(
        records, periods, peak_time, output
    )
    print(f"output: {output}")
    print(f"peak E-field CSV: {peak_output}")
    print(f"peak timestamp: {peak_time:%Y-%m-%d %H:%M:%S} UTC")
    print(f"peak E-field rows: {len(peak_frame)}")
    if not peak_frame.empty:
        print(
            "largest peak E-field: "
            f"{peak_frame.iloc[0]['emtf_site']} "
            f"{peak_frame.iloc[0]['E_magnitude_mV_per_km']:.3f} mV/km"
        )
    print(f"sites: {len(records)}")
    print(f"periods: {', '.join(f'{period:g} s' for period in periods)}")
    print(f"rho period: {rho_period:g} s")


if __name__ == "__main__":
    main()
