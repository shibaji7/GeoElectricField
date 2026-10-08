"""Induced electric-field plots for event studies."""

import glob
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import os

import bezpy
from utils import calculate_GC_distance, detrend_magnetic_field
from utils import compute_Egeo_from_Bgeo
from maps import Map


EMTF_DATA_FOLDER = ".data/EMTF"
CONUS_BOUNDS = (-125.0, -66.9, 24.4, 49.4)


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
                    case_efield_exclude_sites=None):
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

    case_data = data_by_station[case_station]
    magnetic_north = detrend_magnetic_field(case_data["N_geo"].to_numpy(), dT=dT)
    magnetic_east = detrend_magnetic_field(case_data["E_geo"].to_numpy(), dT=dT)
    
    magnetic_north = magnetic_north - np.median(magnetic_north[0:30])
    magnetic_east = magnetic_east - np.median(magnetic_east[0:30])
            
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
    if site is None:
        raise RuntimeError(f"No EMTF response found for case station '{case_station}'.")

    # Calculate a second E-field response at VAQ55 using the same magnetic
    # station and magnetic perturbation time series.
    vaq55_site = None
    e_north_55 = e_east_55 = None
    site_name = str(getattr(site, "name", "")).upper()
    if site_name != "VAQ55":
        available_sites = candidates if case_mt_site else _load_emtf_sites(mt_data_folder)
        vaq55_site = next(
            (candidate for candidate in available_sites
             if str(getattr(candidate, "name", "")).upper() == "VAQ55"),
            None,
        )
        if vaq55_site is None:
            raise RuntimeError(f"EMTF site 'VAQ55' was not found under {mt_data_folder}.")
        e_north_55, e_east_55 = vaq55_site.convolve_fft(
            magnetic_north, magnetic_east, dt=dT
        )

    B = np.sqrt(magnetic_north ** 2 + magnetic_east ** 2)
    B = np.nanmean(B[: min(20, len(B))]) - B

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
        if case_xlim is None:
            case_xlim = (
                event.replace(hour=17, minute=0, second=0, microsecond=0),
                event.replace(hour=20, minute=0, second=0, microsecond=0),
            )
        ax_top.set_xlim(case_xlim)
        mag_ax = ax_top.twinx()
        for station in top_stations:
            data = data_by_station[station]
            north = detrend_magnetic_field(data["N_geo"].to_numpy(), dT=dT)
            east = detrend_magnetic_field(data["E_geo"].to_numpy(), dT=dT)
            field = np.sqrt(north ** 2 + east ** 2)
            field = np.nanmean(field[: min(20, len(field))]) - field
            mag_ax.plot(
                data.tval,
                field,
                lw=0.9,
                ls="-",
                color=station_colors.get(station, "#333333"),
                label=station,
            )
        mag_ax.set_ylabel(r"$\Delta B$ (nT)")
        mag_ax.set_ylim(-150, 50)
        mag_ax.axvline(event, color="0.15", lw=0.7)
        mag_ax.legend(loc="lower right", ncol=3, frameon=False,
                      handlelength=1.6, columnspacing=1.0)
        ax_top.text(0.015, 0.94, "(a)", transform=ax_top.transAxes, fontweight="bold")

        ax_bottom.grid(axis="y", linestyle=":", color="0.55", linewidth=0.45)
        e_magnitude = np.hypot(e_east, e_north)
        primary_site_name = str(getattr(site, "name", None) or case_mt_site or "MT").upper()
        ax_bottom.plot(
            case_data.tval,
            e_east,
            color="#1f77b4",
            lw=1.0,
            ls="-",
            label=rf"{primary_site_name}",
        )
        ax_bottom.plot(
            case_data.tval,
            e_north,
            color="#1fb444",
            lw=1.0,
            ls="-",
            label=rf"{primary_site_name}",
        )
        magnitude_series = [e_magnitude]
        legend_handles = [
            plt.Line2D([], [], color="#1f77b4",
                       label=rf"$|\mathcal{{E}}|$ ({primary_site_name})")
        ]
        # if vaq55_site is not None:
        #     e_magnitude_55 = np.hypot(e_east_55, e_north_55)
        #     ax_bottom.plot(
        #         case_data.tval,
        #         e_magnitude_55,
        #         color="#d62728",
        #         lw=1.0,
        #         ls="-",
        #         label=r"VAQ55",
        #     )
        #     magnitude_series.append(e_magnitude_55)
        #     legend_handles.append(
        #         plt.Line2D([], [], color="#d62728", label=r"$|\mathcal{E}|$ (VAQ55)")
        #     )
        ax_bottom.set_ylabel(r"$|\mathcal{E}|$ (mV/km)")
        finite_magnitude = np.concatenate([
            np.asarray(series)[np.isfinite(series)] for series in magnitude_series
        ])
        ax_bottom.set_ylim(-1000, 2000)
        positive_magnitude = finite_magnitude[finite_magnitude > 0.0]
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
            loc="upper left",
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
            reference = 500.0  # mV/km
            excluded_sites = {
                str(name).upper() for name in (case_efield_exclude_sites or [])
            }
            for site_item, station_name, data in map_fields:
                if data.empty or site_item.rating < 3:
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
                if (case_efield_quiver_limit is not None and
                        magnitude > float(case_efield_quiver_limit)):
                    continue
                site_lon = (float(site_item.longitude) + 180.0) % 360.0 - 180.0
                site_lat = float(site_item.latitude)
                # Compare in physical units before normalizing by the
                # 500 mV/km reference used by the quiver scale.
                lon.append(site_lon)
                lat.append(site_lat)
                ex.append(east_value / reference)
                ey.append(north_value / reference)
                
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
