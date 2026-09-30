"""Publication-style plots for event studies."""

import datetime as dt
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {"GOES-A": "#7b3294", "GOES-B": "#c51b7d", "By": "#4c78a8", "Bz": "#f58518"}


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.linewidth": 0.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "legend.frameon": False,
            "savefig.dpi": 300,
            "figure.dpi": 150,
        }
    )


def _format_axis(ax):
    ax.tick_params(which="both", top=True, right=True, width=0.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)


def _omni_times(omni):
    values = []
    for date, hour in zip(omni.Date, omni.ut):
        values.append(
            dt.datetime.strptime(str(int(date)), "%Y%m%d")
            + dt.timedelta(hours=float(hour))
        )
    return values


def plot_stack(goes, omni, event, dates, output):
    """Plot GOES and OMNI event context in the publication style."""
    apply_style()
    fig, axes = plt.subplots(3, 1, figsize=(6.6, 6.2), sharex=True, constrained_layout=True)
    goes_df = goes.dfs["goes"]

    ax = axes[0]
    ax.semilogy(goes_df.time, goes_df.xrsa, ".", ms=1.0, color=COLORS["GOES-A"], label="GOES A")
    ax.semilogy(goes_df.time, goes_df.xrsb, ".", ms=1.0, color=COLORS["GOES-B"], label="GOES B")
    ax.set_ylabel(r"Irradiance (W m$^{-2}$)")
    ax.set_ylim(1e-8, 1e-2)
    ax.axvline(event, color="0.15", lw=0.7)
    ax.axhline(1e-6, color="0.5", lw=0.4, ls=":")
    ax.axhline(1e-5, color="0.5", lw=0.4, ls=":")
    ax.axhline(1e-4, color="0.5", lw=0.4, ls=":")
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.text(0.015, 0.92, "(A)", transform=ax.transAxes, fontweight="bold")

    utc = _omni_times(omni)
    ax = axes[1]
    ax.plot(utc, omni["ByGSM"], ".", ms=1.0, color=COLORS["By"], label=r"$B_y$")
    ax.plot(utc, omni["BzGSM"], ".", ms=1.0, color=COLORS["Bz"], label=r"$B_z$")
    ax.axhline(0, color="0.3", lw=0.4, ls=":")
    ax.set_ylabel("IMF (nT)")
    ax.set_ylim(-30, 30)
    ax.legend(loc="upper left", fontsize=8, ncol=2)
    ax.text(0.015, 0.88, "(B)", transform=ax.transAxes, fontweight="bold")

    ax = axes[2]
    ax.plot(utc, omni["SymH"], ".", ms=1.0, color="0.2", label="Sym-H")
    ax.plot(utc, omni["AsyH"], ".", ms=1.0, color="#d95f02", label="Asy-H")
    ax.plot(utc, omni["AE"], ".", ms=1.0, color="#1b9e77", label="AE")
    ax.set_ylabel("Indices (nT)")
    ax.set_xlabel("UT")
    ax.set_ylim(-250, 250)
    ax.legend(loc="upper left", fontsize=8, ncol=3)
    ax.text(0.015, 0.88, "(C)", transform=ax.transAxes, fontweight="bold")

    for ax in axes:
        ax.set_xlim(dates)
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        _format_axis(ax)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_magnetic(smag, event, dates, stations, output):
    """Plot ordered SuperMAG magnetic components and |dB|/dt."""
    apply_style()
    fig, axes = plt.subplots(2, 1, figsize=(6.6, 5.2), sharex=True, constrained_layout=True)
    for station in stations:
        data = smag.sm_data[smag.sm_data.iaga == station]
        if data.empty:
            continue
        axes[0].plot(data.tval, data.N_geo, lw=0.7, label=f"{station} N")
        axes[0].plot(data.tval, data.E_geo, lw=0.7, ls="--", label=f"{station} E")
        axes[1].plot(data.tval, data.dBdt, lw=0.7, label=station)

    axes[0].set_ylabel(r"$ΔB_{geo}$ (nT)")
    axes[1].set_ylabel(r"$d|B|/dt$ (nT s$^{-1}$)")
    axes[1].set_xlabel("UT")
    axes[0].text(0.015, 0.9, "(A)", transform=axes[0].transAxes, fontweight="bold")
    axes[1].text(0.015, 0.9, "(B)", transform=axes[1].transAxes, fontweight="bold")
    axes[0].legend(loc="upper left", fontsize=7, ncol=3)
    axes[1].legend(loc="upper left", fontsize=7, ncol=6)
    for ax in axes:
        ax.set_xlim(dates)
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        _format_axis(ax)
    axes[0].axvline(event, color="0.15", lw=0.7)
    axes[1].axvline(event, color="0.15", lw=0.7)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def _map_settings(region):
    if region == "north_polar":
        return ccrs.NorthPolarStereo(central_longitude=0), [-180, 180, 30, 90]
    if region == "europe":
        return ccrs.Stereographic(central_longitude=10, central_latitude=55), [-30, 45, 30, 80]
    return ccrs.Stereographic(central_longitude=-90, central_latitude=45), [-130, -70, 30, 80]


def plot_induced_field_map(smag, dates, map_times, region, output):
    """Plot all available SuperMAG induced-field vectors at selected times."""
    apply_style()
    projection, extent = _map_settings(region)
    fig, axes = plt.subplots(
        2, 3, figsize=(7.2, 5.3), subplot_kw={"projection": projection}, constrained_layout=True
    )
    axes = np.asarray(axes).ravel()
    valid_times = list(map_times)[:6]
    for index, ax in enumerate(axes):
        if index >= len(valid_times):
            ax.set_visible(False)
            continue
        timestamp = valid_times[index]
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5, edgecolor="0.3")
        ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.35, edgecolor="0.45")
        ax.gridlines(draw_labels=False, linewidth=0.3, color="0.7", linestyle=":")
        rows = smag.sm_data[smag.sm_data.tval == timestamp].dropna(
            subset=["glon", "glat", "E_geo", "N_geo"]
        )
        if not rows.empty:
            q = ax.quiver(
                rows.glon.to_numpy(), rows.glat.to_numpy(),
                rows.E_geo.to_numpy() / 500.0, rows.N_geo.to_numpy() / 500.0,
                transform=ccrs.PlateCarree(), color="#d62728", pivot="tail",
                scale=8, scale_units="width", width=0.002, headwidth=4, headlength=5,
            )
            if index == min(len(valid_times), 6) - 1:
                ax.quiverkey(q, 0.88, 1.04, 1.0, r"$E$=500 mV/km", labelpos="N", coordinates="axes")
        ax.set_title(timestamp.strftime("%H:%M UT"), fontsize=9)
        ax.text(0.02, 0.96, f"({chr(65 + index)})", transform=ax.transAxes, fontweight="bold")
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
