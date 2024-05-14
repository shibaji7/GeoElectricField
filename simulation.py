#!/usr/bin/env python3

"""simulation.py: simulate python program for RT"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import numpy as np
import sys
sys.path.extend(["py/", "py/fetch/", "py/process/", "py/plots/"])
import os
os.environ["OMNIDATA_PATH"] = "/home/shibaji/omni/"

import mplstyle
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from get_mag_data import SuperMAG
from goes import FlareTS

from pyomnidata import GetOMNI
import DateTimeTools as TT

def setup(base):
    os.makedirs(f"{base}/figures/mag/", exist_ok=True)
    os.makedirs(f"{base}/figures/elec/", exist_ok=True)
    os.makedirs(f"{base}/figures/map/", exist_ok=True)
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
            "font.weight": "bold"
        }
    )
    return

def stack_plots(ev, g, smag, omni, dates, fname, hr=6):
    fig = plt.figure(figsize=(6, 6*2), dpi=240)
    ax = fig.add_subplot(611)
    ax.set_ylabel(r"Irradiance ($W/m^2$)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.text(0.1, 1.05, dates[0].strftime("%d-")+dates[1].strftime("%d, %b %Y"), ha="left", va="center", transform=ax.transAxes)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.set_xlim(dates)
    ax.semilogy(
        g.dfs["goes"].time,
        g.dfs["goes"].xrsa,
        marker="o",
        color="b",
        ls="None",
        ms=0.3,
        alpha=0.9,
        label=rf"$\lambda_0\sim (0.05-0.4)$ nm",
    )
    ax.semilogy(
        g.dfs["goes"].time,
        g.dfs["goes"].xrsb,
        marker="o",
        color="r",
        ls="None",
        ms=0.3,
        alpha=0.9,
        label=rf"$\lambda_0\sim (0.1-0.8)$ nm",
    )
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_ylim(1e-8, 1e-2)
    ax.axhline(1e-4, color="r", ls="--", lw=0.6, alpha=0.7)
    ax.axhline(1e-5, color="orange", ls="--", lw=0.6, alpha=0.7)
    ax.axhline(1e-6, color="darkgreen", ls="--", lw=0.6, alpha=0.7)
    ax.text(dates[0] + dt.timedelta(minutes=5), 3e-6, "C")
    ax.text(dates[0] + dt.timedelta(minutes=5), 3e-5, "M")
    ax.text(dates[0] + dt.timedelta(minutes=5), 3e-4, "X")
    ax.legend(loc=1)
    ax.set_xticklabels([])

    utc = TT.Datetime(omni.Date,omni.ut)
    ax = fig.add_subplot(612)
    ax.set_ylabel(r"$B_{y,z}$ (nT)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["BzGSM"], "bo", ms=0.4, label=r"$B_z$")
    ax.plot(utc, omni["ByGSM"], "ro", ms=0.4, label=r"$B_y$")
    ax.axhline(0, lw=0.3, ls="--", color="k")
    ax.set_ylim(-20, 20)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_xlim(dates)
    ax.legend(loc=1)
    ax.set_xticklabels([])

    ax = fig.add_subplot(613)
    ax.set_ylabel(r"Proton Density (cm$^{-3}$)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["ProtonDensity"], "ko", ms=0.4)
    ax.set_ylim(0, 20)
    ax.set_xlim(dates)
    ax = ax.twinx()
    ax.set_ylabel(r"Flow Speed (km s$^{-1}$)", fontdict=dict(color="blue"))
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.plot(utc, omni["FlowSpeed"], "bo", ms=0.4)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.set_ylim(400, 1200)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_xticklabels([])

    ax = fig.add_subplot(614)
    ax.set_ylabel(r"Dynamic Pressure (nPa)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["FlowPressure"], "ko", ms=0.4)
    ax.set_ylim(0, 20)
    ax.set_xlim(dates)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_xticklabels([])

    ax = fig.add_subplot(615)
    ax.set_ylabel(r"Sym-H (nT)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["SymH"], "ko", ms=0.4)
    ax.set_ylim(-200, 150)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_xlim(dates)
    ax = ax.twinx()
    ax.set_ylabel(r"Asy-H (nT)", fontdict=dict(color="blue"))
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.plot(utc, omni["AsyH"], "bo", ms=0.4)
    ax.set_ylim(0, 150)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.set_xticklabels([])

    ax = fig.add_subplot(616)
    ax.set_ylabel(r"AL/AU/AE (nT)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["AL"], "ko", ms=0.4, label="AL")
    ax.plot(utc, omni["AU"], "bo", ms=0.4, label="AU")
    ax.plot(utc, omni["AE"], "ro", ms=0.4, label="AE")
    ax.legend(loc=2)
    ax.set_ylim(-2000, 2000)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_xlim(dates)
    ax.set_xlabel("Time (UT)")
    fig.savefig(fname, bbox_inches="tight")

    
    fig = plt.figure(figsize=(6, 6*2), dpi=240)
    for i, stn in enumerate(smag.stations[:6]):
        data = smag.sm_data[smag.sm_data.iaga == stn]
        ax = fig.add_subplot(912+i)
        ax.set_ylabel(r"$\delta$ (nT)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
        ax.set_xlim(dates)
        if i<5: ax.set_xticklabels([])
        ax.axhline(0, ls="--", lw=0.4, color="k")
        for comp, col in zip(["N", "E", "Z"], ["r", "g", "b"]):
            ax.plot(
                data.tval,
                data[comp + "_geo"],
                ls="-",
                color=col,
                lw=0.6,
                label=comp+r"$_{geo}$",
            )
        if i==0: 
            ax.legend(loc=1, fontsize=8)
            ax.text(0.1, 1.05, dates[0].strftime("%d-")+dates[1].strftime("%d, %b %Y"), ha="left", va="center", transform=ax.transAxes)
        ax.set_ylim(-150, 150)
        txt = f"{stn}" + "\n" + rf"$\theta_m,\phi_m$={'%.1f'%data.mlat.tolist()[0]}, {'%.1f'%data.mlon.tolist()[0]}"
                #"\n"+ fr"$\chi$={'%.2f'%np.mean(data.sza)}"
                #rf"$\theta_g,\phi_g$={data.glat.tolist()[0]}, {data.glon.tolist()[0]}" + \
        ax.axvline(ev, ls="-", lw=0.5, color="k")
        ax.text(0.3, 0.2, txt, ha="center", va="center", transform=ax.transAxes, fontdict={"size":8})
        ax = ax.twinx()
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
        ax.set_xlim(dates)
        ax.plot(
            data.tval,
            data["dBdt"],
            ls="-",
            color="k",
            lw=0.4,
        )
        if i==3: ax.set_ylabel(r"$\frac{\partial |B|}{\partial t}$, (nT $s^-1$)")
        if i<5: ax.set_xticklabels([])
        ax.set_ylim(-0.1,0.1)
    ax.set_xlabel("Time (UT)")
    fig.savefig(fname.replace(".png", "-mag.png"), bbox_inches="tight")
    return


def getOMNI(year, Date, ut=[0.0,24.0]):
    data = GetOMNI(year,1)
    use = np.where(((data.Date == Date[0]) & (data.ut >= ut[0])) |
						((data.Date == Date[1]) & (data.ut <= ut[1])) |
						((data.Date > Date[0]) & (data.Date < Date[1])))
    data = data[use]
    return data

def run_create_stack_plots(
    ev, dates, mag_stations=[], 
    omni=True, gef_station_map=dict(),
    calculate_gef=False,
    NERC_benchmark=False,
):
    base = "tmp/{Y}-{m}-{d}-{H}-{M}/".format(
        Y=ev.year,
        m="%02d" % ev.month,
        d="%02d" % ev.day,
        H="%02d" % ev.hour,
        M="%02d" % ev.minute,
    )
    setup(base)
    g = FlareTS(dates)
    smag = SuperMAG(
        base, dates, uid="shibaji7",
        stations=mag_stations
    )
    omni = getOMNI(
        ev.year,
        [
            int(dates[0].strftime("%Y%m%d")),
            int(dates[1].strftime("%Y%m%d")),
        ]
    )
    g.plot_TS()
    g.save("goes.png", f"{base}/figures/")
    g.close()

    stack_plots(ev, g, smag, omni, dates, f"{base}/figures/stack.png")
    return

if __name__ == "__main__":
    ev, dates = (
        dt.datetime(2003,11,4,19,53), 
        [
            dt.datetime(2003,11,4),
            dt.datetime(2003,11,5)
        ]
    )
    ev, dates = (
        dt.datetime(2005,9,7,17,40), 
        [
            dt.datetime(2005,9,7,16),
            dt.datetime(2005,9,7,19)
        ]
    )
    mag_stations = [
        "BOU", "FRD", "FRN", 
        "GLN", "MEA", "NEW", 
        "OTT", "PIN", "VIC"
    ]
    run_create_stack_plots(
        ev, dates, mag_stations
    )