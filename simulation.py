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
    from matplotlib import rc
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

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

def stack_plots(ev, g, smag, omni, dates, fname, hr=1):
    fig = plt.figure(figsize=(6, 6*2), dpi=240)
    ax = fig.add_subplot(611)
    ax.set_ylabel(r"Irradiance ($W/m^2$)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.text(0.1, 1.05, dates[0].strftime("%d, %b %Y"), ha="left", va="center", transform=ax.transAxes)
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
    ax.text(0.1, 0.2, "(a)", ha="left", va="center", transform=ax.transAxes)
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
    ax.text(0.1, 0.2, "(b)", ha="left", va="center", transform=ax.transAxes)
    ax.set_xticklabels([])

    ax = fig.add_subplot(613)
    ax.set_ylabel(r"Proton Density (cm$^{-3}$)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["ProtonDensity"], "ko", ms=0.4)
    ax.set_ylim(0, 10)
    ax.set_xlim(dates)
    ax = ax.twinx()
    ax.set_ylabel(r"Flow Speed (km s$^{-1}$)", fontdict=dict(color="blue"))
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.plot(utc, omni["FlowSpeed"], "bo", ms=0.4)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.set_ylim(200, 500)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.text(0.1, 0.2, "(c)", ha="left", va="center", transform=ax.transAxes)
    ax.set_xticklabels([])

    ax = fig.add_subplot(614)
    ax.set_ylabel(r"Dynamic Pressure (nPa)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["FlowPressure"], "ko", ms=0.4)
    ax.set_ylim(0, 5)
    ax.set_xlim(dates)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.text(0.1, 0.2, "(d)", ha="left", va="center", transform=ax.transAxes)
    ax.set_xticklabels([])

    ax = fig.add_subplot(615)
    ax.set_ylabel(r"Sym-H (nT)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["SymH"], "ko", ms=0.4)
    ax.set_ylim(-50, 50)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_xlim(dates)
    ax = ax.twinx()
    ax.set_ylabel(r"Asy-H (nT)", fontdict=dict(color="blue"))
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.plot(utc, omni["AsyH"], "bo", ms=0.4)
    ax.set_ylim(0, 150)
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.text(0.1, 0.2, "(e)", ha="left", va="center", transform=ax.transAxes)
    ax.set_xticklabels([])

    ax = fig.add_subplot(616)
    ax.set_ylabel(r"AL/AU/AE (nT)")
    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, hr)))
    ax.plot(utc, omni["AL"], "ko", ms=0.4, label="AL")
    ax.plot(utc, omni["AU"], "bo", ms=0.4, label="AU")
    ax.plot(utc, omni["AE"], "ro", ms=0.4, label="AE")
    ax.legend(loc=2)
    ax.set_ylim(-200, 200)
    ax.axvline(ev, ls="-", lw=1.2, color="k")
    ax.set_xlim(dates)
    ax.set_xlabel("Time (UT)")
    ax.text(0.1, 0.2, "(f)", ha="left", va="center", transform=ax.transAxes)
    fig.savefig(fname, bbox_inches="tight")
    fig.savefig(fname.replace(".png", ".eps"), bbox_inches="tight")

    
    fig = plt.figure(figsize=(6, 18), dpi=240)
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
            ax.text(0.1, 1.05, dates[1].strftime("%d, %b %Y"), ha="left", va="center", transform=ax.transAxes)
        ax.set_ylim(-150, 150)
        if i==5:
            ax.set_xlabel("Time (UT)")
        txt = f"({chr(97+i)}) {stn}" + "\n" + rf"$\theta_m,\phi_m$={'%.1f'%data.mlat.tolist()[0]}, {'%.1f'%data.mlon.tolist()[0]}"
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
        if i==3: 
            ax.text(
                1.1, 1.0,
                r"$\frac{\partial |B|}{\partial t}$, (nT $s^-1$)",
                rotation=90, transform=ax.transAxes, ha="center", va="center"
            )
        if i<5: ax.set_xticklabels([])
        ax.set_ylim(-0.1,0.1)
        
    fig.savefig(fname.replace(".png", "-mag.png"), bbox_inches="tight")
    fig.savefig(fname.replace(".png", "-mag.eps"), bbox_inches="tight")
    return


def getOMNI(year, Date, ut=[0.0,24.0]):
    data = GetOMNI(year,1)
    use = np.where(((data.Date == Date[0]) & (data.ut >= ut[0])) |
						((data.Date == Date[1]) & (data.ut <= ut[1])) |
						((data.Date > Date[0]) & (data.Date < Date[1])))
    data = data[use]
    return data

def fetch_USGS_database(elec_stations, base, dates):
    """
    Load USGS dataset
    """
    from utils import get_nearest_station
    from get_efields import (
        download_save_EarthScope_GeoEB,
        load_MT_site,
        load_USArray_one_station_multi_cha,
    )
    from utils import compute_Egeo_from_Bgeo, detrend_magnetic_field

    usgs = {}
    usgs_sites = {}
    return
    

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
    #g = FlareTS(dates)
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
    #g.plot_TS()
    #g.save("goes.png", f"{base}/figures/")
    #g.close()

    #stack_plots(ev, g, smag, omni, dates, f"{base}/figures/stack.png")

    from maps import Map
    dates = [
        dates[0] + dt.timedelta(minutes=i) 
        for i in range(int((dates[1]-dates[0]).total_seconds()/60))
    ]
    dates =[
        dt.datetime(2005, 9, 7, 17, 20),
        dt.datetime(2005, 9, 7, 17, 40),
        dt.datetime(2005, 9, 7, 18),
        dt.datetime(2005, 9, 7, 18, 20),
        dt.datetime(2005, 9, 7, 19),
        dt.datetime(2005, 9, 7, 19, 30),
    ]
    # for i, d in enumerate(dates):
    #     map = Map(d, fig_title=f"({chr(97+i)})")
    #     for i, sta in enumerate(smag.stations):
    #         o = smag.sm_data[
    #             (smag.sm_data.tval == d)
    #             & (smag.sm_data.iaga == sta)
    #         ]
    #         if len(o) > 0:
    #             inst = dict(
    #                 lon=np.mod(o.glon.tolist()[0]+180,360)-180,
    #                 lat=o.glat.tolist()[0],
    #                 code=o.iaga.tolist()[0],
    #                 E_geo=o.E_geo.tolist()[0],
    #                 N_geo=o.N_geo.tolist()[0],
    #                 Z_geo=o.Z_geo.tolist()[0],
    #             )
    #             map.plot_instrument(inst)
    #     map.save(f"{base}/figures/map/map_{d.strftime('%H%M')}.png")
    #     #break
    map = Map(dates[0], nrows=3, ncols=2)
    map.ax.text(
        0.1, 1.05,
        f"Coords: Geo",
        ha="left",
        fontweight="bold",
        fontsize=12,
        transform=map.ax.transAxes,
    )
    ax = map.ax
    for i, d in enumerate(dates):
        ax.text(
            0.1, 0.9, f"({chr(97+i)}) / {d.strftime('%H:%M UT')}", 
            ha="left", va="center",
            fontsize=12,
            transform=ax.transAxes,
        )
        for j, sta in enumerate(smag.stations):
            o = smag.sm_data[
                (smag.sm_data.tval == d)
                & (smag.sm_data.iaga == sta)
            ]
            if len(o) > 0:
                inst = dict(
                    lon=np.mod(o.glon.tolist()[0]+180,360)-180,
                    lat=o.glat.tolist()[0],
                    code=o.iaga.tolist()[0],
                    E_geo=o.E_geo.tolist()[0],
                    N_geo=o.N_geo.tolist()[0],
                    Z_geo=o.Z_geo.tolist()[0],
                )
                map.plot_instrument(inst)
        if i<5: 
            ax = map.add_axes()
            map.ax = ax
        if i==0:
            map.ax.text(
                0.1, 1.05,
                f"Date: 07, Sep 2005",
                ha="left",
                fontweight="bold",
                fontsize=12,
                transform=map.ax.transAxes,
            )
    map.fig.subplots_adjust(wspace=0.05, hspace=0.05)
    map.save(f"{base}/figures/map.png")
    return

def run_mag_e_stack_plots(ev, dates, mag_stations):
    base = "tmp/{Y}-{m}-{d}-{H}-{M}/".format(
        Y=ev.year,
        m="%02d" % ev.month,
        d="%02d" % ev.day,
        H="%02d" % ev.hour,
        M="%02d" % ev.minute,
    )
    setup(base)
    smag = SuperMAG(
        base, dates, uid="shibaji7",
        stations=mag_stations
    )
    from utils import get_nearest_station, compute_Egeo_from_Bgeo
    datasets = {}
    for stn in mag_stations:
        o = (smag.sm_data[
            smag.sm_data.iaga==stn           
        ]).copy()
        glon, glat = (
            o.glon.tolist()[0],
            o.glat.tolist()[0]
        )
        dT = 60
        Ex, Ey, site, station = compute_Egeo_from_Bgeo(
            np.array(o.E_geo), np.array(o.N_geo), 
            glat=glat, glon=glon, dT=dT
        )
        o["mag_stn"], o["GC_km"] = stn, station["GC_km"]
        o["Ex"], o["Ey"]= Ex, Ey
        #print(o.E_geo.tolist()[:10], o.Ex.tolist()[:10])
        datasets[site.name] = o
    
    fig = plt.figure(figsize=(5, 7), dpi=240)
    for i, st in enumerate(datasets.keys()):
        print(st)
        data = datasets[st]
        ax = fig.add_subplot(311+i)
        ax.set_ylabel(r"$E_{geo}$ (mV/km)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        ax.xaxis.set_major_locator(mdates.HourLocator(byhour=range(0, 24, 1)))
        ax.set_xlim(dates)
        for comp, col, lab in zip(["x", "y"], ["r", "b"], [r"$E_x$", r"$E_y$"]):
            ax.plot(
                data.tval,
                data[f"E{comp}"],
                ls="-",
                color=col,
                lw=0.6,
                label=lab,
            )
        if i==0: 
            ax.legend(loc=1, fontsize=12)
            ax.text(0.1, 1.05, dates[1].strftime("%d, %b %Y"), ha="left", va="center", transform=ax.transAxes)
        ax.set_ylim(-50, 100)
        txt = f"({chr(97+i)}) {st} / {data.iaga.tolist()[0]}" #+ rf"$GC_{km},\phi_m$={'%.1f'%data.mlat.tolist()[0]}, {'%.1f'%data.mlon.tolist()[0]}"
                #"\n"+ fr"$\chi$={'%.2f'%np.mean(data.sza)}"
                #rf"$\theta_g,\phi_g$={data.glat.tolist()[0]}, {data.glon.tolist()[0]}" + \
        ax.axvline(ev, ls="-", lw=0.5, color="k")
        ax.text(0.1, 0.8, txt, ha="left", va="center", transform=ax.transAxes, fontdict={"size":8})
    ax.set_xlabel("Time (UT)")
    fig.savefig(base+"figures/Egeo.png", bbox_inches="tight")
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
            dt.datetime(2005,9,7,20)
        ]
    )
    mag_stations = [
        "BOU", "FRD", "FRN", 
        "MEA", "NEW", 
        "OTT", "PIN", "VIC",
        "BLC", "BRW", "BET",
        "CBB", "T16", "CMO",
        "CNL", "DAW", "DLR",
        "EAG", "EDM", "EWA",
        "FCC", "FMC", "FSP",
        "SMI", "FYU", "GAK",
        "GIM", "HOM", "HON",
        "ISL", "KUV", "T40",
        "PGC", "T21", "PKR",
        "PBQ", "T37", "RAL",
        "RES", "SHU", "SIT",
        "TAL", "TEO", "TUC",
        "T25", "VIC", "T38"
    ]
    # mag_stations = [
    #     "RES", "RAL", "PIN",
    #     "FRD", "FRN", "TUC"
    # ]
    # run_create_stack_plots(
    #     ev, dates, mag_stations
    # )
    run_mag_e_stack_plots(
        ev, dates, ["OTT", "FRD", "PIN"]
    )