import glob
import os

import bezpy
import pandas as pd
from loguru import logger


def load_MT_site(station, MT_data_folder=".data/EMTF/"):
    """
    Read MT site infomation from EMTF data file

    Parameters:
      station: name of a MT site
      data_folder: directory of the EMTF files
    """

    xml_file = glob.glob(
        MT_data_folder + "MT_TF_*" + station + "*/*" + station + "*.xml"
    )

    # Reads in an IRIS XML file and returns a Site object
    if len(xml_file) == 1:
        site = bezpy.mt.read_xml(xml_file[0])
        if site.name == station:
            return site
        else:
            logger.info(f"Site other than {station} was found: {site.name}")
    else:
        logger.info(f"No site or multiple sites found: {xml_file}/{station}")
    return


def download_save_EarthScope_GeoEB(stations, save_dir="tmp/geomag_electric/USArray/"):
    """
    Download and save MT geomagnetic and geoelectric field data using bezpy

    Parameters:
      stations: names of MT sites
      save_dir: directory of hdf files to be saved
    """
    os.makedirs(save_dir, exist_ok=True)
    for sta in stations:
        if not os.path.exists(f"{save_dir}{sta}.hdf"):
            site = load_MT_site(sta)
            if site == None:
                continue
            else:
                site.download_waveforms()
                site.save_waveforms(directory=save_dir)
                logger.info(f"Done with station {sta}")
        else:
            logger.info(f"File exists station {sta}.hdf")
    return


def load_USArray_one_station_multi_cha(
    station,
    start_time,
    end_time,
    channels=["BN", "BE", "EN", "EE"],
    data_dir="tmp/geomag_electric/USArray/",
):
    """
    Load USArray geoE and/or geoB data for one station

    Parameters:
      station: name of a USArray station
      start_time: start time of the data to be loaded
      end_time: end time of the data to be loaded
      channels: channels of the USArray data to be loaded
                'BN' is northward geoB; 'BE' is eastward geoB;
                'EN' is northward geoE; 'EE' is eastward geoE.
    """
    o = pd.DataFrame()
    if os.path.exists(f"{data_dir}{station}.hdf"):
        logger.info(f"Load file {station}.hdf")
        site = load_MT_site(station)

        site.load_waveforms(directory=data_dir)
        o = site.waveforms
        o.index.names = ["time"]
        if "time" not in channels:
            channels.append("time")
        o = o.reset_index()
        o = o[channels]
        o = o[(o.time >= start_time) & (o.time <= end_time)]
        o = o.reset_index(drop=True)
    else:
        logger.info(f"File does exists {station}.hdf")
    return o.copy()


def plot_TS_USGS_dataset(
    o,
    ax,
    xlim,
    sta="",
    ylim=[-150, 150],
    comps={
        "ele": {"ls": "-", "lw": 0.5},
        "mag": {"ls": "--", "lw": 0.5},
    },
    xlabel="UT",
    ylabels=[r"$E_{GEO}$, $mv/km$", r"$B_{GEO}$, $nT$"],
    loc=2,
    plot_mag=False,
):
    """
    Overlay station data into axes
    """
    import matplotlib.dates as mdates

    ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
    ax.plot(
        o.time,
        o.EE,
        ls=comps["ele"]["ls"],
        color="r",
        lw=comps["ele"]["lw"],
        label=r"$E_E^{obs}$",
    )
    ax.plot(
        o.time,
        o.EE_sim,
        ls="--",
        color="k",
        lw=comps["ele"]["lw"],
        label=r"$E_E^{sim}$",
    )
    ax.plot(
        o.time,
        o.EN,
        ls=comps["ele"]["ls"],
        color="b",
        lw=comps["ele"]["lw"],
        label=r"$E_N^{obs}$",
    )
    ax.plot(
        o.time,
        o.EN_sim,
        ls="--",
        color="green",
        lw=comps["ele"]["lw"],
        label=r"$E_N^{sim}$",
    )
    ax.set_xlim(xlim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabels[0])
    ax.set_ylim(ylim)
    ax.legend(loc=loc, fontsize=6)
    ax.text(0.9, 0.95, sta.upper(), ha="right", va="center", transform=ax.transAxes)

    if plot_mag:
        ax = ax.twinx()
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r"%H^{%M}"))
        ax.plot(o.time, o.BE, ls=comps["mag"]["ls"], color="r", lw=comps["mag"]["lw"])
        ax.plot(o.time, o.BN, ls=comps["mag"]["ls"], color="b", lw=comps["mag"]["lw"])
        ax.set_ylabel(ylabels[1])
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    return
