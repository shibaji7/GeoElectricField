#!/usr/bin/env python

"""get_sd_data.py: module is dedicated to fetch fitacf data from files."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import bz2
import datetime as dt
import glob
import os
import shutil

import numpy as np
import pandas as pd
import pydarnio as pydarn
from loguru import logger


def dat2fitacf(dates, rad, base="data/{r}/"):
    """
    Convert dat type to fitacf file type
    """
    dirc = base.format(r=rad)
    os.makedirs(dirc, exist_ok=True)
    flist = []
    for d in dates:
        datfname = "/sd-data/%d/dat/%s/%d%02d%02d*" % (
            d.year,
            rad,
            d.year,
            d.month,
            d.day,
        )
        files = glob.glob(datfname)
        files.sort()
        for f in files:
            fndat = f.split("/")[-1]
            fnraw = fndat.replace("dat", "rawacf").replace(".bz2", "")
            fnfit = fnraw.replace("rawacf", "fitacf")
            os.system("cp %s %s/." % (f, dirc))
            os.system("bzip2 -d %s/%s" % (dirc, fndat))
            fndat = fndat.replace(".bz2", "")
            cmd2raw = "dattorawacf %s/%s > %s/%s" % (dirc, fndat, dirc, fnraw)
            os.system(cmd2raw)
            cmd2fit = "make_fit %s/%s > %s/%s" % (dirc, fnraw, dirc, fnfit)
            os.system(cmd2fit)
            os.system("rm %s/%s %s/%s*" % (dirc, fnraw, dirc, fndat))
            os.system("bzip2 %s/%s" % (dirc, fnfit))
            flist.append(fnfit)
    return flist


def to_normal_scan_id(d, key="scan"):
    """Convert to normal scan code"""
    sid = np.array(d[key])
    sid[sid != 0] = 1
    d[key] = sid
    return d


class Gate(object):
    """Class object to hold each range cell value"""

    def __init__(self, bm, i, params=["v", "w_l", "gflg", "p_l", "v_e"], gflg_type=-1):
        """
        initialize the parameters which will be stored
        bm: beam object
        i: index to store
        params: parameters to store
        """
        for p in params:
            if len(getattr(bm, p)) > i:
                setattr(self, p, getattr(bm, p)[i])
            else:
                setattr(self, p, np.nan)
        if gflg_type >= 0 and len(getattr(bm, "gsflg")[gflg_type]) > 0:
            setattr(self, "gflg", getattr(bm, "gsflg")[gflg_type][i])
        return


class Beam(object):
    """Class to hold one beam object"""

    def __init__(self, _rad=None):
        """
        initialize the instance
        """
        self.stid = _rad
        return

    def set(
        self,
        time,
        d,
        s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang"],
        v_params=["pwr0", "v", "w_l", "gflg", "p_l", "slist", "v_e"],
    ):
        """
        Set all parameters
        time: datetime of beam
        d: data dict for other parameters
        s_param: other scalar params
        v_params: other list params
        """
        self.time = time
        for p in s_params:
            if p in d.keys():
                if p == "scan" and d[p] != 0:
                    setattr(self, p, 1)
                else:
                    setattr(self, p, d[p])
            else:
                setattr(self, p, None)
        for p in v_params:
            if p in d.keys():
                setattr(self, p, d[p])
            else:
                setattr(self, p, [])
        self.gs_estimation()
        return

    def set_nc(self, time, d, i, s_params, v_params):
        """
        Set all parameters
        time: datetime of beam
        d: data dict for other parameters
        s_param: other scalar params
        v_params: other list params
        """
        self.time = time
        for p in s_params:
            if p in d.keys():
                setattr(self, p, d[p][i])
            else:
                setattr(self, p, None)
        for p in v_params:
            if p in d.keys():
                setattr(self, p, np.array(d[p])[i, :])
                if "slist" not in v_params and p == "v":
                    setattr(self, "slist", np.argwhere(~np.isnan(getattr(self, "v"))))
                setattr(self, p, getattr(self, p)[~np.isnan(getattr(self, p))])
            else:
                setattr(self, p, [])
        return

    def copy(self, bm):
        """
        Copy all parameters
        """
        for p in bm.__dict__.keys():
            setattr(self, p, getattr(bm, p))
        return

    def gs_estimation(self):
        """
        Estimate GS flag using different criterion
        Cases -
                0. Sundeen et al. |v| + w/3 < 30 m/s
                1. Blanchard et al. |v| + 0.4w < 60 m/s
                2. Blanchard et al. [2009] |v| - 0.139w + 0.00113w^2 < 33.1 m/s
        """
        self.gsflg = {}
        if len(self.v) > 0 and len(self.w_l) > 0:
            self.gsflg[0] = ((np.abs(self.v) + self.w_l / 3.0) < 30.0).astype(int)
        if len(self.v) > 0 and len(self.w_l) > 0:
            self.gsflg[1] = ((np.abs(self.v) + self.w_l * 0.4) < 60.0).astype(int)
        if len(self.v) > 0 and len(self.w_l) > 0:
            self.gsflg[2] = (
                (np.abs(self.v) - 0.139 * self.w_l + 0.00113 * self.w_l**2) < 33.1
            ).astype(int)
        # Modified defination by S. Chakraborty: {W-[50-(0.7*(V+5)**2)]} < 0
        self.gsflg[3] = (
            (np.array(self.w_l) - (50 - (0.7 * (np.array(self.v) + 5) ** 2)) < 0)
        ).astype(int)
        return

    @staticmethod
    def is_gs_estimation(v, w, kind=0):
        """
        Estimate GS flag using different criterion
        Cases -
                0. Sundeen et al. |v| + w/3 < 30 m/s
                1. Blanchard et al. |v| + 0.4w < 60 m/s
                2. Blanchard et al. [2009] |v| - 0.139w + 0.00113w^2 < 33.1 m/s
        """
        if kind == 0:
            gs = ((np.abs(v) + w / 3.0) < 30.0).astype(int)
        if kind == 1:
            gs = ((np.abs(v) + w * 0.4) < 60.0).astype(int)
        if kind == 2:
            gs = ((np.abs(v) - 0.139 * w + 0.00113 * w**2) < 33.1).astype(int)
        if kind == 3:
            gs = ((w_l - (50 - (0.7 * (v + 5) ** 2)) < 0)).astype(int)
        return gs


class Scan(object):
    """Class to hold one scan (multiple beams)"""

    def __init__(self, stime=None, etime=None, stype="normal"):
        """
        initialize the parameters which will be stored
        stime: start time of scan
        etime: end time of scan
        stype: scan type
        """
        self.stime = stime
        self.etime = etime
        self.stype = stype
        self.beams = []
        return

    def update_time(self, up=True):
        """
        Update stime and etime of the scan.
        up: Update average parameters if True
        """
        self.stime = self.beams[0].time
        self.etime = self.beams[-1].time
        if up:
            self._populate_avg_params()
        return

    def _populate_avg_params(self):
        """
        Polulate average parameetrs
        """
        f, nsky = [], []
        for b in self.beams:
            f.append(getattr(b, "tfreq"))
            nsky.append(getattr(b, "noise.sky"))
        self.f, self.nsky = np.mean(f), np.mean(nsky)
        return

    def _estimat_skills(
        self, v_params=["v", "w_l", "p_l", "slist"], s_params=["bmnum"], verbose=False
    ):
        """
        Only used on the median filtered scan data.
        Estimate skills of the median filterd data
        """
        self.skills, labels = {}, {"CONV": [], "KDE": []}
        _u = {key: [] for key in v_params + s_params}
        for b in self.beams:
            labels["CONV"].extend(getattr(b, "gflg_conv"))
            labels["KDE"].extend(getattr(b, "gflg_kde"))
            l = len(getattr(b, "slist"))
            for p in v_params:
                _u[p].extend(getattr(b, p))
            for p in s_params:
                _u[p].extend([getattr(b, p)] * l)
        for name in ["CONV", "KDE"]:
            self.skills[name] = np.nan
            try:
                self.skills[name] = Skills(
                    pd.DataFrame.from_records(_u).values,
                    np.array(labels[name]),
                    name,
                    verbose=verbose,
                )
            except:
                logger.error(
                    " System exception in 'get_sd_data' while extimating skills."
                )
        return


class FetchData(object):
    """Class to fetch data from fitacf files for one radar for atleast a day"""

    def __init__(self, rad, date_range, filetype="fitacf", files=None, verbose=True):
        """
        initialize the vars
        rad = radar code
        date_range = [ start_date, end_date ]
        files = List of files to load the data from
        e.x :   rad = "sas"
                date_range = [
                    datetime.datetime(2017,3,17),
                    datetime.datetime(2017,3,18),
                ]
        """
        self.rad = rad
        self.date_range = date_range
        self.filetype = filetype
        self.files = files
        self.verbose = verbose
        if (rad is not None) and (date_range is not None) and (len(date_range) == 2):
            self._create_files()
        return

    def _create_files(self):
        """
        Create file names from date and radar code
        """
        if self.files is None:
            self.files = []
        reg_ex = "/sd-data/{year}/{ftype}/{rad}/{date}.*.{ftype}.bz2"  # config.get("amgeo.sdloc","reg_ex")
        days = (self.date_range[1] - self.date_range[0]).days + 2
        ent = -1
        for d in range(-1, days):
            e = self.date_range[0] + dt.timedelta(days=d)
            fnames = glob.glob(
                reg_ex.format(
                    year=e.year,
                    ftype=self.filetype,
                    rad=self.rad,
                    date=e.strftime("%Y%m%d"),
                )
            )
            fnames.sort()
            for fname in fnames:
                tm = fname.split(".")[1]
                sc = fname.split(".")[2]
                dus = dt.datetime.strptime(
                    fname.split(".")[0].split("/")[-1] + tm + sc, "%Y%m%d%H%M%S"
                )
                due = dus + dt.timedelta(hours=2)
                if (ent == -1) and (dus <= self.date_range[0] <= due):
                    ent = 0
                if ent == 0:
                    self.files.append(fname)
                if (ent == 0) and (dus <= self.date_range[1] <= due):
                    ent = -1
        return

    def convert2fitacf(self, tmp):
        """Convert fit to fitacf"""
        if self.files is None:
            self.files = []
        reg_ex = "/sd-data/{year}/{ftype}/{rad}/{date}*.{ftype}.bz2"
        days = (self.date_range[1] - self.date_range[0]).days + 2
        ent = -1
        for d in range(-1, days):
            e = self.date_range[0] + dt.timedelta(days=d)
            fnames = glob.glob(
                reg_ex.format(
                    year=e.year,
                    ftype=self.filetype,
                    rad=self.rad,
                    date=e.strftime("%Y%m%d"),
                )
            )
            fnames.sort()
            for fname in fnames:
                tm = (fname.split(".")[0].split("/")[-1])[8:10]
                day = (fname.split(".")[0].split("/")[-1])[:8]
                dus = dt.datetime.strptime(day + tm, "%Y%m%d%H")
                due = dus + dt.timedelta(hours=2)
                if (ent == -1) and (dus <= self.date_range[0] <= due):
                    ent = 0
                if ent == 0:
                    self.files.append(fname)
                if (ent == 0) and (dus <= self.date_range[1] <= due):
                    ent = -1
        nfiles = []
        for f in self.files:
            fnew = tmp + f.split("/")[-1]
            shutil.copy(f, fnew)
            defname = fnew.replace(".bz2", "")
            with open(fnew, "rb") as src, open(defname, "wb") as dest:
                dest.write(bz2.decompress(src.read()))
            os.system("fittofitacf " + defname + " > " + defname + "acf")
            defname = defname + "acf"
            fnew = fnew.replace(".fit", ".fitacf")
            with open(defname, "rb") as src, open(fnew, "wb") as dest:
                dest.write(bz2.compress(src.read()))
            nfiles.append(fnew)
        self.files = nfiles
        return

    def _parse_data(self, data, s_params, v_params, by, scan_prop):
        """
        Parse data by data type
        data: list of data dict
        params: parameter list to fetch
        by: sort data by beam or scan
        scan_prop: provide scan properties if by='scan'
                        {"stype": type of scan, "dur": duration in min}
        """
        _b, _s = [], []
        if self.verbose:
            logger.info("Started converting to beam data.")
        for d in data:
            time = dt.datetime(
                d["time.yr"],
                d["time.mo"],
                d["time.dy"],
                d["time.hr"],
                d["time.mt"],
                d["time.sc"],
                d["time.us"],
            )
            if time >= self.date_range[0] and time <= self.date_range[1]:
                bm = Beam(self.rad)
                bm.set(time, d, s_params, v_params)
                _b.append(bm)
        if self.verbose:
            logger.info("Converted to beam data.")
        if by == "scan":
            if self.verbose:
                logger.info("Started converting to scan data.")
            scan, sc = 0, Scan(None, None, scan_prop["stype"])
            sc.beams.append(_b[0])
            for _ix, d in enumerate(_b[1:]):
                if d.scan == 1 and d.time != _b[_ix].time:
                    sc.update_time()
                    _s.append(sc)
                    sc = Scan(None, None, scan_prop["stype"])
                    sc.beams.append(d)
                else:
                    sc.beams.append(d)
            _s.append(sc)
            if self.verbose:
                logger.info("Converted to scan data.")
        return _b, _s

    def convert_to_pandas(
        self,
        beams,
        s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang", "time"],
        v_params=["v", "w_l", "gflg", "p_l", "slist"],
    ):
        """
        Convert the beam data into dataframe
        """
        _o = dict(zip(s_params + v_params, ([] for _ in s_params + v_params)))
        for b in beams:
            l = len(getattr(b, "slist"))
            for p in v_params:
                _o[p].extend(getattr(b, p))
            for p in s_params:
                _o[p].extend([getattr(b, p)] * l)
        L = len(_o["slist"])
        for p in s_params + v_params:
            if len(_o[p]) < L:
                l = len(_o[p])
                _o[p].extend([np.nan] * (L - l))
        return pd.DataFrame.from_records(_o)

    def scans_to_pandas(
        self,
        scans,
        s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang", "time"],
        v_params=["v", "w_l", "gflg", "p_l", "slist", "elv"],
        start_scnum=0,
    ):
        """
        Convert the scan data into dataframe
        """
        new_cols = ["scnum", "sbnum"]
        _o = dict(
            zip(
                s_params + v_params + new_cols,
                ([] for _ in s_params + v_params + new_cols),
            )
        )
        for idn, s in enumerate(scans):
            for idh, b in enumerate(s.beams):
                l = len(getattr(b, "slist"))
                for p in v_params:
                    _o[p].extend(getattr(b, p))
                for p in s_params:
                    _o[p].extend([getattr(b, p)] * l)
                _o["scnum"].extend([idn + start_scnum] * l)
                _o["sbnum"].extend([idh] * l)
            L = len(_o["slist"])
            for p in s_params + v_params + new_cols:
                if len(_o[p]) < L:
                    l = len(_o[p])
                    _o[p].extend([np.nan] * (L - l))
        return pd.DataFrame.from_records(_o)

    def fetch_data(
        self,
        s_params=[
            "bmnum",
            "noise.sky",
            "tfreq",
            "scan",
            "nrang",
            "intt.sc",
            "intt.us",
            "mppul",
            "nrang",
            "rsep",
            "cp",
            "frang",
            "smsep",
            "lagfr",
            "channel",
            "bmazm",
        ],
        v_params=["pwr0", "v", "w_l", "gflg", "p_l", "slist", "v_e", "elv"],
        by="scan",
        scan_prop={"dur": 1, "stype": "normal"},
    ):
        """
        Fetch data from file list and return the dataset
        params: parameter list to fetch
        by: sort data by beam or scan
        scan_prop: provide scan properties if by='scan'
                   {"stype": type of scan, "dur": duration in min}
        """
        data = []
        for f in self.files:
            with bz2.open(f) as fp:
                fs = fp.read()
            if self.verbose:
                logger.info(f"Read file -  {f}")
            reader = pydarn.SDarnRead(fs, True)
            records = reader.read_fitacf()
            data += records
        if by is not None:
            beams, scans = self._parse_data(data, s_params, v_params, by, scan_prop)
            return beams, scans
        else:
            return data

    @staticmethod
    def FetchSD(base, rad, dates):
        """
        Static method to call the SD
        inevntory functions, pre-process them
        and store to object.
        """
        sd = FetchData(rad, dates)
        sd_fname = base + "SD.{rad}.csv".format(rad=rad)
        if not os.path.exists(sd_fname):
            _, scans = sd.fetch_data()
            o = sd.scans_to_pandas(scans)
            o.to_csv(sd_fname, index=False, header=True, float_format="%g")
        else:
            logger.info(f"Local {sd_fname}")
            o = pd.read_csv(sd_fname, parse_dates=["time"])
        return o


if __name__ == "__main__":
    FetchData.FetchSD(
        "tmp/2015-03-11-16-20/",
        "bks",
        [dt.datetime(2015, 3, 11, 16), dt.datetime(2015, 3, 11, 16, 30)],
    )

    os.system("rm *.log")
    os.system("rm -rf __pycache__/")
