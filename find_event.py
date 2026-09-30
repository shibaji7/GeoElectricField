import sys
import numpy as np
from datetime import datetime

sys.path.extend(["py", "py/fetch", "py/process", "py/plots"])

from event_study import CachedGOES
from get_mag_data import SuperMAG
from efield import _emtf_site_fields

timestamp = datetime(2005, 9, 7, 17, 40)

smag = SuperMAG(
    "tmp/2005-09-07-17-40/",
    (datetime(2005, 9, 7, 16), datetime(2005, 9, 7, 20)),
    uid="shibaji7",
)

fields = _emtf_site_fields(
    smag,
    ".data/EMTF",
    dT=60,
    bounds=(-125, -66.9, 24.4, 49.4),
)

for site, magnetometer_id, data in fields:
    index = (data["tval"] - timestamp).abs().idxmin()
    row = data.loc[index]

    if abs((row["tval"] - timestamp).total_seconds()) > 90:
        continue

    magnitude = np.hypot(
        float(row["induced_E_east"]),
        float(row["induced_E_north"]),
    )

    if magnitude > 1000 and site.rating>=4:
        print(
            f"EMTF site={site.name}, "
            f"magnetometer={magnetometer_id}, "
            f"E={magnitude:.1f} mV/km "
            f"({magnitude / 1000:.2f} V/km)"
            f"R: {site.rating}"
        )