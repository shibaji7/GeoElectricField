from obspy.clients.fdsn import Client
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

client = Client("IRIS")  # resolves to service.iris.edu / EarthScope DMC

net, sta = "EM", "PA878"

# Metadata: fine to request over the full deployment (response doesn't carry waveform bulk)
inv = client.get_stations(network=net, station=sta, level="response",
                           starttime=UTCDateTime("2005-09-02"),
                           endtime=UTCDateTime("2005-09-23"))

# Waveforms: request ONLY the event window (+ padding), not 21 days —
# large multi-week, 8 Hz, 4-channel requests can silently truncate.
t0 = UTCDateTime("2005-09-07T15:30:00")
t1 = UTCDateTime("2005-09-07T20:30:00")

st = client.get_waveforms(network=net, station=sta, location="*",
                           channel="MFE,MFN,MQE,MQN",
                           starttime=t0, endtime=t1)

print(f"Fetched {len(st)} traces:")
for tr in st:
    print(f"  {tr.id}: {tr.stats.starttime} -> {tr.stats.endtime}, "
          f"{tr.stats.npts} pts @ {tr.stats.sampling_rate} Hz")

if len(st) == 0:
    raise RuntimeError("No data returned for this window — check station "
                        "availability (e.g. client.get_stations with "
                        "level='channel') before plotting.")

st_mag = st.select(channel="MF?")
st_ele = st.select(channel="MQ?")

st_mag.remove_response(inventory=inv, output="DEF")
st_ele.remove_response(inventory=inv, output="DEF")

# Magnetic response is in Tesla -> convert to nT
for tr in st_mag:
    tr.data = tr.data * 1e9

for tr in st_ele:
    tr.data = tr.data * 1e6  # V/m -> mV/km
    tr.stats.channel += "_mVkm"

st_plot = st_mag + st_ele

fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

# Zoom window for the actual event: 16:00-20:00 UT
xlim_start = UTCDateTime("2005-09-07T16:00:00").matplotlib_date
xlim_end = UTCDateTime("2005-09-07T20:00:00").matplotlib_date

panels = [
    ("MFE", "H (East), nT", "tab:blue"),
    ("MFN", "H (North), nT", "tab:blue"),
    ("MQE_mVkm", "E (East), mV/km", "tab:red"),
    ("MQN_mVkm", "E (North), mV/km", "tab:red"),
]

any_plotted = False
for ax, (chan, label, color) in zip(axes, panels):
    sel = st_plot.select(channel=chan)
    if not sel:
        ax.text(0.5, 0.5, f"No data for {chan}", ha="center", va="center",
                transform=ax.transAxes, color="gray")
        ax.set_ylabel(label)
        continue
    tr = sel[0]
    ax.plot_date(tr.times("matplotlib"), tr.data, "-", color=color, lw=0.8, marker="")
    ax.set_ylabel(label)
    ax.grid(alpha=0.3)
    ax.set_xlim(xlim_start, xlim_end)
    any_plotted = True

    locator = mdates.AutoDateLocator()
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))

if not any_plotted:
    print("WARNING: nothing plotted — none of the requested channels "
          "had data in this window.")

axes[-1].set_xlabel("Time (UT), 2005-09-07")
fig.suptitle(f"{sta} horizontal H and E fields, 2005-09-07 16:00-20:00 UT (X24 flare)")
fig.autofmt_xdate()
plt.tight_layout()
plt.savefig(f"pub_fig/{sta}_HE_horizontal_event.png", dpi=200)
# plt.show()

print(inv[0][0][0].response)

print(inv[0][0][0].response.instrument_sensitivity)
# inv.select(network="EM", station="PA852", channel="MFE").plot_response(min_freq=1e-3)

for cha_ in inv[0][0]:
    if cha_.code.startswith("MQ"):
        resp = cha_.response
        print(f"--- {cha_.code} ---")
        print(resp)


import numpy as np

# Raw counts, unprocessed
st_raw = client.get_waveforms(network=net, station=sta, location="*",
                               channel="MQE,MQN", starttime=t0, endtime=t1)

for tr_raw in st_raw:
    chan = tr_raw.stats.channel
    cha_meta = inv.select(network=net, station=sta, channel=chan)[0][0][0]
    overall_sens = cha_meta.response.instrument_sensitivity.value  # counts per (V/m), per the stage-product check above

    # Manual deconvolution: counts / (stage1_gain * stage2_gain) = V/m directly
    e_manual_Vm = tr_raw.data / overall_sens

    # ObsPy's automatic deconvolution
    tr_auto = tr_raw.copy()
    tr_auto.remove_response(inventory=inv, output="DEF")

    ratio = np.nanmedian(tr_auto.data / e_manual_Vm)
    print(f"{chan}: median(auto/manual) ratio = {ratio:.4f}")