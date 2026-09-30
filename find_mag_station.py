#!/usr/bin/env python
"""Rank SuperMAG stations by geomagnetic disturbance during an event.

Examples
--------
Peak disturbance near 12:00 UT::

    python find_mag_station.py \
        --event 2017-09-06T12:00 \
        --start 2017-09-06T11:00 \
        --end 2017-09-06T15:00 \
        --max-lat 65 --threshold 200 --metric peak

Peak-to-peak disturbance over the complete window::

    python find_mag_station.py \
        --event 2017-09-06T12:00 \
        --start 2017-09-06T11:00 \
        --end 2017-09-06T15:00 \
        --max-lat 65 --metric peak-to-peak
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
from typing import Optional

import numpy as np
import pandas as pd


PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))


def parse_datetime(value: str) -> dt.datetime:
    """Parse an ISO datetime and return a timezone-naive datetime."""
    return dt.datetime.fromisoformat(value.replace("Z", "+00:00")).replace(
        tzinfo=None
    )


def default_data_file(event: dt.datetime) -> str:
    event_id = event.strftime("%Y-%m-%d-%H-%M")
    return os.path.join(PROJECT_ROOT, "tmp", event_id, "supermag.csv")


def station_quality_mask(group: pd.DataFrame) -> bool:
    """Require finite geographic coordinates and horizontal components."""
    required = ["glat", "glon", "N_geo", "E_geo"]
    if any(column not in group.columns for column in required):
        return False
    values = group[required].to_numpy(dtype=float)
    return len(group) > 1 and np.isfinite(values).all()


def rank_stations(
    data: pd.DataFrame,
    event: dt.datetime,
    start: dt.datetime,
    end: dt.datetime,
    max_latitude: float = 65.0,
    min_latitude: Optional[float] = None,
    min_longitude: Optional[float] = None,
    max_longitude: Optional[float] = None,
    window_minutes: float = 30.0,
    baseline_samples: int = 20,
    metric: str = "peak",
) -> pd.DataFrame:
    """Return stations ranked by the requested geomagnetic metric."""
    data = data.copy()
    data["tval"] = pd.to_datetime(data["tval"])
    data = data[(data.tval >= start) & (data.tval <= end)]
    data = data.sort_values(["iaga", "tval"])

    event_time = pd.Timestamp(event)
    near_start = event_time - pd.Timedelta(minutes=window_minutes)
    near_end = event_time + pd.Timedelta(minutes=window_minutes)
    rows = []

    for station, group in data.groupby("iaga", sort=True):
        group = group.sort_values("tval").reset_index(drop=True)
        if not station_quality_mask(group):
            continue

        latitude = float(group.glat.iloc[0])
        longitude = (float(group.glon.iloc[0]) + 180.0) % 360.0 - 180.0
        if latitude >= max_latitude:
            continue
        if min_latitude is not None and latitude < min_latitude:
            continue
        if min_longitude is not None and longitude < min_longitude:
            continue
        if max_longitude is not None and longitude > max_longitude:
            continue

        horizontal = np.hypot(group["N_geo"], group["E_geo"])
        nbase = min(baseline_samples, len(group))
        baseline = float(horizontal.iloc[:nbase].mean())
        group["delta_h"] = horizontal - baseline

        if metric == "peak":
            sample = group[
                (group.tval >= near_start) & (group.tval <= near_end)
            ]
            if sample.empty:
                continue
            index = sample["delta_h"].abs().idxmax()
            value = float(abs(group.loc[index, "delta_h"]))
            signed_value = float(group.loc[index, "delta_h"])
            peak_time = group.loc[index, "tval"]
        elif metric == "peak-to-peak":
            if group.empty:
                continue
            index_min = group["delta_h"].idxmin()
            index_max = group["delta_h"].idxmax()
            value = float(group["delta_h"].max() - group["delta_h"].min())
            signed_value = float(group.loc[index_max, "delta_h"])
            peak_time = group.loc[index_max, "tval"]
            peak_time_min = group.loc[index_min, "tval"]
        elif metric == "component":
            sample = group[
                (group.tval >= near_start) & (group.tval <= near_end)
            ]
            if sample.empty:
                continue
            north_base = float(group["N_geo"].iloc[:nbase].mean())
            east_base = float(group["E_geo"].iloc[:nbase].mean())
            north_delta = sample["N_geo"] - north_base
            east_delta = sample["E_geo"] - east_base
            component_value = pd.concat(
                [north_delta.abs(), east_delta.abs()]
            )
            index = component_value.idxmax()
            value = float(component_value.loc[index])
            signed_value = float(sample.loc[index, "N_geo"] - north_base)
            peak_time = sample.loc[index, "tval"]
        else:
            raise ValueError(f"Unsupported metric: {metric}")

        result = {
            "station": station,
            "glat": latitude,
            "glon": longitude,
            "metric": metric,
            "value_nT": value,
            "signed_delta_h_nT": signed_value,
            "peak_time": peak_time,
            "baseline_h_nT": baseline,
        }
        if metric == "peak-to-peak":
            result["minimum_time"] = peak_time_min
        rows.append(result)

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(
        "value_nT", ascending=False
    ).reset_index(drop=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Rank SuperMAG stations by geomagnetic disturbance."
    )
    parser.add_argument("--event", required=True, help="Peak time, ISO format")
    parser.add_argument("--start", help="Start of event window, ISO format")
    parser.add_argument("--end", help="End of event window, ISO format")
    parser.add_argument("--data", help="SuperMAG CSV path")
    parser.add_argument("--max-lat", type=float, default=65.0)
    parser.add_argument(
        "--min-lat", type=float, default=None,
        help="Minimum geographic latitude in degrees",
    )
    parser.add_argument(
        "--min-lon", type=float, default=None,
        help="Minimum geographic longitude in degrees (-180 to 180)",
    )
    parser.add_argument(
        "--max-lon", type=float, default=None,
        help="Maximum geographic longitude in degrees (-180 to 180)",
    )
    parser.add_argument("--window-minutes", type=float, default=30.0)
    parser.add_argument("--baseline-samples", type=int, default=20)
    parser.add_argument(
        "--metric",
        choices=("peak", "peak-to-peak", "component"),
        default="peak",
    )
    parser.add_argument("--threshold", type=float, default=200.0)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--output", help="Optional ranked CSV output path")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    event = parse_datetime(args.event)
    start = parse_datetime(args.start) if args.start else event - dt.timedelta(hours=1)
    end = parse_datetime(args.end) if args.end else event + dt.timedelta(hours=3)
    data_file = args.data or default_data_file(event)

    if not os.path.exists(data_file):
        raise FileNotFoundError(f"SuperMAG CSV not found: {data_file}")

    data = pd.read_csv(data_file)
    ranked = rank_stations(
        data,
        event,
        start,
        end,
        max_latitude=args.max_lat,
        min_latitude=args.min_lat,
        min_longitude=args.min_lon,
        max_longitude=args.max_lon,
        window_minutes=args.window_minutes,
        baseline_samples=args.baseline_samples,
        metric=args.metric,
    )
    if ranked.empty:
        print("No usable stations matched the requested criteria.")
        return

    print(f"Data: {data_file}")
    print(f"Stations below latitude {args.max_lat:g}: {len(ranked)}")
    print(f"Metric: {args.metric}; threshold: {args.threshold:g} nT")
    print(f"Stations above threshold: {(ranked.value_nT > args.threshold).sum()}")
    print(ranked.head(args.top).to_string(index=False))

    if args.output:
        ranked.to_csv(args.output, index=False)
        print(f"Saved ranked stations: {args.output}")


if __name__ == "__main__":
    main()
