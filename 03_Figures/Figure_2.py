#!/usr/bin/env python
# ==========================================================
# Plot 3-component waveforms for two OBS stations
# Waveforms + tide overlay, with current speed and direction above
#
# Input SAC filename example:
# ZS.T01.LC.BHE.M.2025.073.000000.SAC
#
# Notes:
#   1. SAC data are acceleration records after response removal.
#   2. Min-max downsampling is applied only for plotting.
#   3. The original waveform amplitudes are not filtered or normalized.
#
# Author : Yang Li
# Date   : 2026-07
# ==========================================================

import gc
import glob
import os
from pathlib import Path
from datetime import datetime, timedelta

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from obspy import UTCDateTime, read
from matplotlib.ticker import FuncFormatter

# Global style
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 7,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "savefig.dpi": 600,
})


# Repository-relative paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = REPO_ROOT.parent
DATA_ROOT = PROJECT_ROOT / "ESS_OBS_Shield_Data_v1.1.0"

DATA_DIR = str(DATA_ROOT / "Pro_Seismic_acceleration_v1.0")
CURRENT_DIR = str(DATA_ROOT / "Raw_Current_Tide_v1.1")
TIDE_FILE = os.path.join(CURRENT_DIR, "Tide_UTC.txt")
SPEED_FILE = os.path.join(CURRENT_DIR, "Flow_Speed_UTC.txt")
DIRECTION_FILE = os.path.join(CURRENT_DIR, "Flow_Direction_UTC.txt")
CURRENT_SAMPLE_SECONDS = 10.0

OUT_DIR = str(REPO_ROOT / "outputs" / "Figure_2")
RMS_FILE = str(DATA_ROOT / "Processed_RMS_v1.0" / "Results_RMS.mat")


# Data configuration
NETWORK = "ZS"
LOCATION = "LC"
STATIONS = ["T01", "T02"]
COMPONENTS = ["BHZ", "BHE", "BHN"]

#  2025-03-14 至 2025-03-31
PLOT_START = datetime(2025, 3, 14, 0, 0, 0)
PLOT_END = datetime(2025, 4, 1, 0, 0, 0)

# Figure configuration
FIG_W = 7.120
FIG_H_TIDE = 1.6
DPI = 600

LABEL_FONTSIZE = 7.0
TICK_FONTSIZE = 6.5

MAX_PLOT_POINTS = 40000

# Min-max downsampling
def minmax_downsample_trace(trace, max_points=40000):
    """
    Downsample one ObsPy Trace for long-duration waveform plotting.
    The data are divided into consecutive blocks. For each block,
    the minimum and maximum samples are retained in their original
    temporal order. This preserves transient amplitudes better than
    simply retaining every nth sample.

    Parameters
    ----------
    trace : obspy.Trace
        Input waveform trace.
    max_points : int
        Approximate maximum number of output points.

    Returns
    -------
    plot_time : numpy.ndarray
        Matplotlib serial date numbers.
    plot_data : numpy.ndarray
        Downsampled waveform values.
    """

    data = trace.data

    if np.ma.isMaskedArray(data):
        data = data.filled(np.nan)

    data = np.asarray(data)
    npts = data.size

    if npts == 0:
        return np.array([]), np.array([])

    # No downsampling needed
    if npts <= max_points:
        selected_indices = np.arange(npts, dtype=np.int64)

    # Min-max downsampling
    else:
        target_blocks = max(1, max_points // 2)
        block_size = int(np.ceil(npts / target_blocks))
        n_complete = (npts // block_size) * block_size
        index_parts = []

        if n_complete > 0:
            blocks = data[:n_complete].reshape(-1, block_size)

            local_min = np.argmin(blocks, axis=1)
            local_max = np.argmax(blocks, axis=1)

            block_start = (np.arange(blocks.shape[0], dtype=np.int64) * block_size )

            index_min = block_start + local_min
            index_max = block_start + local_max

            # Retain the min/max samples in chronological order.
            block_indices = np.column_stack(
                (index_min, index_max)
            )
            block_indices.sort(axis=1)

            index_parts.append(block_indices.ravel())

        # Process the final incomplete block.
        if n_complete < npts:
            remainder = data[n_complete:]

            if remainder.size > 0:
                remainder_indices = np.array([
                    n_complete + int(np.argmin(remainder)),
                    n_complete + int(np.argmax(remainder)),
                ], dtype=np.int64)

                remainder_indices.sort()
                index_parts.append(remainder_indices)

        selected_indices = np.concatenate(index_parts)

        # Explicitly retain the first and last samples.
        selected_indices = np.unique(
            np.concatenate((
                np.array([0], dtype=np.int64),
                selected_indices,
                np.array([npts - 1], dtype=np.int64),
            ))
        )

    # Convert only the selected sample indices to plotting times.
    start_num = mdates.date2num(trace.stats.starttime.datetime)
    delta_day = float(trace.stats.delta) / 86400.0

    plot_time = start_num + selected_indices * delta_day
    plot_data = data[selected_indices]

    return plot_time, plot_data


# Read waveform data
def read_waveforms(
    data_dir,
    stations=STATIONS,
    components=COMPONENTS,
    network=NETWORK,
    location=LOCATION,
    start_plot=PLOT_START,
    end_plot=PLOT_END,
    max_plot_points=MAX_PLOT_POINTS,
):
    """
    Read and downsample the requested OBS waveform records.

    Expected filename:
        NETWORK.STATION.LOCATION.COMPONENT.M.YEAR.JDAY.HHMMSS.SAC

    Example:
        ZS.T01.LC.BHE.M.2025.073.000000.SAC
    """

    data_all = {
        component: {}
        for component in components
    }

    utc_start = UTCDateTime(start_plot)

    # Use a nearly exclusive end time.
    utc_end = UTCDateTime(end_plot) - 1.0e-6

    for component in components:
        for station in stations:

            pattern = os.path.join(
                data_dir,
                (
                    f"{network}.{station}.{location}."
                    f"{component}.M.*.SAC"
                )
            )

            files = sorted(glob.glob(pattern))

            print("-" * 70)
            print(f"Station   : {station}")
            print(f"Component : {component}")
            print(f"Pattern   : {pattern}")
            print(f"Files     : {len(files)}")

            if not files:
                print(
                    f"WARNING: No SAC file found for "
                    f"{station}.{location}.{component}"
                )
                continue

            station_times = []
            station_data = []

            for sac_file in files:
                # Read header first and check temporal overlap.
                try:
                    header_stream = read(
                        sac_file,
                        format="SAC",
                        headonly=True,
                    )
                    header_trace = header_stream[0]

                except Exception as exc:
                    print(f"WARNING: Cannot read header: {sac_file}")
                    print(f"         {exc}")
                    continue

                file_start = header_trace.stats.starttime
                file_end = header_trace.stats.endtime

                if file_end < utc_start or file_start > utc_end:
                    print(
                        "Skip outside plotting period: "
                        f"{os.path.basename(sac_file)}"
                    )
                    del header_stream
                    continue

                print(f"Reading: {os.path.basename(sac_file)}")
                print(
                    f"  File range: "
                    f"{file_start.datetime} -- "
                    f"{file_end.datetime}"
                )

                del header_stream
                try:
                    stream = read(
                        sac_file,
                        format="SAC",
                        starttime=utc_start,
                        endtime=utc_end,
                    )

                except Exception as exc:
                    print(f"WARNING: Cannot read waveform: {sac_file}")
                    print(f"         {exc}")
                    continue

                for trace in stream:
                    # Apply trimming again to ensure exact limits.
                    trace.trim(
                        starttime=utc_start,
                        endtime=utc_end,
                        nearest_sample=False,
                        pad=False,
                    )

                    if trace.stats.npts == 0:
                        continue

                    original_npts = trace.stats.npts

                    plot_time, plot_data = minmax_downsample_trace(
                        trace,
                        max_points=max_plot_points,
                    )

                    if plot_data.size == 0:
                        continue

                    station_times.append(plot_time)
                    station_data.append(plot_data)

                    print(
                        f"  Samples: {original_npts:,} "
                        f"-> plotting points: {plot_data.size:,}"
                    )
                    print(
                        f"  Sampling interval: "
                        f"{trace.stats.delta:.6f} s"
                    )

                # Release the large waveform array before
                # reading the next SAC file.
                del stream
                gc.collect()

            if not station_times:
                print(
                    f"WARNING: No usable data for "
                    f"{station}.{location}.{component}"
                )
                continue

            combined_time = np.concatenate(station_times)
            combined_data = np.concatenate(station_data)

            # Ensure chronological order if several SAC files exist.
            sort_index = np.argsort(combined_time)
            combined_time = combined_time[sort_index]
            combined_data = combined_data[sort_index]

            data_all[component][station] = {
                "time": combined_time,
                "data": combined_data,
            }

            print(
                f"Loaded {station}.{location}.{component}: "
                f"{combined_data.size:,} plotting points"
            )

    return data_all


# Read tide data
def read_tide(
    tide_file,
    start_plot=PLOT_START,
    end_plot=PLOT_END,
):
    """
    Read tide level data.

    Expected columns:
        YYYY-MM-DD HH:MM:SS tide_level_cm

    Tide levels are converted from cm to m.
    """

    tide_time = []
    tide_level = []

    with open(tide_file, "r", encoding="utf-8") as file_object:
        # Skip header.
        next(file_object, None)

        for line_number, line in enumerate(
            file_object,
            start=2,
        ):
            parts = line.strip().split()

            if len(parts) < 3:
                continue

            time_string = " ".join(parts[:2])

            try:
                current_time = datetime.strptime(
                    time_string,
                    "%Y-%m-%d %H:%M:%S",
                )
                current_level = float(parts[2]) / 100.0

            except ValueError:
                print(
                    f"WARNING: Invalid tide record "
                    f"at line {line_number}: {line.strip()}"
                )
                continue

            tide_time.append(current_time)
            tide_level.append(current_level)

    tide_time = np.asarray(tide_time, dtype=object)
    tide_level = np.asarray(tide_level, dtype=float)

    # Half-open interval: [start_plot, end_plot)
    mask = np.array([
        start_plot <= current_time < end_plot
        for current_time in tide_time
    ])

    tide_time = tide_time[mask]
    tide_level = tide_level[mask]

    print("-" * 70)
    print(f"Tide records loaded: {tide_level.size:,}")

    if tide_level.size > 0:
        print(
            f"Tide range: "
            f"{tide_time[0]} -- {tide_time[-1]}"
        )
        print(
            f"Tide level: "
            f"{np.nanmin(tide_level):.3f} -- "
            f"{np.nanmax(tide_level):.3f} m"
        )

    return tide_time, tide_level


# Read environmental observations without interpolation or averaging.
def read_current(file_path, is_direction=False, start_plot=PLOT_START, end_plot=PLOT_END):
    records = {}
    conflicts = set()
    with open(file_path, "r", encoding="utf-8-sig") as handle:
        next(handle, None)
        for line_number, line in enumerate(handle, start=2):
            parts = line.split()
            if not parts:
                continue
            try:
                time = datetime.strptime(" ".join(parts[:2]), "%Y-%m-%d %H:%M:%S")
                value = float(parts[2])
            except (ValueError, IndexError):
                print(f"WARNING: Invalid current record: {file_path}:{line_number}")
                continue
            if not start_plot <= time < end_plot:
                continue
            if not np.isfinite(value) or value < 0 or (is_direction and value > 360):
                value = np.nan
            elif is_direction:
                value %= 360.0
            if time in conflicts:
                continue
            previous = records.get(time, np.nan)
            if np.isfinite(previous) and np.isfinite(value):
                if not np.isclose(previous, value, rtol=0, atol=1e-10):
                    conflicts.add(time)
                    records[time] = np.nan
            elif time not in records or np.isfinite(value):
                records[time] = value
    times = sorted(records)
    values = np.array([records[t] for t in times], dtype=float)
    if not np.any(np.isfinite(values)):
        raise ValueError(f"No valid current observations in plotting interval: {file_path}")
    print(f"Current records: {file_path}: {len(times)}; conflicting timestamps: {len(conflicts)}")
    return mdates.date2num(times), values


def break_current_gaps(times, values):
    gaps = np.flatnonzero(np.diff(times) * 86400 > 1.5 * CURRENT_SAMPLE_SECONDS) + 1
    return np.insert(times, gaps, times[gaps]), np.insert(values, gaps, np.nan)


# RMS utilities
def matlab_datenum_to_datetime(values):
    values = np.asarray(values, dtype=float).ravel()
    result = []
    for value in values:
        day = int(value)
        frac = value - day
        result.append(
            datetime.fromordinal(day)
            + timedelta(days=frac)
            - timedelta(days=366)
        )
    return np.asarray(result, dtype=object)


def load_rms_results(mat_file):
    """Read Results_RMS.mat generated by the MATLAB RMS script."""
    if not os.path.isfile(mat_file):
        raise FileNotFoundError(f"RMS MAT file not found: {mat_file}")

    mat = loadmat(mat_file, simplify_cells=True)
    if "Results" not in mat:
        raise KeyError(f"'Results' is absent from {mat_file}")

    raw = mat["Results"]
    if isinstance(raw, dict):
        raw = [raw]
    elif isinstance(raw, np.ndarray):
        raw = raw.tolist()

    rms = {}
    for item in raw:
        if not isinstance(item, dict):
            raise TypeError("Unexpected Results structure in RMS MAT file.")

        station = str(item["station"]).strip()
        time_dt = matlab_datenum_to_datetime(item["TIME"])
        rms[station] = {
            "time": time_dt,
            "time_num": mdates.date2num(time_dt.tolist()),
            "curr": np.asarray(item["CURR"], dtype=float).ravel(),
            "rms_h": np.asarray(item["BH_RMS"], dtype=float).ravel(),
            "rms_z": np.asarray(item["BZ_RMS"], dtype=float).ravel(),
        }

    for station in STATIONS:
        if station not in rms:
            raise KeyError(f"RMS data for {station} are absent from {mat_file}")

    if "plot_start_date" in mat:
        rms_start = matlab_datenum_to_datetime([float(np.asarray(mat["plot_start_date"]).squeeze())])[0]
    else:
        rms_start = datetime(2025, 3, 15, 0, 0, 0)

    if "plot_end_date" in mat:
        rms_end = matlab_datenum_to_datetime([float(np.asarray(mat["plot_end_date"]).squeeze())])[0]
    else:
        rms_end = datetime(2025, 3, 30, 23, 0, 0)

    print("-" * 70)
    print(f"RMS file: {mat_file}")
    print(f"RMS display range: {rms_start} -- {rms_end}")
    for station in STATIONS:
        print(f"{station}: {rms[station]['time'].size} RMS windows")

    return rms, rms_start, rms_end


def get_common_rms_ylim(rms, start_dt, end_dt):
    values = []
    for station in STATIONS:
        time = rms[station]["time"]
        mask = np.array([(start_dt <= t <= end_dt) for t in time])
        for key in ("rms_h", "rms_z"):
            y = rms[station][key][mask]
            y = y[np.isfinite(y) & (y > 0)]
            if y.size:
                values.append(y)

    if not values:
        raise ValueError("No positive finite RMS values in the selected display interval.")

    y = np.concatenate(values)
    lower_exp = np.floor(2 * np.log10(np.min(y))) / 2
    upper_exp = np.ceil(2 * np.log10(np.max(y))) / 2
    if upper_exp <= lower_exp:
        upper_exp = lower_exp + 0.5
    return 10.0 ** np.array([lower_exp, upper_exp])


def style_rms_axis(ax, y_lim):
    ax.set_yscale("log")
    ax.set_ylim(y_lim)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE, length=2.2, pad=0.5, direction="in")
    ax.tick_params(axis="y", right=True)
    ax.grid(True, which="major", linestyle="--", linewidth=0.35, alpha=0.35)
    ax.grid(False, which="minor")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)


# Plot Figure 2
def plot_figure2(
    data_all,
    tide_time,
    tide_level,
    speed_time,
    speed_values,
    direction_time,
    direction_values,
    rms,
    rms_start,
    rms_end,
    fig_w=FIG_W,
    fig_h=FIG_H_TIDE,
    start_plot=PLOT_START,
    end_plot=PLOT_END,
):
    stations = STATIONS
    components = COMPONENTS
    colors = {"T01": "k", "T02": "r"}

    # Preserve the height of the existing waveform/environment block, then
    # append two RMS rows using the original Figure 2 arrangement.
    original_row_height = fig_h * (0.99 - 0.05) / (len(components) + 0.10 * (len(components) - 1))
    row_gap = original_row_height * 0.10
    # top_height = fig_h + 3 * original_row_height + 2 * row_gap
    top_height = (fig_h + 3 * original_row_height + 2 * row_gap) * 0.85
    rms_row_height = 1.55
    total_height = top_height + 2 * rms_row_height + 0.34

    fig = plt.figure(figsize=(fig_w, total_height), dpi=DPI, facecolor="white")
    outer = fig.add_gridspec(
        2, 1,
        height_ratios=[top_height, 2 * rms_row_height + 0.20],
        hspace=0.1,
    )

    top_grid = outer[0].subgridspec(
        5, 1,
        height_ratios=[1.3, 1.3, 1, 1, 1],
        hspace=0.10,
    )

    rms_grid = outer[1].subgridspec(
        2, 1,
        height_ratios=[1, 1],
        hspace=0.1,
    )
    mid_grid = rms_grid[0].subgridspec(1, 2, wspace=0.135)
    bottom_grid = rms_grid[1].subgridspec(1, 2, wspace=0.135)

    # Top block: current speed, direction and three-component waveforms.
    all_axes = []
    speed_axis = fig.add_subplot(top_grid[0, 0])
    all_axes.append(speed_axis)
    direction_axis = fig.add_subplot(top_grid[1, 0], sharex=speed_axis)
    all_axes.append(direction_axis)
    axes = []
    for i in range(3):
        ax = fig.add_subplot(top_grid[i + 2, 0], sharex=speed_axis)
        axes.append(ax)
        all_axes.append(ax)

    plot_time, plot_speed = break_current_gaps(speed_time, speed_values)
    speed_axis.plot(plot_time, plot_speed, color="#3B5B92", linewidth=0.45)
    speed_axis.set_ylabel("Speed (m/s)", fontsize=LABEL_FONTSIZE, labelpad=5)
    speed_axis.set_ylim(0, max(0.1, np.nanmax(speed_values) * 1.05))
    speed_axis.yaxis.set_major_locator(plt.MaxNLocator(nbins=3))

    direction_axis.scatter(
        direction_time, direction_values, s=0.5, color="#A65E2E",
        linewidths=0, rasterized=True,
    )
    direction_axis.set_ylabel("Direction (°)", fontsize=LABEL_FONTSIZE, labelpad=5)
    direction_axis.set_ylim(0, 360)
    direction_axis.set_yticks([0, 120, 300])

    for env_axis in (speed_axis, direction_axis):
        env_axis.grid(True, axis="x", linestyle="--", linewidth=0.35, alpha=0.3)
        env_axis.tick_params(axis="both", labelsize=TICK_FONTSIZE, length=2.2, pad=0.8, direction="in")
        env_axis.tick_params(axis="y", right=True)
        for spine in env_axis.spines.values():
            spine.set_linewidth(0.8)

    tide_axes = [axis.twinx() for axis in axes]
    component_labels = {"BHZ": "HHZ", "BHE": "HH1", "BHN": "HH2"}

    for index, component in enumerate(components):
        axis = axes[index]
        tide_axis = tide_axes[index]
        available = [sta for sta in stations if sta in data_all.get(component, {})]

        if available:
            maxima = []
            for station in available:
                waveform = data_all[component][station]["data"]
                finite = waveform[np.isfinite(waveform)]
                if finite.size:
                    maxima.append(np.max(np.abs(finite)))
            ymax = max(maxima) if maxima else 1.0
            if not np.isfinite(ymax) or ymax <= 0:
                ymax = 1.0

            for station in available:
                axis.plot(
                    data_all[component][station]["time"],
                    data_all[component][station]["data"],
                    color=colors[station], linewidth=0.15, alpha=0.95,
                    solid_capstyle="butt",
                )
            axis.set_ylim(-ymax, ymax)
        else:
            axis.text(0.5, 0.5, f"No {component} data", transform=axis.transAxes,
                      ha="center", va="center", fontsize=TICK_FONTSIZE)

        axis.set_yticks([])
        axis.set_ylabel(component_labels.get(component, component), fontsize=LABEL_FONTSIZE,
                        rotation=90, labelpad=5, va="center")
        axis.grid(True, axis="x", linestyle="--", linewidth=0.35, alpha=0.3)
        axis.tick_params(axis="both", which="major", labelsize=TICK_FONTSIZE,
                         length=2.2, pad=0.8, direction="in")
        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.6)

        if tide_level.size > 0:
            tide_axis.plot(tide_time, tide_level, color="0.55", linewidth=0.75, alpha=0.9)
        tide_axis.set_ylim(0, 4.5)
        tide_axis.set_yticks([1, 2, 3, 4])
        tide_axis.set_yticklabels(["1", "2", "3", "4"], fontsize=6, color="0.35")
        tide_axis.tick_params(axis="y", labelsize=TICK_FONTSIZE, length=2.2,
                              colors="0.35", pad=0.8, direction="in")
        tide_axis.set_ylabel("Tide level (m)" if index == 1 else "", fontsize=LABEL_FONTSIZE,
                             color="0.35", labelpad=2)
        for spine in tide_axis.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)
        tide_axis.spines["right"].set_color("0.45")

    for axis in all_axes:
        axis.set_xlim(start_plot, end_plot)
        axis.tick_params(axis="x", which="major", length=2.2, direction="in", labelbottom=False)
        axis.tick_params(axis="x", which="minor", length=1.2, direction="in")

    axes[-1].xaxis.set_major_locator(mdates.DayLocator(interval=1))
    axes[-1].xaxis.set_minor_locator(mdates.HourLocator(byhour=[12]))

    def format_top_date(x, pos=None):
        date = mdates.num2date(x).date()
        days = (date - start_plot.date()).days
        if start_plot.date() <= date < end_plot.date() and days % 2 == 0:
            return date.strftime("%b.%d")
        return ""

    axes[-1].xaxis.set_major_formatter(FuncFormatter(format_top_date))
    axes[-1].tick_params(axis="x", labelbottom=True, labelsize=TICK_FONTSIZE,
                         direction="in", pad=1.0, length=2.2)

    for ax in (speed_axis, direction_axis):
        ax.yaxis.set_label_coords(-0.045, 0.5)

    # RMS panels: (b) and (c) time evolution; (d) and (e) versus current speed.
    ax_b = fig.add_subplot(mid_grid[0, 0])
    ax_c = fig.add_subplot(mid_grid[0, 1])
    ax_d = fig.add_subplot(bottom_grid[0, 0])
    ax_e = fig.add_subplot(bottom_grid[0, 1])

    y_lim = get_common_rms_ylim(rms, rms_start, rms_end)
    tide_num = mdates.date2num(tide_time.tolist()) if tide_time.size else np.array([])
    rms_start_num = mdates.date2num(rms_start)
    rms_end_num = mdates.date2num(rms_end)
    tide_mask = (tide_num >= rms_start_num) & (tide_num <= rms_end_num) if tide_num.size else np.array([], dtype=bool)

    # (b) RMS_H versus time.
    handles_b = []
    for station in stations:
        t = rms[station]["time_num"]
        y = rms[station]["rms_h"]
        mask = (t >= rms_start_num) & (t <= rms_end_num) & np.isfinite(y) & (y > 0)
        h = ax_b.scatter(t[mask], y[mask], s=6, c=colors[station], edgecolors="none", zorder=3)
        handles_b.append(h)
    style_rms_axis(ax_b, y_lim)
    ax_b.set_xlim(rms_start_num, rms_end_num)
    ax_b.set_ylabel(r"$\mathrm{RMS}_{\mathrm{H}}$ (m/s$^2$)", fontsize=LABEL_FONTSIZE)

    tide_b = ax_b.twinx()
    if tide_num.size:
        h_tide, = tide_b.plot(tide_num[tide_mask], tide_level[tide_mask], color="0.55", linewidth=0.7)
    else:
        h_tide, = tide_b.plot([], [], color="0.55", linewidth=0.7)
    tide_b.set_ylim(0, 4.2)
    tide_b.set_yticks([1, 2, 3, 4])
    tide_b.tick_params(axis="y",right=True, labelright=False, length=2.2, colors="0.35", direction="in")
    tide_b.spines["right"].set_color("0.45")

    legend = ax_b.legend(
        [handles_b[0], handles_b[1], h_tide], ["T01", "T02", "Tide level"],
        loc="upper center", ncol=3, frameon=True, fontsize=6.2,
        handlelength=1.0, handletextpad=0.3, columnspacing=0.6, borderpad=0.25,
    )
    legend.get_frame().set_linewidth(0.5)

    # (c) RMS_Z versus time.
    for station in stations:
        t = rms[station]["time_num"]
        y = rms[station]["rms_z"]
        mask = (t >= rms_start_num) & (t <= rms_end_num) & np.isfinite(y) & (y > 0)
        ax_c.scatter(t[mask], y[mask], s=6, c=colors[station], edgecolors="none", zorder=3)
    style_rms_axis(ax_c, y_lim)
    ax_c.set_xlim(rms_start_num, rms_end_num)
    ax_c.set_ylabel(r"$\mathrm{RMS}_{\mathrm{Z}}$ (m/s$^2$)", fontsize=LABEL_FONTSIZE)

    def format_rms_date(x, pos=None):
        date = mdates.num2date(x)
        if date.month == 3 and 16 <= date.day <= 30 and (date.day - 16) % 3 == 0:
            return date.strftime("Mar.%d")
        return ""

    for ax in (ax_b, ax_c):
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax.xaxis.set_major_formatter(FuncFormatter(format_rms_date))

    tide_c = ax_c.twinx()
    if tide_num.size:
        tide_c.plot(tide_num[tide_mask], tide_level[tide_mask], color="0.55", linewidth=0.7)
    tide_c.set_ylim(0, 4.2)
    tide_c.set_yticks([1, 2, 3, 4])
    tide_c.set_ylabel("Tide level (m)", fontsize=LABEL_FONTSIZE, color="0.35", labelpad=2)
    tide_c.tick_params(axis="y", labelsize=TICK_FONTSIZE, length=2.2, colors="0.35", direction="in")
    tide_c.spines["right"].set_color("0.45")

    # (d) RMS_H versus current speed.
    for station in stations:
        x = rms[station]["curr"]
        y = rms[station]["rms_h"]
        t = rms[station]["time_num"]
        mask = ((t >= rms_start_num) & (t <= rms_end_num) &
                np.isfinite(x) & np.isfinite(y) & (y > 0))
        ax_d.scatter(x[mask], y[mask], s=6, c=colors[station], edgecolors="none")
    style_rms_axis(ax_d, y_lim)
    ax_d.set_xlim(0, 0.8)
    ax_d.set_xlabel("Current speed (m/s)", fontsize=LABEL_FONTSIZE)
    ax_d.set_ylabel(r"$\mathrm{RMS}_{\mathrm{H}}$ (m/s$^2$)", fontsize=LABEL_FONTSIZE)

    # (e) RMS_Z versus current speed.
    for station in stations:
        x = rms[station]["curr"]
        y = rms[station]["rms_z"]
        t = rms[station]["time_num"]
        mask = ((t >= rms_start_num) & (t <= rms_end_num) &
                np.isfinite(x) & np.isfinite(y) & (y > 0))
        ax_e.scatter(x[mask], y[mask], s=6, c=colors[station], edgecolors="none")
    style_rms_axis(ax_e, y_lim)
    ax_e.set_xlim(0, 0.8)
    ax_e.set_xlabel("Current speed (m/s)", fontsize=LABEL_FONTSIZE)
    ax_e.set_ylabel(r"$\mathrm{RMS}_{\mathrm{Z}}$ (m/s$^2$)", fontsize=LABEL_FONTSIZE)

    def format_current_tick(x, pos=None):
        if any(np.isclose(x, value) for value in (0.1, 0.3, 0.5, 0.7)):
            return f"{x:.1f}"
        return ""

    for ax in (ax_d, ax_e):
        ax.set_xlim(0, 0.8)
        ax.set_xticks(np.arange(0, 0.81, 0.1))
        ax.xaxis.set_major_formatter(FuncFormatter(format_current_tick))

    fig.subplots_adjust(left=0.065, right=0.975, top=0.990, bottom=0.055)

    os.makedirs(OUT_DIR, exist_ok=True)
    output_tiff = os.path.join(OUT_DIR, "Figure_2.tiff")
    plt.savefig(
        output_tiff,
        format="tiff",
        dpi=DPI,
        bbox_inches="tight",
        pad_inches=0.02,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)
    print("-" * 70)
    print(f"Saved: {output_tiff}")


# Main
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    data_all = read_waveforms(
        data_dir=DATA_DIR,
        stations=STATIONS,
        components=COMPONENTS,
        network=NETWORK,
        location=LOCATION,
        start_plot=PLOT_START,
        end_plot=PLOT_END,
        max_plot_points=MAX_PLOT_POINTS,
    )

    tide_time, tide_level = read_tide(
        tide_file=TIDE_FILE,
        start_plot=PLOT_START,
        end_plot=PLOT_END,
    )

    speed_time, speed_values = read_current(SPEED_FILE)
    direction_time, direction_values = read_current(DIRECTION_FILE, is_direction=True)
    rms, rms_start, rms_end = load_rms_results(RMS_FILE)

    plot_figure2(
        data_all=data_all,
        tide_time=tide_time,
        tide_level=tide_level,
        speed_time=speed_time,
        speed_values=speed_values,
        direction_time=direction_time,
        direction_values=direction_values,
        rms=rms,
        rms_start=rms_start,
        rms_end=rms_end,
        fig_w=FIG_W,
        fig_h=FIG_H_TIDE,
        start_plot=PLOT_START,
        end_plot=PLOT_END,
    )


if __name__ == "__main__":
    main()

