#!/usr/bin/env python3
"""
Fig2_Waveform.py

Plot the three-component continuous acceleration waveforms for T01 and T02
with tide-level overlays. This script generates the waveform panel used in
Figure 2.

Expected repository layout
--------------------------
Test_SUSTech/
├── ESS_OBS_Shield_Data_v1.0/
│   ├── Pro_Seismic_acceleration_v1.0/
│   └── Raw_Current_Tide_v1.0/
│       └── Tide_UTC.txt
│
└── ESS_OBS_Shield_Scripts_v1.0/
    └── Fig2_Waveform.py

Requirements
------------
Python 3
ObsPy
NumPy
Matplotlib
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from obspy import read


plt.rcParams.update(
    {
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
    }
)

FIG_W = 7.120
FIG_H_TIDE = 1.6
DPI = 600

LABEL_FONTSIZE = 8.0
TICK_FONTSIZE = 7.5

START_PLOT = datetime(2025, 3, 14)
END_PLOT = datetime(2025, 3, 31)


def repository_paths() -> tuple[Path, Path, Path]:
    """Return the waveform input, tide input, and figure output paths."""

    script_root = Path(__file__).resolve().parent
    project_root = script_root.parent

    data_root = project_root / "ESS_OBS_Shield_Data_v1.0"

    waveform_dir = data_root / "Pro_Seismic_acceleration_v1.0"
    tide_file = data_root / "Raw_Current_Tide_v1.0" / "Tide_UTC.txt"
    output_dir = script_root / "outputs" / "Figure2"

    if not waveform_dir.is_dir():
        raise FileNotFoundError(
            f"Continuous acceleration directory not found: {waveform_dir}"
        )

    if not tide_file.is_file():
        raise FileNotFoundError(f"Tide file not found: {tide_file}")

    output_dir.mkdir(parents=True, exist_ok=True)

    return waveform_dir, tide_file, output_dir


def find_waveform_file(
    data_dir: Path,
    station: str,
    component: str,
) -> Path:
    """Find the archived continuous SAC file for one station/component."""

    exact_pattern = f"ZS.{station}.LC.{component}*.SAC"
    files = sorted(data_dir.glob(exact_pattern))

    if not files:
        fallback_pattern = f"*{station}*{component}*.SAC"
        files = sorted(data_dir.glob(fallback_pattern))

    if not files:
        raise FileNotFoundError(
            f"No continuous SAC file found for {station} {component} "
            f"under {data_dir}"
        )

    if len(files) > 1:
        print(
            f"WARNING: multiple files found for {station} {component}; "
            f"using {files[0].name}"
        )

    return files[0]


def read_waveforms(
    data_dir: Path,
    stations: tuple[str, ...] = ("T01", "T02"),
    components: tuple[str, ...] = ("BHE", "BHN", "BHZ"),
    start_plot: datetime = START_PLOT,
    end_plot: datetime = END_PLOT,
) -> dict[str, dict[str, dict[str, np.ndarray]]]:
    """Read and trim the continuous acceleration waveforms."""

    data_all: dict[str, dict[str, dict[str, np.ndarray]]] = {}

    for component in components:
        data_all[component] = {}

        for station in stations:
            waveform_file = find_waveform_file(
                data_dir=data_dir,
                station=station,
                component=component,
            )

            stream = read(str(waveform_file))
            trace = stream[0].copy()

            trace.trim(
                starttime=trace.stats.starttime.__class__(start_plot),
                endtime=trace.stats.starttime.__class__(end_plot),
                nearest_sample=True,
                pad=False,
            )

            if trace.stats.npts == 0:
                raise ValueError(
                    f"No samples remain after trimming {waveform_file}"
                )

            data_all[component][station] = {
                "time": np.asarray(trace.times("matplotlib")),
                "data": np.asarray(trace.data),
            }

    return data_all


def read_tide(
    tide_file: Path,
    start_plot: datetime = START_PLOT,
    end_plot: datetime = END_PLOT,
) -> tuple[np.ndarray, np.ndarray]:
    """Read UTC tide level and convert centimetres to metres."""

    tide_time: list[datetime] = []
    tide_level: list[float] = []

    with tide_file.open("r", encoding="utf-8") as file_handle:
        next(file_handle)

        for line in file_handle:
            parts = line.strip().split()

            if len(parts) < 3:
                continue

            time_text = " ".join(parts[:2])
            level_text = parts[2]

            tide_time.append(
                datetime.strptime(
                    time_text,
                    "%Y-%m-%d %H:%M:%S",
                )
            )
            tide_level.append(float(level_text) / 100.0)

    tide_time_array = np.asarray(tide_time)
    tide_level_array = np.asarray(tide_level)

    mask = (
        (tide_time_array >= start_plot)
        & (tide_time_array <= end_plot)
    )

    return tide_time_array[mask], tide_level_array[mask]


def plot_waveforms_with_tide(
    data_all: dict[str, dict[str, dict[str, np.ndarray]]],
    tide_time: np.ndarray,
    tide_level: np.ndarray,
    output_dir: Path,
    fig_w: float = FIG_W,
    fig_h: float = FIG_H_TIDE,
) -> Path:
    """Generate the Figure 2 waveform panel."""

    stations = ("T01", "T02")
    colors = {"T01": "k", "T02": "r"}
    components = ("BHZ", "BHE", "BHN")

    figure, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(fig_w, fig_h),
        dpi=DPI,
        sharex=True,
    )

    figure.patch.set_facecolor("white")
    tide_axes = [axis.twinx() for axis in axes]

    for index, component in enumerate(components):
        axis = axes[index]
        tide_axis = tide_axes[index]

        available_stations = [
            station
            for station in stations
            if station in data_all[component]
        ]

        if not available_stations:
            continue

        ymax = max(
            np.max(np.abs(data_all[component][station]["data"]))
            for station in available_stations
        )

        if ymax == 0:
            ymax = 1.0

        for station in available_stations:
            axis.plot(
                data_all[component][station]["time"],
                data_all[component][station]["data"],
                color=colors[station],
                linewidth=0.15,
                alpha=0.95,
            )

        axis.set_ylim(-ymax, ymax)
        axis.set_yticks([])
        axis.set_ylabel(
            component,
            fontsize=LABEL_FONTSIZE,
            rotation=90,
            labelpad=5,
            va="center",
        )
        axis.grid(
            True,
            axis="x",
            linestyle="--",
            linewidth=0.35,
            alpha=0.3,
        )
        axis.tick_params(
            axis="both",
            which="major",
            labelsize=TICK_FONTSIZE,
            length=2.2,
            pad=0.8,
            direction="in",
        )

        for spine in axis.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.6)

        axis.spines["left"].set_position(("outward", 0))

        tide_axis.plot(
            tide_time,
            tide_level,
            color="0.55",
            linewidth=0.75,
            alpha=0.9,
        )
        tide_axis.set_ylim(0, 4.5)
        tide_axis.set_yticks([1, 2, 3, 4])
        tide_axis.set_yticklabels(
            ["1", "2", "3", "4"],
            fontsize=6,
            color="0.35",
        )
        tide_axis.tick_params(
            axis="y",
            labelsize=TICK_FONTSIZE,
            length=2.2,
            colors="0.35",
            pad=0.8,
            direction="in",
        )

        if index == 1:
            tide_axis.set_ylabel(
                "Tide (m)",
                fontsize=LABEL_FONTSIZE,
                color="0.35",
                labelpad=2,
            )
        else:
            tide_axis.set_ylabel("")

        for spine in tide_axis.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(0.8)

        tide_axis.spines["right"].set_color("0.45")
        tide_axis.spines["right"].set_position(("outward", 0))

    xmin = START_PLOT
    xmax = END_PLOT

    axes[-1].set_xlim([xmin, xmax])

    xticks = np.linspace(
        mdates.date2num(xmin),
        mdates.date2num(datetime(2025, 3, 30)),
        9,
    )

    axes[-1].set_xticks(xticks)
    axes[-1].xaxis.set_major_formatter(
        mdates.DateFormatter("%b.%d")
    )
    axes[-1].tick_params(
        axis="x",
        labelsize=TICK_FONTSIZE,
        direction="in",
        pad=1.0,
        length=2.2,
    )

    plt.subplots_adjust(
        left=0.01,
        right=0.99,
        top=0.99,
        bottom=0.05,
        hspace=0.1,
    )

    output_file = output_dir / "Figure2_Waveform_Tide.png"

    plt.savefig(
        output_file,
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(figure)

    print(f"Saved: {output_file}")
    return output_file


def main() -> None:
    waveform_dir, tide_file, output_dir = repository_paths()

    data_all = read_waveforms(waveform_dir)
    tide_time, tide_level = read_tide(tide_file)

    plot_waveforms_with_tide(
        data_all=data_all,
        tide_time=tide_time,
        tide_level=tide_level,
        output_dir=output_dir,
        fig_w=FIG_W,
        fig_h=FIG_H_TIDE,
    )


if __name__ == "__main__":
    main()
