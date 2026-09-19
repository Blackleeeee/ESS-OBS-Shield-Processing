#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Final ESS-style BHZ PSD and persistent high-frequency peak diagnostics."""

from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Repository-relative paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = REPO_ROOT.parent
DATA_ROOT = PROJECT_ROOT / "ESS_OBS_Shield_Data_v1.1.0"

PSD_T01_FILE = DATA_ROOT / "Processed_PSD_v1.0" / "PSD_T01_BHZ_20250314_20250331.PSD"
PSD_T02_FILE = DATA_ROOT / "Processed_PSD_v1.0" / "PSD_T02_BHZ_20250314_20250331.PSD"
DIFF_FILE = DATA_ROOT / "Processed_PSDDiff_v1.0" / "DiffPSD_BHZ_20250314_20250331.PSD"
FLOW_FILE = DATA_ROOT / "Raw_Current_Tide_v1.1" / "Flow_Speed_UTC.txt"
TIDE_FILE = DATA_ROOT / "Raw_Current_Tide_v1.1" / "Tide_UTC.txt"
OUT_DIR = REPO_ROOT / "outputs" / "Figure_S5"

# Analysis window: complete UTC days 15–30 March 2025
ANALYSIS_START = pd.Timestamp("2025-03-15 00:00:00")
ANALYSIS_END = pd.Timestamp("2025-03-31 00:00:00")

# Persistent high-frequency peak
PEAK_FMIN = 10.0
PEAK_FMAX = 40.0
PEAK_TRACK_HALF_WIDTH_HZ = 8.0

# Colors
COLOR_T01 = "#0072BD"
COLOR_T02 = "#D95319"
COLOR_DIFF = "#D91A1A"
COLOR_TIDE = "0.50"

# Figure style
FIG_W_CM = 18.0
P1_H_CM = 5.2
P2_H_CM = 15.0
DPI = 600

FONT_NAME = "Times New Roman"
FS_TICK = 7.0
FS_LABEL = 7.5
FS_PANEL = 8.0
AX_LW = 0.70
TICK_W = 0.60

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": [FONT_NAME, "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": FS_TICK,
    "axes.labelsize": FS_LABEL,
    "xtick.labelsize": FS_TICK,
    "ytick.labelsize": FS_TICK,
    "legend.fontsize": FS_TICK,
    "axes.linewidth": AX_LW,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": TICK_W,
    "ytick.major.width": TICK_W,
    "xtick.minor.width": 0.45,
    "ytick.minor.width": 0.45,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.7,
    "ytick.minor.size": 1.7,
    "savefig.dpi": DPI,
})


def matlab_datenum_to_datetime(dn):
    day = int(dn)
    frac = float(dn) - day
    return datetime.fromordinal(day) + timedelta(days=frac) - timedelta(days=366)


def read_psd(path):
    arr = np.loadtxt(path)
    times = pd.DatetimeIndex([matlab_datenum_to_datetime(v) for v in arr[0, 1:]])
    freq = np.asarray(arr[1:, 0], dtype=float)
    psd = np.asarray(arr[1:, 1:], dtype=float)

    order = np.argsort(times.values)
    return times[order], freq, psd[:, order]


def subset_period(data):
    times, freq, matrix = data
    mask = (times >= ANALYSIS_START) & (times < ANALYSIS_END)
    return times[mask], freq, matrix[:, mask]


def check_psd_grids(reference, other, name_ref, name_other):
    tref, fref, _ = reference
    toth, foth, _ = other

    if len(tref) != len(toth) or not np.array_equal(tref.values, toth.values):
        raise ValueError(f"Time grids differ between {name_ref} and {name_other}.")
    if len(fref) != len(foth) or not np.allclose(fref, foth, equal_nan=True):
        raise ValueError(f"Frequency grids differ between {name_ref} and {name_other}.")


def read_hourly_flow(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]

    time_col = "Time" if "Time" in df.columns else df.columns[0]
    speed_col = "Speed_mps" if "Speed_mps" in df.columns else df.columns[1]

    df[time_col] = pd.to_datetime(df[time_col], format="%Y-%m-%d %H:%M:%S", errors="coerce")
    df[speed_col] = pd.to_numeric(df[speed_col], errors="coerce")
    df = df.dropna(subset=[time_col, speed_col]).sort_values(time_col)
    df = df.drop_duplicates(subset=time_col).set_index(time_col)

    return df[speed_col].resample("1h").mean()


def read_tide(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).strip() for c in df.columns]

    time_col = "Time" if "Time" in df.columns else df.columns[0]
    level_col = "waterlevel" if "waterlevel" in df.columns else df.columns[1]

    df[time_col] = pd.to_datetime(df[time_col], format="%Y-%m-%d %H:%M:%S", errors="coerce")
    df[level_col] = pd.to_numeric(df[level_col], errors="coerce")
    df = df.dropna(subset=[time_col, level_col]).sort_values(time_col)

    df[level_col] = df[level_col] / 100.0
    mask = (df[time_col] >= ANALYSIS_START) & (df[time_col] < ANALYSIS_END)

    return df.loc[mask, time_col].to_numpy(), df.loc[mask, level_col].to_numpy()


def extract_persistent_peak(freq, psd):
    band = (freq >= PEAK_FMIN) & (freq <= PEAK_FMAX)
    fb = freq[band]
    sb = psd[band, :]

    if len(fb) < 3:
        raise ValueError("Too few frequency bins in the 10–40 Hz peak-search band.")

    median_psd = np.nanmedian(sb, axis=1)
    ref_freq = float(fb[np.nanargmax(median_psd)])

    track = np.abs(fb - ref_freq) <= PEAK_TRACK_HALF_WIDTH_HZ
    peak_freq = np.full(sb.shape[1], np.nan)
    peak_psd = np.full(sb.shape[1], np.nan)

    for i in range(sb.shape[1]):
        y = sb[:, i]
        finite = np.isfinite(y)
        if finite.sum() < 3:
            continue

        local = np.zeros(len(y), dtype=bool)
        local[1:-1] = (
            finite[1:-1]
            & finite[:-2]
            & finite[2:]
            & (y[1:-1] >= y[:-2])
            & (y[1:-1] >= y[2:])
        )

        candidates = np.where(local & track)[0]
        if len(candidates):
            j = candidates[np.argmin(np.abs(fb[candidates] - ref_freq))]
        else:
            idx = np.where(track & finite)[0]
            if len(idx) == 0:
                continue
            j = idx[np.nanargmax(y[idx])]

        peak_freq[i] = fb[j]
        peak_psd[i] = y[j]

    return peak_freq, peak_psd, ref_freq


def set_left_ylabel(ax, text):
    ax.set_ylabel(text)
    ax.yaxis.set_label_coords(-0.060, 0.5)


def add_panel_label(fig, ax, label):
    pos = ax.get_position()
    fig.text(pos.x0 - 0.07, pos.y1, label, ha="left", va="bottom", fontsize=FS_PANEL)


def style_frequency_axis(ax):
    ax.set_xscale("log")
    ax.set_xlim(0.01, 40.0)
    ticks = [0.01, 0.05, 0.1, 0.5, 1, 5, 10, 20, 40]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0.01", "0.05", "0.1", "0.5", "1", "5", "10", "20", "40"])
    ax.grid(True, which="major", linestyle="--", linewidth=0.35, alpha=0.28)
    ax.tick_params(top=True, right=True)


def style_time_axis(ax):
    ax.set_xlim(ANALYSIS_START, ANALYSIS_END)
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b.%d"))
    ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
    ax.grid(True, which="major", linestyle="--", linewidth=0.35, alpha=0.28)
    ax.tick_params(top=True, right=True, labelbottom=True)


def style_speed_axis(ax):
    ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))
    ax.grid(True, which="major", linestyle="--", linewidth=0.35, alpha=0.28)
    ax.tick_params(top=True, right=True)


def plot_hourly_psd_single(data, color, ylabel, label, filename, ylim=None, zero_line=False):
    _, freq, matrix = data

    fig, ax = plt.subplots(
        figsize=(FIG_W_CM / 2.54, P1_H_CM / 2.54),
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.065, right=0.992, bottom=0.22, top=0.965)

    for i in range(matrix.shape[1]):
        ax.plot(freq, matrix[:, i], color=color, alpha=0.10, linewidth=0.35, rasterized=True)

    ax.plot(freq, np.nanmedian(matrix, axis=1), color=color, linewidth=1.20, label="Median")
    style_frequency_axis(ax)
    ax.set_xlabel("Frequency (Hz)")
    set_left_ylabel(ax, ylabel)
    ax.text(0.008, 0.965, label, transform=ax.transAxes, ha="left", va="top", fontsize=FS_LABEL)

    if zero_line:
        ax.axhline(0, color="0.35", linewidth=0.6, linestyle="--")
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.legend(loc="upper right", frameon=False)

    out = OUT_DIR / filename
    fig.savefig(out, dpi=DPI, format="tiff", pil_kwargs={"compression": "tiff_lzw"}, facecolor="white")
    plt.close(fig)
    return out


def plot_final_composite(
    times, freq, psd_t02,
    peak_freq, peak_psd, ref_freq,
    flow_speed_hourly,
    tide_time, tide_level,
):
    valid_peak = np.isfinite(peak_freq) & np.isfinite(peak_psd)
    valid_flow = valid_peak & np.isfinite(flow_speed_hourly)

    fig = plt.figure(
        figsize=(FIG_W_CM / 2.54, P2_H_CM / 2.54),
        facecolor="white",
    )

    # Compact layout; every panel shows its own x-axis tick labels.
    gs = fig.add_gridspec(
        4, 1,
        left=0.08,
        right=0.955,
        bottom=0.075,
        top=0.975,
        hspace=0.3,
        height_ratios=[1.3, 1.3, 1.3, 1.3],
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[2, 0])
    ax_d = fig.add_subplot(gs[3, 0])

    # (a) T02-BHZ hourly PSD curves
    for i in range(psd_t02.shape[1]):
        ax_a.plot(
            freq, psd_t02[:, i],
            color=COLOR_T02,
            linewidth=0.42,
            alpha=0.20,
            rasterized=True,
        )

    style_frequency_axis(ax_a)
    ax_a.axvline(20.0, color="0.35", linestyle="--", linewidth=0.70)
    ax_a.set_xlabel("Frequency (Hz)", labelpad=2)
    set_left_ylabel(ax_a, "PSD (dB)")
    add_panel_label(fig, ax_a, "(a)")

    # (b) Peak frequency vs time
    ax_b.scatter(
        times[valid_peak], peak_freq[valid_peak],
        s=15,
        color=COLOR_T02,
        edgecolors="none",
        alpha=0.82,
        rasterized=True,
    )
    ax_b.axhline(ref_freq, color="0.45", linewidth=0.65, linestyle="--")
    ax_b.set_ylim(PEAK_FMIN, PEAK_FMAX)
    ax_b.set_yticks(np.arange(10, 41, 5))
    style_time_axis(ax_b)
    # ax_b.set_xlabel("Time (UTC)", labelpad=2)
    set_left_ylabel(ax_b, "Peak frequency (Hz)")
    add_panel_label(fig, ax_b, "(b)")

    # (c) Peak PSD vs time + tide level
    ax_c_tide = ax_c.twinx()
    ax_c_tide.plot(
        tide_time, tide_level,
        color=COLOR_TIDE,
        linewidth=0.85,
        alpha=0.95,
        zorder=1,
    )
    ax_c_tide.set_ylim(0, 4.5)
    ax_c_tide.set_yticks(np.arange(0, 5, 1))
    ax_c_tide.set_ylabel("Tide level (m)", color=COLOR_TIDE, fontsize=FS_LABEL, labelpad=3)
    ax_c_tide.tick_params(
        axis="y", colors=COLOR_TIDE,
        labelsize=FS_TICK,
        width=TICK_W,
        length=3.0,
        direction="in",
    )
    ax_c_tide.spines["right"].set_color(COLOR_TIDE)
    ax_c_tide.spines["right"].set_linewidth(AX_LW)

    ax_c.scatter(
        times[valid_peak], peak_psd[valid_peak],
        s=15,
        color=COLOR_T02,
        edgecolors="none",
        alpha=0.82,
        zorder=3,
        rasterized=True,
    )
    style_time_axis(ax_c)
    # ax_c.set_xlabel("Time (UTC)", labelpad=2)
    set_left_ylabel(ax_c, "Peak PSD (dB)")
    add_panel_label(fig, ax_c, "(c)")

    # (d) Peak PSD vs current speed
    ax_d.scatter(
        flow_speed_hourly[valid_flow],
        peak_psd[valid_flow],
        s=16,
        color=COLOR_T02,
        edgecolors="none",
        alpha=0.82,
        rasterized=True,
    )
    style_speed_axis(ax_d)
    ax_d.set_xlabel(r"Current speed (m s$^{-1}$)", labelpad=2)
    set_left_ylabel(ax_d, "Peak PSD (dB)")
    add_panel_label(fig, ax_d, "(d)")

    out = OUT_DIR / "Figure_S5.tiff"
    fig.savefig(
        out,
        dpi=DPI,
        format="tiff",
        pil_kwargs={"compression": "tiff_lzw"},
        facecolor="white",
    )
    plt.close(fig)
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    data_t01_full = read_psd(PSD_T01_FILE)
    data_t02_full = read_psd(PSD_T02_FILE)
    data_diff_full = read_psd(DIFF_FILE)

    check_psd_grids(data_t01_full, data_t02_full, "T01 BHZ", "T02 BHZ")
    check_psd_grids(data_t01_full, data_diff_full, "T01 BHZ", "DiffPSD BHZ")

    data_t01 = subset_period(data_t01_full)
    data_t02 = subset_period(data_t02_full)
    data_diff = subset_period(data_diff_full)

    # P1: retain three independent hourly-spectra figures.
    p1a = plot_hourly_psd_single(
        data_t01, COLOR_T01, "PSD (dB)", "T01 BHZ",
        "P1a_T01_BHZ_hourly_PSD_20250315_20250330.tiff",
    )
    p1b = plot_hourly_psd_single(
        data_t02, COLOR_T02, "PSD (dB)", "T02 BHZ",
        "P1b_T02_BHZ_hourly_PSD_20250315_20250330.tiff",
    )
    p1c = plot_hourly_psd_single(
        data_diff, COLOR_DIFF, r"PSD$_{\mathrm{diff}}$ (dB)", r"BHZ PSD$_{\mathrm{diff}}$",
        "P1c_BHZ_hourly_PSDdiff_20250315_20250330.tiff",
        ylim=(-30, 30),
        zero_line=True,
    )

    # P2: persistent T02-BHZ high-frequency peak.
    times, freq, psd_t02 = data_t02
    peak_freq, peak_psd, ref_freq = extract_persistent_peak(freq, psd_t02)

    hourly_flow = read_hourly_flow(FLOW_FILE)
    psd_hours = times.floor("h")
    flow_speed_hourly = hourly_flow.reindex(psd_hours).to_numpy(dtype=float)

    tide_time, tide_level = read_tide(TIDE_FILE)

    result = pd.DataFrame({
        "time_UTC": times,
        "peak_frequency_Hz": peak_freq,
        "peak_PSD_dB": peak_psd,
        "current_speed_hourly_mean_mps": flow_speed_hourly,
    })
    out_csv = OUT_DIR / "T02_BHZ_10_40Hz_persistent_peak_20250315_20250330.csv"
    result.to_csv(out_csv, index=False, float_format="%.6f")

    p2 = plot_final_composite(
        times, freq, psd_t02,
        peak_freq, peak_psd, ref_freq,
        flow_speed_hourly,
        tide_time, tide_level,
    )

    print(f"Analysis interval: {ANALYSIS_START} to {ANALYSIS_END} (end excluded)")
    print(f"Peak-search band: {PEAK_FMIN:.1f}-{PEAK_FMAX:.1f} Hz")
    print(f"Reference persistent peak: {ref_freq:.3f} Hz")
    print(f"Median tracked peak frequency: {np.nanmedian(peak_freq):.3f} Hz")
    print(f"Valid hourly peaks: {np.isfinite(peak_freq).sum()} / {len(peak_freq)}")
    print(f"Valid peak/current pairs: {np.sum(np.isfinite(peak_psd) & np.isfinite(flow_speed_hourly))}")
    print(f"Saved: {p1a}")
    print(f"Saved: {p1b}")
    print(f"Saved: {p1c}")
    print(f"Saved: {p2}")
    print(f"Saved: {out_csv}")


if __name__ == "__main__":
    main()

