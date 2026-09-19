#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 5 + P/S-wave SNR analysis for the paired OBS experiment.

SNR analysis
------------
- Main SNR band: 0.01-0.1 Hz.
- Input: instrument-response-corrected velocity in m/s.
- Supplement: two frequency bands, station-separated three-component waveforms.
- Tables: nine SNR bands plus three representative bands, all using 60-s windows.
- P-wave SNR: HHZ.
- S-wave SNR: rotation-invariant total horizontal RMS from HH1 and HH2:
      RMS_H = sqrt(mean(HH1^2 + HH2^2))
- Common pre-P ambient-noise window for both P and S.
- T01/T02 use exactly the same windows.
- SNR_dB = 20 * log10(RMS_signal / RMS_noise)

Requires:
    numpy
    matplotlib
    obspy
"""

from pathlib import Path
import csv
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from obspy import read, UTCDateTime
from obspy.taup import TauPyModel


# Repository-relative paths
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = REPO_ROOT.parent
DATA_ROOT = PROJECT_ROOT / "ESS_OBS_Shield_Data_v1.1.0"

INDIR = DATA_ROOT / "Processed_Earthquake_velocity_v1.0" / "202503280620"
OUTDIR = REPO_ROOT / "outputs" / "Figure_5_S4"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_FIG5 = OUTDIR / "Figure_5.tiff"
OUT_SNR_CSV = OUTDIR / "SNR_bp.csv"
OUT_SNR_TXT = OUTDIR / "SNR_summary.txt"


# Event and geometry
ORIGIN = UTCDateTime("2025-03-28T06:20:52.709")
MW = 7.7
EVENT_DEPTH_KM = 10.0
DIST_DEG = 24.8398

# If MANUAL_P_SEC / MANUAL_S_SEC are not None, they override TauP predictions.
# Values are seconds after ORIGIN.
MANUAL_P_SEC = None
MANUAL_S_SEC = None
TAUP_MODEL = "iasp91"


# Plot settings
PLOT_START_SEC = 0.0
PLOT_END_SEC = 3600.0

FIG_FILTER = (0.01, 0.1)
SNR_FILTER = (0.01, 0.1)

FILTER_CORNERS = 4
TAPER_FRACTION = 0.01
TAPER_MAX_LENGTH_SEC = 20.0

# SNR windows, seconds relative to predicted/manual arrivals
NOISE_REL_P = (-120.0, -60.0)  # common noise window relative to P
P_SIGNAL_REL = (-5.0, 55.0)
S_SIGNAL_REL = (-5.0, 55.0)

# Figure style
FIG_WIDTH_MM = 180.0
FIG_HEIGHT_MM = 100.0
DPI = 600
FONT = "Liberation Serif"
FS = 7.0
FS_PANEL = 7.0
FS_TITLE = 7.0
LW_WAVE = 0.32
LW_AXIS = 0.60
LW_PHASE = 0.45
PHASE_COLOR = "0.65"

COLOR_T01 = "black"
COLOR_T02 = "red"


# Mapping requested for revised manuscript
CHANNELS = {
    "HHZ": "BHZ",
    "HH1": "BHE",
    "HH2": "BHN",
}

STATIONS = ("T01", "T02")


def mm_to_inch(mm):
    return mm / 25.4


def sac_path(station, sac_channel):
    return INDIR / f"202503280620.ZS.{station}.{sac_channel}.sac"


def load_trace(station, sac_channel):
    path = sac_path(station, sac_channel)
    if not path.exists():
        raise FileNotFoundError(f"Missing SAC file: {path}")

    st = read(str(path))
    if len(st) != 1:
        raise RuntimeError(f"Expected one trace in {path}, found {len(st)}")

    tr = st[0]
    tr.data = np.asarray(tr.data, dtype=np.float64)

    # SAC reference times can differ from the catalog origin by a few
    # milliseconds because of header precision / sample alignment.  Accept
    # offsets smaller than one sample; they are physically negligible here.
    dt = float(tr.stats.delta)
    start_offset = float(tr.stats.starttime - ORIGIN)
    end_shortfall = float((ORIGIN + PLOT_END_SEC) - tr.stats.endtime)

    if end_shortfall > dt:
        raise RuntimeError(
            f"{path.name} does not cover the requested plot interval.\n"
            f"Trace: {tr.stats.starttime} -- {tr.stats.endtime}\n"
            f"Need through: {ORIGIN + PLOT_END_SEC}"
        )

    if start_offset > dt:
        raise RuntimeError(
            f"{path.name} starts more than one sample after earthquake origin.\n"
            f"Trace start: {tr.stats.starttime}, origin: {ORIGIN}, "
            f"offset={start_offset:.6f} s, dt={dt:.6f} s"
        )
    elif start_offset > 0:
        print(
            f"  Note: {path.name} starts {start_offset:.6f} s after origin "
            f"(< one sample, dt={dt:.6f} s); accepted without time shifting."
        )

    return tr


def preprocess_unfiltered(tr):
    """Demean + linear detrend only; no amplitude normalization."""
    out = tr.copy()
    out.detrend("demean")
    out.detrend("linear")
    return out


def preprocess_filtered(tr, freqmin, freqmax):
    """
    Filter the full 7200-s trace first, then slice later.
    This minimizes edge effects in the 3600-s displayed interval.
    """
    out = tr.copy()
    out.detrend("demean")
    out.detrend("linear")
    out.taper(max_percentage=TAPER_FRACTION, max_length=TAPER_MAX_LENGTH_SEC, type="cosine")
    out.filter(
        "bandpass",
        freqmin=freqmin,
        freqmax=freqmax,
        corners=FILTER_CORNERS,
        zerophase=True,
    )
    return out


def extract_relative(tr, start_sec, end_sec):
    """
    Extract data using the trace's true sample times relative to ORIGIN.

    This avoids forcing the first sample to exactly 0 s when the SAC header
    starts a few milliseconds after the catalog origin.
    """
    t = tr.times(reftime=ORIGIN)
    x = np.asarray(tr.data, dtype=np.float64)

    # Half-sample tolerance keeps the nearest boundary sample without changing
    # the SAC timing.
    half_dt = 0.5 * float(tr.stats.delta)
    mask = (t >= start_sec - half_dt) & (t <= end_sec + half_dt)

    if not np.any(mask):
        raise RuntimeError(
            f"No samples found in requested window {start_sec:.3f}-{end_sec:.3f} s "
            f"for {tr.id}"
        )

    return t[mask], x[mask]


def _valid_sac_marker(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return False
    return np.isfinite(v) and abs(v + 12345.0) > 1e-3


def get_arrival_seconds(raw_traces):
    """
    Arrival priority:
    1) MANUAL_P_SEC / MANUAL_S_SEC, if both are set.
    2) SAC picks: A for P, T0 for S, converted to seconds relative to
       the SAC origin marker O as (A-O) and (T0-O).
    3) TauP iasp91 fallback.

    Using A-O/T0-O is preferable here because the SAC reference time is
    millisecond-quantized and O is not exactly zero.
    """
    if MANUAL_P_SEC is not None and MANUAL_S_SEC is not None:
        return float(MANUAL_P_SEC), float(MANUAL_S_SEC), "manual"

    sac_p = []
    sac_s = []

    for sta in STATIONS:
        for label in CHANNELS:
            tr = raw_traces[sta][label]
            sac = getattr(tr.stats, "sac", None)
            if sac is None:
                continue

            a = getattr(sac, "a", None)
            t0 = getattr(sac, "t0", None)
            o = getattr(sac, "o", 0.0)

            if not _valid_sac_marker(o):
                o = 0.0
            else:
                o = float(o)

            if _valid_sac_marker(a):
                sac_p.append(float(a) - o)
            if _valid_sac_marker(t0):
                sac_s.append(float(t0) - o)

    if sac_p and sac_s:
        p_sec = float(np.median(sac_p))
        s_sec = float(np.median(sac_s))

        # Warn if header picks are not consistent among channels/stations.
        p_spread = float(np.ptp(sac_p)) if len(sac_p) > 1 else 0.0
        s_spread = float(np.ptp(sac_s)) if len(sac_s) > 1 else 0.0
        if p_spread > 0.05 or s_spread > 0.05:
            print(
                "Warning: SAC arrival markers differ among traces: "
                f"P spread={p_spread:.3f} s, S spread={s_spread:.3f} s"
            )

        return p_sec, s_sec, "SAC A/T0 markers relative to O"

    # Fallback only if SAC markers are unavailable.
    model = TauPyModel(model=TAUP_MODEL)

    parr = model.get_travel_times(
        source_depth_in_km=EVENT_DEPTH_KM,
        distance_in_degree=DIST_DEG,
        phase_list=["P", "p", "Pn"],
    )
    sarr = model.get_travel_times(
        source_depth_in_km=EVENT_DEPTH_KM,
        distance_in_degree=DIST_DEG,
        phase_list=["S", "s", "Sn"],
    )

    if not parr:
        raise RuntimeError("TauP returned no P arrival.")
    if not sarr:
        raise RuntimeError("TauP returned no S arrival.")

    p_sec = min(a.time for a in parr)
    s_sec = min(a.time for a in sarr)
    return float(p_sec), float(s_sec), f"TauP {TAUP_MODEL} fallback"


def rms(x):
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return np.nan
    return float(np.sqrt(np.mean(x ** 2)))


def horizontal_rms(h1, h2):
    """
    Rotation-invariant horizontal total RMS:
        sqrt(mean(HH1^2 + HH2^2))
    """
    h1 = np.asarray(h1, dtype=float)
    h2 = np.asarray(h2, dtype=float)
    n = min(h1.size, h2.size)
    if n == 0:
        return np.nan
    return float(np.sqrt(np.mean(h1[:n] ** 2 + h2[:n] ** 2)))


def safe_snr_db(signal_rms, noise_rms):
    if (
        not np.isfinite(signal_rms)
        or not np.isfinite(noise_rms)
        or signal_rms <= 0.0
        or noise_rms <= 0.0
    ):
        return np.nan
    return float(20.0 * np.log10(signal_rms / noise_rms))


def window_data(tr, start_sec, end_sec):
    _, x = extract_relative(tr, start_sec, end_sec)
    return x


def compute_snr(filtered_snr, p_sec, s_sec, band=None):
    band = SNR_FILTER if band is None else band
    noise0 = p_sec + NOISE_REL_P[0]
    noise1 = p_sec + NOISE_REL_P[1]
    p0 = p_sec + P_SIGNAL_REL[0]
    p1 = p_sec + P_SIGNAL_REL[1]
    s0 = s_sec + S_SIGNAL_REL[0]
    s1 = s_sec + S_SIGNAL_REL[1]

    if noise0 < 0:
        raise RuntimeError(
            "The selected pre-P noise window starts before the earthquake origin. "
            "Adjust NOISE_REL_P."
        )

    results = []

    for sta in STATIONS:
        z = filtered_snr[sta]["HHZ"]
        h1 = filtered_snr[sta]["HH1"]
        h2 = filtered_snr[sta]["HH2"]

        # P wave on vertical component
        z_noise = window_data(z, noise0, noise1)
        z_signal = window_data(z, p0, p1)
        p_noise_rms = rms(z_noise)
        p_signal_rms = rms(z_signal)
        p_snr = safe_snr_db(p_signal_rms, p_noise_rms)

        results.append({
            "station": sta,
            "phase": "P",
            "component": "HHZ",
            "band_hz": f"{band[0]:g}-{band[1]:g}",
            "noise_start_s": noise0,
            "noise_end_s": noise1,
            "signal_start_s": p0,
            "signal_end_s": p1,
            "noise_rms": p_noise_rms,
            "signal_rms": p_signal_rms,
            "snr_db": p_snr,
        })

        # S wave on rotation-invariant total horizontal energy
        h1_noise = window_data(h1, noise0, noise1)
        h2_noise = window_data(h2, noise0, noise1)
        h1_signal = window_data(h1, s0, s1)
        h2_signal = window_data(h2, s0, s1)

        s_noise_rms = horizontal_rms(h1_noise, h2_noise)
        s_signal_rms = horizontal_rms(h1_signal, h2_signal)
        s_snr = safe_snr_db(s_signal_rms, s_noise_rms)

        results.append({
            "station": sta,
            "phase": "S",
            "component": "HH=sqrt(HH1^2+HH2^2)",
            "band_hz": f"{band[0]:g}-{band[1]:g}",
            "noise_start_s": noise0,
            "noise_end_s": noise1,
            "signal_start_s": s0,
            "signal_end_s": s1,
            "noise_rms": s_noise_rms,
            "signal_rms": s_signal_rms,
            "snr_db": s_snr,
        })

    return results


def save_snr_csv(results):
    fields = [
        "station", "phase", "component", "band_hz",
        "noise_start_s", "noise_end_s",
        "signal_start_s", "signal_end_s",
        "noise_rms", "signal_rms", "snr_db",
    ]

    with OUT_SNR_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def result_lookup(results, station, phase):
    for row in results:
        if row["station"] == station and row["phase"] == phase:
            return row
    raise KeyError((station, phase))


def save_snr_summary(results, p_sec, s_sec, arrival_source):
    p1 = result_lookup(results, "T01", "P")["snr_db"]
    p2 = result_lookup(results, "T02", "P")["snr_db"]
    s1 = result_lookup(results, "T01", "S")["snr_db"]
    s2 = result_lookup(results, "T02", "S")["snr_db"]

    lines = [
        "P/S-wave SNR summary",
        "",
        f"Event origin: {ORIGIN}",
        f"Arrival source: {arrival_source}",
        f"P arrival: {p_sec:.2f} s after origin ({ORIGIN + p_sec})",
        f"S arrival: {s_sec:.2f} s after origin ({ORIGIN + s_sec})",
        f"SNR band: {SNR_FILTER[0]:g}-{SNR_FILTER[1]:g} Hz",
        "",
        f"Noise window: P{NOISE_REL_P[0]:+.0f} to P{NOISE_REL_P[1]:+.0f} s",
        f"P signal window: P{P_SIGNAL_REL[0]:+.0f} to P{P_SIGNAL_REL[1]:+.0f} s",
        f"S signal window: S{S_SIGNAL_REL[0]:+.0f} to S{S_SIGNAL_REL[1]:+.0f} s",
        "",
        f"P SNR, T01: {p1:.3f} dB",
        f"P SNR, T02: {p2:.3f} dB",
        f"P Delta SNR (T02-T01): {p2-p1:.3f} dB",
        "",
        f"S SNR, T01: {s1:.3f} dB",
        f"S SNR, T02: {s2:.3f} dB",
        f"S Delta SNR (T02-T01): {s2-s1:.3f} dB",
        "",
        "P-wave SNR uses HHZ.",
        "S-wave SNR uses rotation-invariant total horizontal RMS:",
        "RMS_H = sqrt(mean(HH1^2 + HH2^2)).",
        "",
        "Important: T01 and T02 use exactly the same filters and time windows.",
    ]
    OUT_SNR_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def common_ylim(data_a, data_b, pad=1.08):
    m = max(
        float(np.nanmax(np.abs(data_a))),
        float(np.nanmax(np.abs(data_b))),
    )
    if not np.isfinite(m) or m == 0.0:
        m = 1.0
    return -pad * m, pad * m


def configure_wave_axis(ax, is_bottom, show_xticklabels, show_xlabel):
    ax.set_xlim(PLOT_START_SEC, PLOT_END_SEC)
    ax.xaxis.set_major_locator(MultipleLocator(600))

    # Keep only the left y axis; amplitudes are not normalized.
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_linewidth(LW_AXIS)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.set_yticks([])
    ax.tick_params(axis="y", left=False, right=False)

    if is_bottom:
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_linewidth(LW_AXIS)
        ax.tick_params(
            axis="x", which="major", direction="in",
            length=2.6, width=LW_AXIS,
            bottom=True, top=False, labelbottom=show_xticklabels,
            labelsize=FS, pad=1.5,
        )
        if show_xlabel:
            ax.set_xlabel("Time since origin (s)", fontsize=FS_TITLE, labelpad=2.5)
        else:
            ax.set_xlabel("")
    else:
        ax.spines["bottom"].set_visible(False)
        ax.tick_params(
            axis="x", which="both",
            bottom=False, top=False, labelbottom=False,
        )


def draw_phase_lines(ax, p_sec, s_sec, add_labels=False):
    for x in (p_sec, s_sec):
        ax.axvline(
            x, color=PHASE_COLOR, lw=LW_PHASE,
            ls="-", dashes=(3, 2), zorder=1
        )

    if add_labels:
        ymin, ymax = ax.get_ylim()
        y = ymax - 0.08 * (ymax - ymin)
        ax.text(
            p_sec + 18, y, "P",
            ha="left", va="top",
            fontsize=FS_TITLE, fontfamily=FONT
        )
        ax.text(
            s_sec + 18, y, "S",
            ha="left", va="top",
            fontsize=FS_TITLE, fontfamily=FONT
        )


def prepare_plot_data(raw_traces):
    raw_proc = {sta: {} for sta in STATIONS}
    fig_filt = {sta: {} for sta in STATIONS}
    snr_filt = {sta: {} for sta in STATIONS}

    for sta in STATIONS:
        for label in CHANNELS:
            tr = raw_traces[sta][label]
            raw_proc[sta][label] = preprocess_unfiltered(tr)
            fig_filt[sta][label] = preprocess_filtered(
                tr, FIG_FILTER[0], FIG_FILTER[1]
            )
            snr_filt[sta][label] = preprocess_filtered(
                tr, SNR_FILTER[0], SNR_FILTER[1]
            )

    return raw_proc, fig_filt, snr_filt


def collect_display_data(processed):
    out = {sta: {} for sta in STATIONS}
    for sta in STATIONS:
        for label in CHANNELS:
            t, x = extract_relative(
                processed[sta][label],
                PLOT_START_SEC,
                PLOT_END_SEC,
            )
            out[sta][label] = (t, x)
    return out


def make_figure5(raw_display, filt_display, p_sec, s_sec, row_labels=None, outfile=None):
    if row_labels is None:
        row_labels = ("Unfiltered", f"Band-pass {FIG_FILTER[0]:g}-{FIG_FILTER[1]:g} Hz")
    outfile = OUT_FIG5 if outfile is None else outfile
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": [FONT],
        "font.size": FS,
        "axes.linewidth": LW_AXIS,
        "xtick.major.width": LW_AXIS,
        "ytick.major.width": LW_AXIS,
        "savefig.dpi": DPI,
    })

    fig = plt.figure(
        figsize=(mm_to_inch(FIG_WIDTH_MM), mm_to_inch(FIG_HEIGHT_MM))
    )

    # Compact ESS double-column layout.
    # Reduced left margin and overall height; moderate gaps between panels.
    outer = fig.add_gridspec(
        2, 2,
        left=0.032, right=0.98,
        bottom=0.085, top=0.97,
        wspace=0.06, hspace=0.15,
    )

    panel_defs = [
        ("a", "T01", row_labels[0], raw_display, COLOR_T01, 0, 0),
        ("b", "T02", row_labels[0], raw_display, COLOR_T02, 0, 1),
        ("c", "T01", row_labels[1],
         filt_display, COLOR_T01, 1, 0),
        ("d", "T02", row_labels[1],
         filt_display, COLOR_T02, 1, 1),
    ]

    # Matched y limits between T01 and T02 for each component within each row.
    row_ylims = {}
    for row_name, dataset in (("raw", raw_display), ("filt", filt_display)):
        row_ylims[row_name] = {}
        for label in CHANNELS:
            x1 = dataset["T01"][label][1]
            x2 = dataset["T02"][label][1]
            row_ylims[row_name][label] = common_ylim(x1, x2)

    for panel_letter, sta, proc_text, dataset, color, r, c in panel_defs:
        # Larger internal spacing makes the individual waveform axes slightly
        # shorter and visually separates HHZ/HH1/HH2 more clearly.
        inner = outer[r, c].subgridspec(3, 1, hspace=0.16)
        axes = [fig.add_subplot(inner[i, 0]) for i in range(3)]

        row_key = "raw" if r == 0 else "filt"

        for i, label in enumerate(("HHZ", "HH1", "HH2")):
            ax = axes[i]
            t, x = dataset[sta][label]
            ax.plot(t, x, color=color, lw=LW_WAVE, rasterized=False)
            ax.set_ylim(row_ylims[row_key][label])

            # Only HH2 carries x ticks. The x-axis label itself is shown
            # only for the lower row, i.e. panels (c) and (d).
            configure_wave_axis(
                ax,
                is_bottom=(i == 2),
                show_xticklabels=(i == 2 and r == 1),
                show_xlabel=(i == 2 and r == 1),
            )

            draw_phase_lines(ax, p_sec, s_sec, add_labels=(i == 0))

            if c == 0:
                ax.set_ylabel(
                    label,
                    rotation=90,
                    fontsize=FS_TITLE,
                    labelpad=5,
                    va="center",
                )
            else:
                ax.set_ylabel("")

        # Put station/filter information inside the HHZ axis at upper right,
        # rather than above the panel, to avoid crowding.
        axes[0].text(
            0.985, 0.9,
            f"{proc_text}; {sta}",
            transform=axes[0].transAxes,
            ha="right", va="top",
            fontsize=FS_TITLE,
            fontfamily=FONT,
        )

        # Panel label outside upper-left.
        axes[0].text(
            -0.05, 1.01, f"({panel_letter})",
            transform=axes[0].transAxes,
            ha="left", va="bottom",
            fontsize=FS_PANEL,
            fontfamily=FONT,
            clip_on=False,
        )

    fig.savefig(
        outfile,
        dpi=DPI,
        format="tiff",
        pil_kwargs={"compression": "tiff_lzw"},
        bbox_inches=None,
        facecolor="white",
    )
    plt.close(fig)


# Supplementary waveforms and frequency sensitivity
SUPP_BANDS = ((0.05, 0.1), (0.1, 1.0))
SNR_BANDS = ((0.01, 0.1), (0.02, 0.1), (0.03, 0.1), (0.05, 0.1),
             (0.05, 0.2), (0.05, 0.5), (0.1, 0.2), (0.1, 0.5), (0.1, 1.0))
TABLE_BANDS = ((0.01, 0.1), (0.05, 0.1), (0.1, 1.0))
OUT_SUPP_FIG = OUTDIR / "Figure.S4.tiff"
OUT_ALL_CSV = OUTDIR / "SNR_all_9bands_60s.csv"
OUT_TABLE_CSV = OUTDIR / "SNR_representative_3bands_60s.csv"
OUT_TABLE_TXT = OUTDIR / "SNR_representative_3bands_60s.txt"


def validate_event_traces(raw_traces, p_sec, s_sec):
    required_start = p_sec + NOISE_REL_P[0]
    required_end = s_sec + S_SIGNAL_REL[1]
    for sta in STATIONS:
        h1, h2 = raw_traces[sta]["HH1"], raw_traces[sta]["HH2"]
        if (len(h1.data) != len(h2.data) or
                abs(float(h1.stats.starttime - h2.stats.starttime)) > 1e-6 or
                abs(float(h1.stats.delta - h2.stats.delta)) > 1e-9):
            raise ValueError(f"{sta}: HH1 and HH2 must have matching sample times.")
        for ch in CHANNELS:
            tr = raw_traces[sta][ch]
            t, x = extract_relative(tr, required_start, required_end)
            dt = float(tr.stats.delta)
            if (abs(t[0] - required_start) > dt or abs(t[-1] - required_end) > dt
                    or not np.all(np.isfinite(tr.data))):
                raise ValueError(f"{sta} {ch}: incomplete event interval or nonfinite samples.")


def make_supplement(raw_traces, p_sec, s_sec, arrival_source):
    validate_event_traces(raw_traces, p_sec, s_sec)
    all_results, table_rows, display = [], [], {}
    view = (PLOT_START_SEC, PLOT_END_SEC)
    for band in SNR_BANDS:
        filtered = {sta: {ch: preprocess_filtered(raw_traces[sta][ch], *band)
                          for ch in CHANNELS} for sta in STATIONS}
        results = compute_snr(filtered, p_sec, s_sec, band=band)
        if not all(np.isfinite(r["snr_db"]) for r in results):
            raise ValueError(f"Invalid SNR in {band}; inspect signal/background RMS.")
        for r in results:
            delta = (result_lookup(results, "T02", r["phase"])["snr_db"] -
                     result_lookup(results, "T01", r["phase"])["snr_db"])
            all_results.append(dict(r, delta_snr_db=delta, arrival_source=arrival_source))
        if band in TABLE_BANDS:
            row = {"band_hz": f"{band[0]:g}-{band[1]:g}"}
            for phase in ("P", "S"):
                v1 = result_lookup(results, "T01", phase)["snr_db"]
                v2 = result_lookup(results, "T02", phase)["snr_db"]
                row.update({f"{phase}_T01_dB": v1, f"{phase}_T02_dB": v2,
                            f"{phase}_delta_dB": v2-v1})
            table_rows.append(row)
        if band in SUPP_BANDS:
            display[band] = {sta: {ch: extract_relative(filtered[sta][ch], *view)
                                  for ch in CHANNELS} for sta in STATIONS}

    for path, rows in ((OUT_ALL_CSV, all_results), (OUT_TABLE_CSV, table_rows)):
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    lines = ["P/S SNR of instrument-response-corrected velocity records",
             f"Event origin: {ORIGIN}; Mw {MW}", f"Arrival source: {arrival_source}",
             f"P: {p_sec:.3f} s; S: {s_sec:.3f} s after origin",
             f"Common background: P{NOISE_REL_P[0]:+g} to P{NOISE_REL_P[1]:+g} s",
             f"P window: {P_SIGNAL_REL}; S window: {S_SIGNAL_REL} s relative to arrivals",
             "P uses HHZ; S uses sqrt(mean(HH1^2 + HH2^2)).",
             "SNR = 20 log10(RMS_phase_window / RMS_preP_background).",
             "Phase-window RMS includes both earthquake signal and background noise.",
             "Delta SNR = T02 - T01. All table values are dB.", "",
             "Band (Hz)     P:T01   P:T02   P:Delta   S:T01   S:T02   S:Delta"]
    for row in table_rows:
        vals = [row[k] for k in list(row)[1:]]
        lines.append(f"{row['band_hz']:<12}" + " ".join(f"{v:8.2f}" for v in vals))
    OUT_TABLE_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    labels = tuple(f"Band-pass {lo:g}-{hi:g} Hz" for lo, hi in SUPP_BANDS)
    make_figure5(display[SUPP_BANDS[0]], display[SUPP_BANDS[1]], p_sec, s_sec,
                 row_labels=labels, outfile=OUT_SUPP_FIG)
    for path in (OUT_SUPP_FIG, OUT_ALL_CSV, OUT_TABLE_CSV, OUT_TABLE_TXT):
        print(f"Saved: {path}")


def main():
    print("Loading SAC files...")
    raw_traces = {sta: {} for sta in STATIONS}

    for sta in STATIONS:
        for label, sac_channel in CHANNELS.items():
            tr = load_trace(sta, sac_channel)
            raw_traces[sta][label] = tr
            sac = getattr(tr.stats, "sac", None)
            if sac is not None:
                b = getattr(sac, "b", np.nan)
                o = getattr(sac, "o", np.nan)
                a = getattr(sac, "a", np.nan)
                t0 = getattr(sac, "t0", np.nan)
                marker_text = f", B={b:.6f}, O={o:.6f}, A(P)={a:.3f}, T0(S)={t0:.3f}"
            else:
                marker_text = ""

            print(
                f"{sta} {label} <- {sac_channel}: "
                f"Fs={tr.stats.sampling_rate:g} Hz, "
                f"start={tr.stats.starttime}, end={tr.stats.endtime}"
                f"{marker_text}"
            )

    # Validate sampling rates within each station and across stations.
    rates = [
        raw_traces[sta][label].stats.sampling_rate
        for sta in STATIONS
        for label in CHANNELS
    ]
    if not np.allclose(rates, rates[0]):
        raise RuntimeError(f"Sampling rates are inconsistent: {rates}")

    p_sec, s_sec, arrival_source = get_arrival_seconds(raw_traces)
    print("")
    print(f"Arrival source: {arrival_source}")
    print(f"P = {p_sec:.2f} s after origin -> {ORIGIN + p_sec}")
    print(f"S = {s_sec:.2f} s after origin -> {ORIGIN + s_sec}")

    if not (PLOT_START_SEC < p_sec < s_sec < PLOT_END_SEC):
        raise RuntimeError(
            "Predicted/manual P/S arrivals are outside the 0-3600 s plot "
            "window or are not ordered P < S."
        )

    print("")
    validate_event_traces(raw_traces, p_sec, s_sec)
    print("Preprocessing...")
    raw_proc, fig_filt, snr_filt = prepare_plot_data(raw_traces)
    raw_display = collect_display_data(raw_proc)
    filt_display = collect_display_data(fig_filt)

    print("Making revised Figure 5...")
    make_figure5(raw_display, filt_display, p_sec, s_sec)

    print("Computing P/S SNR...")
    results = compute_snr(snr_filt, p_sec, s_sec)
    save_snr_csv(results)
    save_snr_summary(results, p_sec, s_sec, arrival_source)
    make_supplement(raw_traces, p_sec, s_sec, arrival_source)

    print("")
    for phase in ("P", "S"):
        r1 = result_lookup(results, "T01", phase)
        r2 = result_lookup(results, "T02", phase)
        delta = r2["snr_db"] - r1["snr_db"]
        print(
            f"{phase}: T01={r1['snr_db']:.3f} dB, "
            f"T02={r2['snr_db']:.3f} dB, "
            f"Delta(T02-T01)={delta:+.3f} dB"
        )

    print("")
    print("Done.")
    print(f"Figure 5:   {OUT_FIG5}")
    print(f"SNR CSV:    {OUT_SNR_CSV}")
    print(f"SNR summary:{OUT_SNR_TXT}")


if __name__ == "__main__":
    main()

