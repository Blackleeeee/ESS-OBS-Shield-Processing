#!/usr/bin/env python3
"""
02_Split_Acceleration_SAC_to_Daily.py

Split continuous response-corrected acceleration SAC files into complete
UTC calendar-day SAC files.

Default input:
    ESS_OBS_Shield_Data_v1.0/Pro_Seismic_acceleration_v1.0

Default output:
    ESS_OBS_Shield_Scripts_v1.0/working_data/SAC_Day

Processing:
    - discard an incomplete first UTC day;
    - write only complete 24-hour UTC days;
    - skip an incomplete final day;
    - preserve waveform values;
    - do not detrend, taper, filter, or resample.

Requirements:
    Python 3
    ObsPy
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from obspy import UTCDateTime, read
from obspy.core import Trace
from obspy.core.util.attribdict import AttribDict


SECONDS_PER_DAY = 86400


def parse_arguments() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    project_root = repo_root.parent

    default_input = (
        project_root
        / "ESS_OBS_Shield_Data_v1.0"
        / "Pro_Seismic_acceleration_v1.0"
    )
    default_output = repo_root / "working_data" / "SAC_Day"

    parser = argparse.ArgumentParser(
        description=(
            "Split continuous acceleration SAC files into complete UTC "
            "calendar-day SAC files."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input,
        help=f"Continuous acceleration SAC directory. Default: {default_input}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output,
        help=f"Daily SAC output directory. Default: {default_output}",
    )
    parser.add_argument(
        "--network",
        default="ZS",
        help="Fallback network code for output filenames. Default: ZS",
    )
    return parser.parse_args()


def first_complete_midnight(start_time: UTCDateTime) -> UTCDateTime:
    midnight = UTCDateTime(
        year=start_time.year,
        month=start_time.month,
        day=start_time.day,
    )
    if abs(start_time - midnight) < 1.0e-6:
        return midnight
    return midnight + SECONDS_PER_DAY


def update_sac_headers(trace: Trace, start_time: UTCDateTime) -> None:
    if not hasattr(trace.stats, "sac") or trace.stats.sac is None:
        trace.stats.sac = AttribDict()

    trace.stats.sac.nzyear = start_time.year
    trace.stats.sac.nzjday = start_time.julday
    trace.stats.sac.nzhour = 0
    trace.stats.sac.nzmin = 0
    trace.stats.sac.nzsec = 0
    trace.stats.sac.nzmsec = 0
    trace.stats.sac.b = 0.0
    trace.stats.sac.e = (trace.stats.npts - 1) * trace.stats.delta


def split_one_file(
    sac_file: Path,
    output_root: Path,
    fallback_network: str,
) -> int:
    try:
        stream = read(str(sac_file))
    except Exception as exc:
        print(f"ERROR: failed to read {sac_file}: {exc}", file=sys.stderr)
        return 0

    if not stream:
        print(f"ERROR: no trace found in {sac_file}", file=sys.stderr)
        return 0

    trace = stream[0]
    sampling_rate = float(trace.stats.sampling_rate)

    if sampling_rate <= 0:
        print(
            f"ERROR: invalid sampling rate in {sac_file}: {sampling_rate}",
            file=sys.stderr,
        )
        return 0

    samples_per_day_float = SECONDS_PER_DAY * sampling_rate
    samples_per_day = int(round(samples_per_day_float))

    if not np.isclose(
        samples_per_day_float,
        samples_per_day,
        rtol=0.0,
        atol=1.0e-6,
    ):
        print(
            f"ERROR: non-integer samples per day for {sac_file}. "
            f"Sampling rate: {sampling_rate}",
            file=sys.stderr,
        )
        return 0

    station = str(trace.stats.station).strip()
    channel = str(trace.stats.channel).strip()
    network = str(trace.stats.network).strip() or fallback_network

    if not station or not channel:
        print(
            f"ERROR: missing station/channel code in {sac_file}",
            file=sys.stderr,
        )
        return 0

    original_start = trace.stats.starttime
    current_day = first_complete_midnight(original_start)
    total_npts = int(trace.stats.npts)
    written = 0

    while True:
        start_index = int(
            round((current_day - original_start) * sampling_rate)
        )
        end_index_exclusive = start_index + samples_per_day

        if start_index < 0:
            current_day += SECONDS_PER_DAY
            continue

        if end_index_exclusive > total_npts:
            break

        daily_data = trace.data[start_index:end_index_exclusive].copy()

        if daily_data.size != samples_per_day:
            print(
                f"WARNING: incomplete day skipped at {current_day} "
                f"for {sac_file}",
                file=sys.stderr,
            )
            current_day += SECONDS_PER_DAY
            continue

        month_text = current_day.strftime("%Y%m")
        day_text = current_day.strftime("%Y%m%d")

        output_dir = (
            output_root
            / month_text
            / day_text
            / station
        )
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / (
            f"{network}.{station}.{channel}.{day_text}.SAC"
        )

        daily_trace = trace.copy()
        daily_trace.data = daily_data
        daily_trace.stats.starttime = current_day
        daily_trace.stats.npts = daily_data.size
        update_sac_headers(daily_trace, current_day)

        print(
            f"WRITE: {output_file} | "
            f"start={daily_trace.stats.starttime} | "
            f"npts={daily_trace.stats.npts}"
        )

        try:
            daily_trace.write(str(output_file), format="SAC")
        except Exception as exc:
            print(
                f"ERROR: failed to write {output_file}: {exc}",
                file=sys.stderr,
            )
            current_day += SECONDS_PER_DAY
            continue

        written += 1
        current_day += SECONDS_PER_DAY

    return written


def main() -> None:
    args = parse_arguments()

    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()

    if not input_dir.is_dir():
        print(
            f"FATAL: input directory not found: {input_dir}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    sac_files = sorted(
        path
        for path in input_dir.rglob("*")
        if path.is_file() and path.suffix.upper() == ".SAC"
    )

    if not sac_files:
        print(
            f"FATAL: no SAC files found under: {input_dir}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    total_written = 0

    for sac_file in sac_files:
        print("\n" + "=" * 72)
        print(f"PROCESS: {sac_file}")
        print("=" * 72)

        total_written += split_one_file(
            sac_file=sac_file,
            output_root=output_dir,
            fallback_network=args.network,
        )

    print("\nDaily splitting completed.")
    print(f"Input directory : {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Daily files written: {total_written}")


if __name__ == "__main__":
    main()
