#!/usr/bin/env bash
# ============================================================================
# 03_Split_Daily_SAC_to_Hourly.sh
#
# Split complete 24-hour acceleration SAC files into one-hour SAC files for
# input to MyPSD_main.m.
#
# Processing retained from the original workflow:
#   rmean
#   rtr
#   taper type cosine width 0.1
#
# Default input:
#   ./working_data/SAC_Day
#
# Default output:
#   ./working_data/Hourly_Seismic_acceleration_v1.0
#
# Output structure:
#   YYYYMM/YYYYMMDD/STATION/STATION.COMPONENT/
#
# Usage:
#   ./01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh
#
# Optional explicit paths:
#   ./01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh \
#       DAILY_SAC_DIR HOURLY_SAC_DIR
# ============================================================================

set -euo pipefail

export SAC_DISPLAY_COPYRIGHT=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DEFAULT_INPUT="${REPO_ROOT}/working_data/SAC_Day"
DEFAULT_OUTPUT="${REPO_ROOT}/working_data/Hourly_Seismic_acceleration_v1.0"

INPUT_DIR="${1:-${DEFAULT_INPUT}}"
OUTPUT_DIR="${2:-${DEFAULT_OUTPUT}}"

if ! command -v sac >/dev/null 2>&1; then
    echo "ERROR: SAC executable was not found in PATH." >&2
    exit 1
fi

if ! command -v saclst >/dev/null 2>&1; then
    echo "ERROR: saclst executable was not found in PATH." >&2
    exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: daily SAC directory not found: ${INPUT_DIR}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
shopt -s nullglob

daily_count=0
hourly_count=0

for month_dir in "${INPUT_DIR}"/20????; do
    [[ -d "${month_dir}" ]] || continue
    month="$(basename "${month_dir}")"

    for day_dir in "${month_dir}"/20??????; do
        [[ -d "${day_dir}" ]] || continue
        day="$(basename "${day_dir}")"

        for station_dir in "${day_dir}"/*; do
            [[ -d "${station_dir}" ]] || continue
            station="$(basename "${station_dir}")"

            daily_files=("${station_dir}"/*.SAC)

            for input_file in "${daily_files[@]}"; do
                [[ -f "${input_file}" ]] || continue

                channel="$(saclst kcmpnm f "${input_file}" | awk '{print $2}')"

                if [[ -z "${channel}" ]]; then
                    echo "WARNING: failed to read channel from ${input_file}" >&2
                    continue
                fi

                component_dir="${OUTPUT_DIR}/${month}/${day}/${station}/${station}.${channel}"
                mkdir -p "${component_dir}"

                echo
                echo "PROCESS: ${input_file}"
                echo "STATION: ${station}"
                echo "CHANNEL: ${channel}"

                for ((hour_index=0; hour_index<=23; hour_index++)); do
                    begin_seconds=$((hour_index * 3600))
                    end_seconds=$((begin_seconds + 3600))
                    printf -v hour_text "%02d" "${hour_index}"

                    output_file="${component_dir}/ZS.${station}.${channel}.${day}_${hour_text}.SAC"

                    echo "  HOUR ${hour_text}: ${begin_seconds}-${end_seconds} s"

                    sac >/dev/null << EOF
setbb cb ${begin_seconds}
setbb ce ${end_seconds}
cut %cb% %ce%
r ${input_file}
rmean
rtr
taper type cosine width 0.1
write ${output_file}
q
EOF

                    if [[ ! -s "${output_file}" ]]; then
                        echo "ERROR: hourly SAC not created: ${output_file}" >&2
                        exit 1
                    fi

                    hourly_count=$((hourly_count + 1))
                done

                daily_count=$((daily_count + 1))
            done
        done
    done
done

if (( daily_count == 0 )); then
    echo "ERROR: no daily SAC files were processed under ${INPUT_DIR}" >&2
    exit 1
fi

echo
echo "Hourly splitting completed."
echo "Input directory       : ${INPUT_DIR}"
echo "Output directory      : ${OUTPUT_DIR}"
echo "Daily files processed : ${daily_count}"
echo "Hourly files generated: ${hourly_count}"
