#!/usr/bin/env bash
# ============================================================================
# 01_Remove_Instrument_Response.sh
#
# Remove the instrument response from uncorrected continuous SAC waveforms
# and convert the seismic channels to acceleration.
#
# Processing retained from the original workflow:
#   rmean
#   rtr
#   taper
#   trans from polszeros ... to acc freq 0.001 0.005 40 45
#
# IMPORTANT:
#   The input must be UNCORRECTED continuous SAC data.
#   Do not use Pro_Seismic_acceleration_v1.0 as the input because those files
#   have already been corrected and converted to acceleration.
#
# Default response file:
#   ../ESS_OBS_Shield_Data_v1.0/
#   Raw_Metadata_v1.0/Instrument-response.SACPZ
#
# Default input:
#   ../ESS_OBS_Shield_Data_v1.0/Raw_Continuous_SAC_v1.0
#
# Default output:
#   ./working_data/Pro_Seismic_acceleration_v1.0
#
# Usage:
#   ./01_Preprocessing/01_Remove_Instrument_Response.sh
#
# Optional explicit paths:
#   ./01_Preprocessing/01_Remove_Instrument_Response.sh \
#       INPUT_SAC_DIR OUTPUT_SAC_DIR RESPONSE_SACPZ
# ============================================================================

set -euo pipefail

export SAC_DISPLAY_COPYRIGHT=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"

DEFAULT_INPUT="${PROJECT_ROOT}/ESS_OBS_Shield_Data_v1.0/Raw_Continuous_SAC_v1.0"
DEFAULT_OUTPUT="${REPO_ROOT}/working_data/Pro_Seismic_acceleration_v1.0"
DEFAULT_PZ="${PROJECT_ROOT}/ESS_OBS_Shield_Data_v1.0/Raw_Metadata_v1.0/Instrument-response.SACPZ"

INPUT_DIR="${1:-${DEFAULT_INPUT}}"
OUTPUT_DIR="${2:-${DEFAULT_OUTPUT}}"
PZ_FILE="${3:-${DEFAULT_PZ}}"

if ! command -v sac >/dev/null 2>&1; then
    echo "ERROR: SAC executable was not found in PATH." >&2
    exit 1
fi

if ! command -v saclst >/dev/null 2>&1; then
    echo "ERROR: saclst executable was not found in PATH." >&2
    exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: uncorrected continuous SAC directory not found:" >&2
    echo "       ${INPUT_DIR}" >&2
    echo >&2
    echo "The current archived Pro_Seismic_acceleration_v1.0 directory is" >&2
    echo "already response-corrected and must not be used as this input." >&2
    exit 1
fi

if [[ ! -f "${PZ_FILE}" ]]; then
    echo "ERROR: SACPZ response file not found: ${PZ_FILE}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

processed_count=0
skipped_count=0

while IFS= read -r -d '' input_file; do
    channel="$(saclst kcmpnm f "${input_file}" | awk '{print $2}')"

    case "${channel}" in
        BHE|BHN|BHZ|HHE|HHN|HHZ)
            ;;
        *)
            echo "SKIP: unsupported/non-seismic channel ${channel}: ${input_file}"
            skipped_count=$((skipped_count + 1))
            continue
            ;;
    esac

    relative_path="${input_file#${INPUT_DIR}/}"
    output_file="${OUTPUT_DIR}/${relative_path}"
    output_parent="$(dirname "${output_file}")"
    mkdir -p "${output_parent}"

    echo "PROCESS: ${input_file}"
    echo "OUTPUT : ${output_file}"
    echo "CHANNEL: ${channel}"

    sac >/dev/null << EOF
r ${input_file}
rmean
rtr
taper
trans from polszeros s ${PZ_FILE} to acc freq 0.001 0.005 40 45
w ${output_file}
q
EOF

    if [[ ! -s "${output_file}" ]]; then
        echo "ERROR: output SAC file was not created: ${output_file}" >&2
        exit 1
    fi

    processed_count=$((processed_count + 1))
done < <(find "${INPUT_DIR}" -type f \( -iname "*.SAC" -o -iname "*.sac" \) -print0 | sort -z)

if (( processed_count == 0 )); then
    echo "ERROR: no supported seismic SAC files were processed." >&2
    exit 1
fi

echo
echo "Instrument-response removal completed."
echo "Input directory : ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Response file   : ${PZ_FILE}"
echo "Processed files : ${processed_count}"
echo "Skipped files   : ${skipped_count}"
