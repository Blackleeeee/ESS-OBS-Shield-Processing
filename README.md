# ESS OBS Shield Processing and Figure Scripts

## Overview

This repository contains the MATLAB, Python, and shell scripts used to process continuous ocean-bottom seismic records and reproduce the RMS, PSD, PSD-difference, and current-speed analyses used to evaluate a double-spherical-shell shield for shallow-water ocean-bottom seismometers.

The archived analysis workflow begins with continuous seismic waveforms that have already been corrected for instrument response and converted to acceleration.

## Repository structure

```text
ESS_OBS_Shield_Scripts_v1.0/
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── .gitignore
│
├── 01_Preprocessing/
│   ├── 01_Remove_Instrument_Response.sh
│   ├── 02_Split_Acceleration_SAC_to_Daily.py
│   └── 03_Split_Daily_SAC_to_Hourly.sh
│
├── 02_Calculation/
│   └── CalPSD_main.m
│
└── 03_Figures/
    ├── Fig2_Waveform.py
    ├── Fig2_RMS.m
    ├── Fig3_PSD_PSDDiff.m
    └── Fig4_PSDDiff_CurrentSpeed.m
```

The functions required for PSD calculation are included as local functions inside `CalPSD_main.m`; no separate MATLAB function directory is required.

## Associated data repository

The software and data directories are expected to be adjacent:

```text
Test_SUSTech/
├── ESS_OBS_Shield_Data_v1.0/
└── ESS_OBS_Shield_Scripts_v1.0/
```

The scripts use the following archived data products:

```text
ESS_OBS_Shield_Data_v1.0/
├── Pro_Seismic_acceleration_v1.0/
├── Processed_PSD_v1.0/
├── Processed_PSDDiff_v1.0/
├── Processed_RMS_v1.0/
├── Raw_Current_Tide_v1.0/
├── Raw_Metadata_v1.0/
└── Raw_Seismic_mseed_v1.0/
```

Insert the identifiers after publication:

```text
DATA DOI: <10.5281/zenodo.21501794>
SOFTWARE DOI: <10.5281/zenodo.21511732>
```

## Important data-level distinction

`Pro_Seismic_acceleration_v1.0` contains continuous seismic records that have already been corrected for instrument response and converted to acceleration. These data must not be processed again with `01_Remove_Instrument_Response.sh`.

The directly reproducible PSD workflow distributed with this release begins with:

```text
Pro_Seismic_acceleration_v1.0
```

The response-removal script is retained as a provenance script. It requires uncorrected continuous SAC files as input. The archived raw miniSEED files cannot be passed directly to this SAC script without an additional miniSEED-to-continuous-SAC conversion and merging step.

## Software requirements

### Python

- Python 3.9 or later
- ObsPy
- NumPy
- Matplotlib

Install the direct Python dependencies from the repository root:

```bash
python -m pip install -r requirements.txt
```

### MATLAB

MATLAB R2016b or later is recommended because the scripts use local functions at the end of script files. The Statistics and Machine Learning Toolbox may be required for functions such as `prctile`.

### SAC

The preprocessing shell scripts require:

- Seismic Analysis Code (`sac`)
- `saclst`
- Bash

Ensure that both commands are available in `PATH`:

```bash
command -v sac
command -v saclst
```

## Processing workflow

Run all commands from the root of `ESS_OBS_Shield_Scripts_v1.0`.

### Step 1: optional instrument-response removal

Script:

```text
01_Preprocessing/01_Remove_Instrument_Response.sh
```

Purpose:

- removes the instrument response from uncorrected continuous SAC files;
- converts the seismic channels to acceleration;
- applies the original processing sequence:

```text
rmean
rtr
taper
trans from polszeros ... to acc freq 0.001 0.005 40 45
```

Response file:

```text
../ESS_OBS_Shield_Data_v1.0/
└── Raw_Metadata_v1.0/
    └── Instrument-response.SACPZ
```

This step is not required when using the archived `Pro_Seismic_acceleration_v1.0` files.

```bash
chmod +x 01_Preprocessing/01_Remove_Instrument_Response.sh
```

Run only when uncorrected continuous SAC files are available:

```bash
./01_Preprocessing/01_Remove_Instrument_Response.sh \
    /path/to/uncorrected_continuous_SAC \
    ./working_data/Pro_Seismic_acceleration_v1.0 \
    ../ESS_OBS_Shield_Data_v1.0/Raw_Metadata_v1.0/Instrument-response.SACPZ
```

### Step 2: split continuous acceleration SAC into daily files

Script:

```text
01_Preprocessing/02_Split_Acceleration_SAC_to_Daily.py
```

Default input:

```text
../ESS_OBS_Shield_Data_v1.0/
└── Pro_Seismic_acceleration_v1.0/
```

Default output:

```text
working_data/SAC_Day/
└── YYYYMM/YYYYMMDD/STATION/
```

Only complete UTC calendar days are written. Incomplete edge days are excluded.

```bash
python 01_Preprocessing/02_Split_Acceleration_SAC_to_Daily.py
```

### Step 3: split daily SAC into hourly PSD-input files

Script:

```text
01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh
```

Default input:

```text
working_data/SAC_Day/
```

Default output:

```text
working_data/Hourly_Seismic_acceleration_v1.0/
└── YYYYMM/YYYYMMDD/STATION/STATION.COMPONENT/
```

Each hourly segment retains the original operations:

```text
rmean
rtr
taper type cosine width 0.1
```

The original SAC cutting intervals are retained:

```text
0–3600 s
3600–7200 s
...
82800–86400 s
```

```bash
chmod +x 01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh
./01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh
```

### Step 4: calculate hourly PSD

Script:

```text
02_Calculation/CalPSD_main.m
```

Input:

```text
working_data/Hourly_Seismic_acceleration_v1.0/
```

The calculation follows the original workflow:

1. loop over dates, stations, components, and hourly SAC files;
2. calculate the one-sided FFT amplitude spectrum;
3. apply logarithmic frequency smoothing;
4. interpolate to the prescribed frequency vector;
5. convert the amplitude spectrum to PSD;
6. save the PSD matrices.

```bash
matlab -batch "run('02_Calculation/CalPSD_main.m')"
```

Verify that the generated products agree with the archived files in:

```text
../ESS_OBS_Shield_Data_v1.0/Processed_PSD_v1.0/
```

## Figure-generation workflow

### Figure 2: continuous waveforms and tide

Script:

```text
03_Figures/Fig2_Waveform.py
```

Input:

```text
../ESS_OBS_Shield_Data_v1.0/
├── Pro_Seismic_acceleration_v1.0/
└── Raw_Current_Tide_v1.0/Tide_UTC.txt
```

The plotted seismic waveforms are response-corrected acceleration records.

```bash
python 03_Figures/Fig2_Waveform.py
```

### Figure 2: RMS calculation and plots

Script:

```text
03_Figures/Fig2_RMS.m
```

Input:

```text
../ESS_OBS_Shield_Data_v1.0/
├── Pro_Seismic_acceleration_v1.0/
└── Raw_Current_Tide_v1.0/
    ├── dataforzzh.mat
    └── Tide_UTC.txt
```

Principal settings:

- RMS window length: 1800 s;
- expected sampling rate: 100 Hz;
- horizontal resultant: `sqrt(E^2 + N^2)`;
- linear detrending within each window;
- no additional taper before RMS calculation;
- current speed interpolated to the start time of each RMS window.

```bash
matlab -batch "run('03_Figures/Fig2_RMS.m')"
```

### Figure 3: PSD and PSD difference

Script:

```text
03_Figures/Fig3_PSD_PSDDiff.m
```

Input:

```text
../ESS_OBS_Shield_Data_v1.0/Processed_PSD_v1.0/
```

The right panel uses:

```text
PSDdiff = PSD(T02) - PSD(T01)
```

Negative values indicate lower spectral power at the shielded OBS.

```bash
matlab -batch "run('03_Figures/Fig3_PSD_PSDDiff.m')"
```

### Figure 4: PSD difference versus current speed

Script:

```text
03_Figures/Fig4_PSDDiff_CurrentSpeed.m
```

Input:

```text
../ESS_OBS_Shield_Data_v1.0/
├── Processed_PSDDiff_v1.0/
└── Raw_Current_Tide_v1.0/
    ├── Flow_Speed_UTC.txt
    ├── FloodEvents.txt
    └── EbbEvents.txt
```

The script compares flood and ebb conditions in three frequency bands:

- 0.01–0.1 Hz
- 0.1–10 Hz
- 10–40 Hz

```bash
matlab -batch "run('03_Figures/Fig4_PSDDiff_CurrentSpeed.m')"
```

## Units and time standard

- Seismic acceleration: `m/s^2`
- PSD: `dB re (m/s^2)^2/Hz`
- Current speed: `m/s`
- Tide level: `m`
- Time standard: UTC

The current record stored in `dataforzzh.mat` is converted from Beijing time to UTC inside `Fig2_RMS.m`. Text current and tide products whose filenames contain `_UTC` are treated as UTC inputs.

## Intermediate and generated files

The following directories are created locally and are intentionally excluded from version control:

```text
working_data/
outputs/
```

They contain reproducible intermediate waveforms, calculated products, and figure files.

## Reproducibility checks before release

Check for remaining local absolute paths:

```bash
grep -R "/home/lyang" .
grep -R "/Disk1\|/Disk2" .
```

Check Python syntax:

```bash
python -m py_compile \
    01_Preprocessing/02_Split_Acceleration_SAC_to_Daily.py \
    03_Figures/Fig2_Waveform.py
```

Check shell syntax:

```bash
bash -n 01_Preprocessing/01_Remove_Instrument_Response.sh
bash -n 01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh
```

Check MATLAB scripts:

```bash
matlab -batch "checkcode('02_Calculation/CalPSD_main.m')"
matlab -batch "checkcode('03_Figures/Fig2_RMS.m')"
matlab -batch "checkcode('03_Figures/Fig3_PSD_PSDDiff.m')"
matlab -batch "checkcode('03_Figures/Fig4_PSDDiff_CurrentSpeed.m')"
```

Because the scripts are stored in subdirectories, confirm that each script resolves the data directory as a sibling of `ESS_OBS_Shield_Scripts_v1.0`, rather than as a child of the script repository.

## Citation

After publication of the Zenodo software record, cite the archived software version using the metadata in `CITATION.cff`.

Suggested citation template:

```text
Li, Y. (2026). Processing and Figure-Generation Scripts for Evaluating a Double-Spherical-Shell Shield for Shallow-Water Ocean-Bottom Seismic Noise Reduction (Version 1.0) [Software]. Zenodo. <10.5281/zenodo.21508388>
```

The associated dataset should be cited separately using its dataset DOI.

## License

This software is released under the MIT License. See `LICENSE`.
