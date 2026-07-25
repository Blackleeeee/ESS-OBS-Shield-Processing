ESS OBS Shield Processing and Figure-Generation Scripts

1. Overview

This repository contains the MATLAB, Python, and Bash scripts used to process continuous ocean-bottom seismic records and reproduce the root mean square (RMS), power spectral density (PSD), PSD-difference, and current-dependent analyses used to evaluate a double-spherical-shell shield for shallow-water ocean-bottom seismometers.

The directly reproducible analysis workflow begins with continuous seismic waveforms that have already been corrected for instrument response and converted to acceleration.

2. Repository Structure

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

The functions required for PSD calculation are included as local functions inside CalPSD_main.m; no separate MATLAB function directory is required.

3. Persistent Identifiers and Source Repository

Dataset: https://doi.org/10.5281/zenodo.21501794

Software release:  https://doi.org/10.5281/zenodo.21543745

GitHub repository: https://github.com/Blackleeeee/ESS-OBS-Shield-Processing

Software version: 1.0.0

Software license: MIT License

4. Associated Data Repository

The software and expanded data directories are expected to be adjacent:

Test_SUSTech/
├── ESS_OBS_Shield_Data_v1.0/
└── ESS_OBS_Shield_Scripts_v1.0/

The scripts use the following archived data products:

ESS_OBS_Shield_Data_v1.0/
├── Pro_Seismic_acceleration_v1.0/
├── Processed_PSD_v1.0/
├── Processed_PSDDiff_v1.0/
├── Processed_RMS_v1.0/
├── Raw_Current_Tide_v1.0/
├── Raw_Metadata_v1.0/
└── Raw_Seismic_mseed_v1.0/

For Zenodo distribution, each data directory is provided as a separate .tar.gz archive. The archives must be extracted before running the scripts.

5. Important Data-Level Distinction

Pro_Seismic_acceleration_v1.0 contains continuous seismic records that have already been corrected for instrument response and converted to acceleration. These data must not be processed again with 01_Remove_Instrument_Response.sh.

The directly reproducible workflow distributed with this release begins with:

Pro_Seismic_acceleration_v1.0

The response-removal script is retained as a provenance script. It requires uncorrected continuous SAC files as input. The archived raw miniSEED files cannot be passed directly to this SAC script without an additional miniSEED-to-continuous-SAC conversion and merging step, which is not included in this release.

6. Software Requirements

6.1 Python

Python 3.9 or later

ObsPy

NumPy

Matplotlib

Install the direct Python dependencies from the repository root:

python -m pip install -r requirements.txt

6.2 MATLAB

MATLAB R2016b or later is recommended because the MATLAB scripts contain local functions at the end of script files. The Statistics and Machine Learning Toolbox may be required for functions such as prctile.

6.3 SAC

The preprocessing shell scripts require:

Seismic Analysis Code (sac);

saclst;

Bash.

Confirm that the SAC executables are available:

command -v sac
command -v saclst

7. Processing Workflow

Run the following commands from the root directory of ESS_OBS_Shield_Scripts_v1.0.

7.1 Optional Instrument-Response Removal

Script:

01_Preprocessing/01_Remove_Instrument_Response.sh

Purpose:

remove the instrument response from uncorrected continuous SAC files;

convert the seismic records to acceleration;

apply the original processing sequence:

rmean
rtr
taper
trans from polszeros ... to acc freq 0.001 0.005 40 45

Response file:

../ESS_OBS_Shield_Data_v1.0/
└── Raw_Metadata_v1.0/
    └── Instrument-response.SACPZ

This step is not required when using the archived Pro_Seismic_acceleration_v1.0 files.

Make the script executable:

chmod +x 01_Preprocessing/01_Remove_Instrument_Response.sh

Run only when uncorrected continuous SAC files are available:

./01_Preprocessing/01_Remove_Instrument_Response.sh     /path/to/uncorrected_continuous_SAC     ./working_data/Pro_Seismic_acceleration_v1.0     ../ESS_OBS_Shield_Data_v1.0/Raw_Metadata_v1.0/Instrument-response.SACPZ

7.2 Split Continuous Acceleration SAC Into Daily Files

Script:

01_Preprocessing/02_Split_Acceleration_SAC_to_Daily.py

Default input:

../ESS_OBS_Shield_Data_v1.0/
└── Pro_Seismic_acceleration_v1.0/

Default output:

working_data/SAC_Day/
└── YYYYMM/YYYYMMDD/STATION/

Only complete UTC calendar days are written. Incomplete edge days are excluded.

Run:

python 01_Preprocessing/02_Split_Acceleration_SAC_to_Daily.py

7.3 Split Daily SAC Into Hourly PSD-Input Files

Script:

01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh

Default input:

working_data/SAC_Day/

Default output:

working_data/Hourly_Seismic_acceleration_v1.0/
└── YYYYMM/YYYYMMDD/STATION/STATION.COMPONENT/

Each hourly segment retains the original operations:

rmean
rtr
taper type cosine width 0.1

The original SAC cutting intervals are retained:

0–3600 s
3600–7200 s
...
82800–86400 s

Run:

chmod +x 01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh
./01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh

7.4 Calculate Hourly PSD

Script:

02_Calculation/CalPSD_main.m

Input:

working_data/Hourly_Seismic_acceleration_v1.0/

The calculation follows the original workflow:

loop over dates, stations, components, and hourly SAC files;

calculate the one-sided FFT amplitude spectrum;

apply logarithmic frequency smoothing;

interpolate to the prescribed frequency vector;

convert the amplitude spectrum to PSD;

save the PSD matrices.

Run:

matlab -batch "run('02_Calculation/CalPSD_main.m')"

The generated PSD products can be compared with:

../ESS_OBS_Shield_Data_v1.0/Processed_PSD_v1.0/

8. Figure-Generation Workflow

8.1 Figure 2: Continuous Waveforms and Tide

Script:

03_Figures/Fig2_Waveform.py

Inputs:

../ESS_OBS_Shield_Data_v1.0/
├── Pro_Seismic_acceleration_v1.0/
└── Raw_Current_Tide_v1.0/Tide_UTC.txt

The plotted seismic waveforms are response-corrected acceleration records.

Run:

python 03_Figures/Fig2_Waveform.py

8.2 Figure 2: RMS Calculation and Plots

Script:

03_Figures/Fig2_RMS.m

Inputs:

../ESS_OBS_Shield_Data_v1.0/
├── Pro_Seismic_acceleration_v1.0/
└── Raw_Current_Tide_v1.0/
    ├── dataforzzh.mat
    └── Tide_UTC.txt

Principal settings:

RMS window length: 1,800 s;

window overlap: none;

expected sampling rate: 100 Hz;

horizontal resultant: sqrt(E^2 + N^2);

linear detrending within each window;

no additional taper before RMS calculation;

minimum valid-data fraction: 99%;

current speed interpolated to the start time of each RMS window.

Run:

matlab -batch "run('03_Figures/Fig2_RMS.m')"

8.3 Figure 3: PSD and PSD Difference

Script:

03_Figures/Fig3_PSD_PSDDiff.m

Input:

../ESS_OBS_Shield_Data_v1.0/Processed_PSD_v1.0/

The PSD difference is defined as:

PSDdiff = PSD(T02) - PSD(T01)

Negative values indicate lower spectral power at the shielded OBS.

Run:

matlab -batch "run('03_Figures/Fig3_PSD_PSDDiff.m')"

8.4 Figure 4: PSD Difference Versus Current Speed

Script:

03_Figures/Fig4_PSDDiff_CurrentSpeed.m

Required inputs:

../ESS_OBS_Shield_Data_v1.0/
├── Processed_PSDDiff_v1.0/
└── Raw_Current_Tide_v1.0/
    ├── Flow_Speed_UTC.txt
    ├── FloodEvents.txt
    └── EbbEvents.txt

FloodEvents.txt and EbbEvents.txt contain the start and end times of the flood- and ebb-tide intervals used by the plotting script. The script reads these interval files directly and does not derive them internally.

The script compares flood and ebb conditions in three frequency bands:

0.01–0.1 Hz;

0.1–10 Hz;

10–40 Hz.

Run:

matlab -batch "run('03_Figures/Fig4_PSDDiff_CurrentSpeed.m')"

9. Units and Time Standard

Seismic acceleration: m/s^2;

PSD: dB re (m/s^2)^2/Hz;

PSD difference: dB;

Current speed: m/s;

Tide level: m;

Time standard: UTC.

The current record stored in dataforzzh.mat is converted from Beijing time (UTC+8) to UTC inside Fig2_RMS.m. Text current and tide products whose filenames contain _UTC are treated as UTC inputs. The UTC conversion must not be applied twice.

10. Intermediate and Generated Files

The following directories are created locally and are intentionally excluded from version control:

working_data/
outputs/

They contain reproducible intermediate waveforms, calculated products, and figure files.

11. Validation Checks

Check for remaining local absolute paths:

grep -R "/home/lyang" .
grep -R "/Disk1\|/Disk2" .

Check Python syntax:

python -m py_compile     01_Preprocessing/02_Split_Acceleration_SAC_to_Daily.py     03_Figures/Fig2_Waveform.py

Check shell syntax:

bash -n 01_Preprocessing/01_Remove_Instrument_Response.sh
bash -n 01_Preprocessing/03_Split_Daily_SAC_to_Hourly.sh

Check MATLAB scripts:

matlab -batch "checkcode('02_Calculation/CalPSD_main.m')"
matlab -batch "checkcode('03_Figures/Fig2_RMS.m')"
matlab -batch "checkcode('03_Figures/Fig3_PSD_PSDDiff.m')"
matlab -batch "checkcode('03_Figures/Fig4_PSDDiff_CurrentSpeed.m')"

Because the scripts are stored in subdirectories, each script must resolve the data directory as a sibling of ESS_OBS_Shield_Scripts_v1.0, rather than as a child of the software repository.

12. Citation

Please cite this software as:

Li, Y. (2026). Processing and figure-generation scripts for evaluating a double-spherical-shell shield for shallow-water ocean-bottom seismic noise reduction (Version 1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.21511732

The associated dataset should be cited separately:

Li, Y. (2026). Seismic and hydrodynamic data for evaluating a double-spherical-shell shield in shallow-water ocean-bottom seismic observations (Version 1.0.0) [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.21501794

When the dataset and software are cited together in a reference list, year suffixes such as 2026a and 2026b may be assigned according to the journal's reference-ordering rules.

13. License

This software is released under the MIT License. See LICENSE.
