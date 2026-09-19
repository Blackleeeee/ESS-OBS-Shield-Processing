# ESS OBS Shield Processing and Figure Scripts

## Overview

This repository contains the MATLAB, Python, and Bash scripts used to process paired shallow-water ocean-bottom seismic observations and reproduce the analyses in the revised manuscript *Tidal Modulation of Flow-Induced Noise in Shallow-Water Ocean-Bottom Seismic Records: Evidence From a Paired Shielding Experiment*.

Version **v1.1.0** retains the preprocessing and PSD-calculation workflow from the previous release and updates the figure-generation workflow to match the final revised manuscript. The revised release covers Figures 2–5, Figures S1–S5, and the P-/S-wave SNR calculations reported in Table S1.

## Repository structure

```text
ESS-OBS-Shield-Processing/
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── .gitignore
├── 01_Preprocessing/
│   ├── 01_Remove_Instrument_Response.sh
│   ├── 02_Split_Acceleration_SAC_to_Daily.py
│   └── 03_Split_Daily_SAC_to_Hourly.sh
├── 02_Calculation/
│   └── CalPSD_main.m
└── 03_Figures/
    ├── Figure_2.py
    ├── Figure_3.m
    ├── Figure_4.m
    ├── Figure_5&S4.py
    ├── Figure_S1.m
    ├── Figure_S2.m
    ├── Figure_S3.m
    └── Figure_S5.py
```

Generated figures and intermediate products are written under `outputs/` and are excluded from version control.

## Associated data repository

The final figure scripts expect the software and expanded dataset directories to be adjacent:

```text
workspace/
├── ESS-OBS-Shield-Processing/
└── ESS_OBS_Shield_Data_v1.1.0/
```

The dataset contains:

```text
ESS_OBS_Shield_Data_v1.1.0/
├── Raw_Seismic_mseed_v1.0/
├── Pro_Seismic_acceleration_v1.0/
├── Processed_Earthquake_velocity_v1.0/
├── Raw_Current_Tide_v1.1/
├── Raw_Metadata_v1.0/
├── Processed_RMS_v1.0/
├── Processed_PSD_v1.0/
└── Processed_PSDDiff_v1.0/
```

## Associated data repository

The final figure scripts are designed to reproduce the analyses using the associated dataset release:

Li, Y. et al. (2026).  
Seismic and Hydrodynamic Data From a Paired Shielding Experiment on Flow-Induced Noise in Shallow-Water Ocean-Bottom Seismic Records (Version v1.1.0) [Dataset]. Zenodo.

https://doi.org/10.5281/zenodo.22843942

Previous dataset version:
https://doi.org/10.5281/zenodo.21501794

## Final manuscript mapping

| Manuscript item | Script | Direct data inputs |
|---|---|---|
| Figure 2 | `03_Figures/Figure_2.py` | continuous acceleration SAC; `Results_RMS.mat`; tide level; current speed; current direction |
| Figure 3 | `03_Figures/Figure_3.m` | T01/T02/ZSD02 PSD matrices; PSDdiff calculated internally |
| Figure 4 | `03_Figures/Figure_4.m` | T01/T02 PSD matrices; current speed; flood/ebb interval files |
| Figure 5 | `03_Figures/Figure_5&S4.py` | response-corrected earthquake velocity SAC files |
| Figure S1 | `03_Figures/Figure_S1.m` | PSD matrices; PSDdiff matrices; tide level |
| Figure S2 | `03_Figures/Figure_S2.m` | continuous acceleration SAC; vertical PSDdiff |
| Figure S3 | `03_Figures/Figure_S3.m` | PSD matrices; current speed; current direction; flood/ebb intervals |
| Figure S4 | `03_Figures/Figure_5&S4.py` | response-corrected earthquake velocity SAC files |
| Figure S5 | `03_Figures/Figure_S5.py` | T01/T02 vertical PSD; vertical PSDdiff; current speed; tide level |
| Table S1 | `03_Figures/Figure_5&S4.py` | response-corrected earthquake velocity SAC files |

Figure 1 is a site/deployment figure and has no dedicated analysis script in this repository.

## Important processing definitions

### Horizontal total PSD

Where the manuscript uses horizontal-total power, the two recorded horizontal PSDs are summed in the linear domain before returning to dB:

```text
PSD_HH = 10 log10(10^(PSD_HH1/10) + 10^(PSD_HH2/10))
```

This quantity is invariant to horizontal-axis rotation.

### PSD difference

```text
PSDdiff = PSD_T02 - PSD_T01
```

Negative values indicate lower spectral power at the shielded OBS T02 than at the unshielded OBS T01.

### P- and S-wave SNR

`Figure_5&S4.py` uses response-corrected ground-velocity records. P-wave SNR is calculated on HHZ. S-wave SNR uses the rotation-invariant total horizontal RMS:

```text
RMS_H = sqrt(mean(HH1^2 + HH2^2))
SNR_dB = 20 log10(RMS_signal / RMS_noise)
```

T01 and T02 use identical filters and time windows.

## Figure scripts

### Figure 2

`Figure_2.py` combines current speed and direction, three-component continuous acceleration waveforms, tide level, precomputed 30-min RMS values, and RMS-current relationships in the final Figure 2 layout.

Direct inputs:

```text
ESS_OBS_Shield_Data_v1.1.0/
├── Pro_Seismic_acceleration_v1.0/
├── Processed_RMS_v1.0/Results_RMS.mat
└── Raw_Current_Tide_v1.1/
    ├── Tide_UTC.txt
    ├── Flow_Speed_UTC.txt
    └── Flow_Direction_UTC.txt
```

### Figure 3

`Figure_3.m` calculates whole-period PSD medians and IQRs for T01, T02, and ZSD02. Horizontal-total PSD is formed from the two horizontal channels in linear power. The T02–T01 PSD difference is calculated internally from the paired PSD matrices.

### Figure 4

`Figure_4.m` compares PSDdiff with measured surface-current speed during flood and ebb tides in the 0.01–0.1, 0.1–10, and 10–40 Hz bands for HHZ and HH.

Direct environmental inputs:

```text
Raw_Current_Tide_v1.1/
├── Flow_Speed_UTC.txt
├── FloodEvents.txt
└── EbbEvents.txt
```

### Figure 5, Figure S4, and Table S1

`Figure_5&S4.py` uses the six response-corrected velocity SAC files for the 28 March 2025 Mw 7.7 earthquake. It generates the main earthquake comparison, additional filtered waveform comparisons, and the P-/S-wave SNR products used for Table S1.

### Figure S1

`Figure_S1.m` generates ambient-noise PSD and PSDdiff spectrograms with tide-level overlays.

### Figure S2

`Figure_S2.m` calculates magnitude-squared coherence within T01 and T02 and vertical coherence between the two OBSs using 1-h windows divided into nonoverlapping 5-min subwindows. The final panel overlays vertical PSDdiff.

### Figure S3

`Figure_S3.m` compares absolute T01/T02 PSD distributions for jointly constrained tidal-stage and current-direction groups in common current-speed bins. Flood-tide samples are associated with the ~300° current-direction sector and ebb-tide samples with the ~120° sector. Frequency-time PSD samples are pooled within each accepted speed bin, while the minimum-sample gate is based on independent hourly records.

### Figure S5

`Figure_S5.py` characterizes the persistent high-frequency vertical spectral peak at T02, including hourly spectra, tracked peak frequency, peak PSD with tide level, and peak PSD versus measured current speed.

## Software requirements

### Python

- Python 3.9 or later
- NumPy
- Matplotlib
- ObsPy
- SciPy
- pandas

Install the direct Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

### MATLAB

MATLAB R2016b or later is recommended. The Statistics and Machine Learning Toolbox may be required for percentile calculations.

### SAC

The preprocessing shell scripts require Seismic Analysis Code (`sac`), `saclst`, and Bash.

## Preprocessing and PSD calculation

The scripts in `01_Preprocessing/` and `02_Calculation/` are preserved from the previous software release without scientific modification. They document the original response-removal, segmentation, and PSD-generation workflow. The final revised figure scripts consume the publicly released processed products directly.

## Reproducibility checks

Before release, verify that no local absolute paths remain in `03_Figures/`:

```bash
grep -R "/home/lyang\|/Disk1\|/Disk2" 03_Figures
```

Check Python syntax:

```bash
python -m py_compile \
    03_Figures/Figure_2.py \
    '03_Figures/Figure_5&S4.py' \
    03_Figures/Figure_S5.py
```

## Version 1.1.0 changes

- retained the existing preprocessing and PSD-calculation scripts;
- replaced the previous Figure 2–4 plotting scripts with the final revised versions;
- added Figure 5 and Figure S4 earthquake waveform and SNR analysis;
- added Figure S1 spectrogram and PSDdiff analysis;
- added Figure S2 coherence analysis;
- added Figure S3 flood–ebb/current-direction-constrained PSD comparison;
- added Figure S5 persistent high-frequency peak analysis;
- updated Python dependencies to include SciPy and pandas;
- replaced local figure-script input paths with repository-relative paths.

## Citation

Li, Y. et al. (2026).  
Processing and Figure Generation Scripts for a Paired Shielding Experiment on Flow-Induced Noise in Shallow-Water Ocean-Bottom Seismic Records (Version v1.1.0) [Software]. Zenodo.

https://doi.org/xxxxx

## License

This software is released under the MIT License. See `LICENSE`.
