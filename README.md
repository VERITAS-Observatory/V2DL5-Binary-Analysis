# V2DL5 - high-level analysis for VERITAS with gammapy

This is a collection of simple scripts for the high-level analysis of VERITAS data with gammapy. 

This includes for a given list of runs or for a cone search around a given direction or named source. 

- source detection analysis including integral flux (or flux upper limits), reflection region model
- spectral analysis, reflected region model

Reminder on data levels (as defined e.g., by CTAO):
- DL3: event lists for gamma-ray like events
- DL4: binned data (in time, space, energy)
- DL5: science data products (e.g., sky maps, energy spectra, light curves)

This repository proivdes scripts to generate DL5 science data products.

## Installation

To install the required python packages, run:

```
conda env create -f environment.yml
```

Activate the environment to start the analysis:
```
conda activate v2dl5
```
