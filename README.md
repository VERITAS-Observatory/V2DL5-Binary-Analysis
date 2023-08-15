# V2DL5 - high-level analysis for VERITAS with gammapy

[![gammapy](https://img.shields.io/badge/powered%20by-gammapy-orange.svg?style=flat)](https://www.gammapy.org/)
[![GitHub Super-Linter](https://github.com/GernotMaier/V2DL5/actions/workflows/linter.yml/badge.svg)](https://github.com/marketplace/actions/super-linter)

This is a collection of simple scripts for the high-level analysis of VERITAS data with gammapy.

This includes for a given list of runs or for a cone search around a given direction or named source.

- source detection analysis including integral flux (or flux upper limits), reflection region model
- spectral analysis, reflected region model

Reminder on data levels (as defined e.g., by CTAO):

- DL3: event lists for gamma-ray like events
- DL4: binned data (in time, space, energy)
- DL5: science data products (e.g., sky maps, energy spectra, light curves)

This repository provides scripts to generate DL5 science data products.

## Acknowledgement

This work relies heavily on the [gammapy](https://gammapy.org/) development and especially on the excellent [tutorials](https://docs.gammapy.org/1.1/tutorials/index.html) provided by the gammapy team.

## Installation

To install the required python packages, run:

```bash
conda env create -f environment.yml
```

Activate the environment to start the analysis:

```bash
conda activate v2dl5
```
