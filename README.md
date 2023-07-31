# V2DL5 - high-level analysis for VERITAS with gammapy

This is a collection of simple scripts for the high-level analysis of VERITAS data with gammapy. 

This includes for a given list of runs or for a cone search around a given direction or named source. 

- source detection analysis including integral flux (or flux upper limits), reflection region model
- spectral analysis, reflected region model

## Installation

To install the required python packages, run:

```
conda env create -f environment.yml
```

Activate the environment to start the analysis:
```
conda activate v2dl5
```
