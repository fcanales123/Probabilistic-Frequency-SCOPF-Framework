# Probabilistic-Frequency-SCOPF-Framework
This repository contains a probabilistic Security-Constrained Optimal Power Flow (SCOPF) framework submitted to: PMAPS 2026.

The source code implements: 

- Models **stochastic wind power output** (including wake effects, availability, and wind variability)
- Captures **correlated flexibility** between wind farms and co-located electrolyzers (WEL plants)
- Simulates **frequency response dynamics** across:
  - Pre-contingency dispatch
  - Frequency Containment Reserves (FCR)
  - Tertiary redispatch (RD)
- Quantifies system security using a **probabilistic congestion-based risk index (ECI)**

---

## ⚙️ Requirements

- MATLAB (tested with recent versions)
- **MATPOWER v8.0 or later** (required)

Make sure MATPOWER is installed and added to your MATLAB path:

```matlab
addpath(genpath('path_to_matpower'));

## License
This work is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
