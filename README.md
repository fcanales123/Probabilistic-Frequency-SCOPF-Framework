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

- **MATPOWER v8.0 or later** (required)
- **Gurobi v13.0.2 o later**  (optional but recommended)

Make sure MATPOWER is installed and added to your MATLAB path:

## License
This work is licensed under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
