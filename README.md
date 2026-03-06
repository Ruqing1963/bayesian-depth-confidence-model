# Bayesian Depth–Confidence Model for Deep Orogenic Gold Exploration

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A Bayesian posterior probability model that quantifies how drill depth transforms the diagnostic power of a gold mineralization intercept. The model uses **logistic survival functions** calibrated to the fault-valve mechanism (Sibson, 1988) to demonstrate that lithostatic pressure acts as a natural high-pass filter, selectively eliminating sub-economic fluid systems while preserving crustal-scale orogenic gold conduits.

## Key Result

At a targeted exploration site with prior P₀(Large) = 0.01, a single drill intercept at **2500 m depth** yields **65–100% posterior confidence** (central: **97.3%**) of connecting to a super-large (≥100 t Au) deposit system — compared to only 2% at the surface.

<p align="center">
  <img src="figures/fig1_logistic_model.png" width="750">
</p>

## Model Overview

### Logistic Survival Functions

Unlike exponential decay formulations, the logistic model naturally encodes the brittle–ductile transition (BDT) as a **threshold** rather than a gradient:

$$P(\text{Hit} \mid \text{Type}) = \frac{A_{\text{Type}}}{1 + \exp\!\left(k_{\text{Type}} \cdot (H - H_{c,\text{Type}})\right)}$$

| Parameter | Large System | Small System | Physical Basis |
|-----------|:---:|:---:|---|
| *A* (geometric prefactor) | 0.40 | 0.20 | Drill intercept probability at surface |
| *H<sub>c</sub>* (critical closure depth) | 4000 m | 1000 m | Fault-valve capacity / ambient permeability closure |
| *k* (steepness) | 0.005 | 0.003–0.008 | BDT transition zone thickness (~880 m for *k* = 0.005) |

### Paleo-Depth Calibration

The critical closure depth is derived from first principles:

$$H_c = H_{\text{paleo,crit}} - \Delta_{\text{erosion}}$$

For the Abitibi belt: *H*<sub>paleo,crit</sub> ≈ 6000 m, Δ<sub>erosion</sub> ≈ 5000 m → *H*<sub>c,S</sub> ≈ 1000 m.

### Posterior Probability Table

| *k<sub>S</sub>* | H = 0 | H = 1000 m | H = 1500 m | H = 2000 m | H = 2500 m |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.003 (gentle) | 2.0% | 2.5% | 5.5% | 20% | 65% |
| **0.005 (central)** | **2.0%** | **3.9%** | **21%** | **76%** | **97.3%** |
| 0.008 (steep) | 2.0% | 9.6% | 72% | 99.3% | ~100% |

## Repository Structure

```
bayesian-depth-confidence-model/
├── README.md
├── LICENSE
├── paper/
│   ├── paper.tex              # LaTeX source
│   ├── paper.pdf              # Compiled PDF
│   └── fig1_latex.pdf         # Figure (vector, for LaTeX)
├── figures/
│   └── fig1_logistic_model.png  # Figure (300 DPI PNG)
├── scripts/
│   ├── logistic_model.py      # Core Bayesian computation
│   └── generate_figure.py     # Figure generation script
└── data/
    └── parameters.json        # Model parameters
```

## Quick Start

```bash
# Clone the repository
git clone https://github.com/Ruqing1963/bayesian-depth-confidence-model.git
cd bayesian-depth-confidence-model

# Run the model (requires Python 3.7+ with numpy)
python scripts/logistic_model.py

# Regenerate the figure (requires matplotlib)
python scripts/generate_figure.py
```

## Dependencies

- Python ≥ 3.7
- NumPy
- Matplotlib (for figure generation only)
- LaTeX distribution (for paper compilation only)

## Case Study: LaRonde Mine, Abitibi Belt

LaRonde is an Au-rich VMS deposit overprinted by orogenic deformation, mined to ~3100 m. At depth, economic gold is confined to **2 structural-stratigraphic corridors** out of >20 prospective horizons identified at surface — a ~90% elimination rate by depth-filtering. The model yields P(Large | Hit, H = 3100 m) = 99.5%.

## Citation

If you use this model in your research, please cite:

```bibtex
@article{Chen2026depth,
  author  = {Ruqing Chen},
  title   = {A {Bayesian} Depth--Confidence Model for Deep Orogenic Gold
             Exploration Using Logistic Survival Functions Calibrated to
             the Fault-Valve Mechanism},
  year    = {2026},
  note    = {GUT Geoservice Inc., Montreal, Canada},
  url     = {https://github.com/Ruqing1963/bayesian-depth-confidence-model}
}
```

The number-theoretic investigation that inspired this geological model:

```bibtex
@article{Chen2026primes,
  author  = {Ruqing Chen},
  title   = {Prime Constellations of $n^k - (n-1)^k$: Algebraic Obstructions,
             {Bateman--Horn} Verification, and the Non-Monotonicity of Maximum
             Constellation Lengths},
  year    = {2026},
  doi     = {10.5281/zenodo.18819849},
  url     = {https://zenodo.org/records/18819849}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Author

**Ruqing Chen**
GUT Geoservice Inc., Montreal, Canada
ruqing@hotmail.com
