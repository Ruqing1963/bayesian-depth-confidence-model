#!/usr/bin/env python3
"""
Bayesian Depth-Confidence Model for Deep Orogenic Gold Exploration
Using Logistic Survival Functions Calibrated to the Fault-Valve Mechanism

Author: Ruqing Chen, GUT Geoservice Inc., Montreal, Canada
Date:   March 2026

This script computes the posterior probability P(Large | Hit, H) that a
drill-core intercept of gold mineralization at depth H belongs to a
super-large (>=100 t Au) orogenic deposit system.

Usage:
    python logistic_model.py
    python logistic_model.py --depth 2500 --ks 0.005
    python logistic_model.py --sweep
"""

import argparse
import json
import numpy as np
import os


# ============================================================
# Model Parameters (default values)
# ============================================================
DEFAULT_PARAMS = {
    "prior": {
        "P_large": 0.01,
        "description": "Targeted-site base rate: fraction of sites with surface anomalies that prove to host super-large systems"
    },
    "large_system": {
        "A_L": 0.40,
        "H_c_L": 4000,
        "k_L": 0.005,
        "description": "Crustal-scale orogenic system: 70-deg dip orebody, fault-valve sustained conduit"
    },
    "small_system": {
        "A_S": 0.20,
        "H_c_S": 1000,
        "k_S": 0.005,
        "k_S_range": [0.003, 0.008],
        "description": "Sub-economic passive system: ambient permeability dependent, sealed by BDT"
    },
    "paleo_depth": {
        "H_paleo_crit": 6000,
        "Delta_erosion": 5000,
        "H_c_S_derived": 1000,
        "formula": "H_c_S = H_paleo_crit - Delta_erosion",
        "description": "Critical closure depth derived from Archean paleo-depth minus post-orogenic erosion"
    }
}


def logistic_survival(H, A, H_c, k):
    """
    Logistic survival function.

    P(Hit | Type) = A / [1 + exp(k * (H - H_c))]

    Parameters
    ----------
    H : float or np.ndarray
        Present-day drill depth (meters).
    A : float
        Geometric intercept probability at H = 0.
    H_c : float
        Critical closure depth (meters) — half-survival depth.
    k : float
        Steepness parameter. Transition zone width: Delta_Z ≈ 4.4 / k.

    Returns
    -------
    float or np.ndarray
        Survival probability at depth H.
    """
    return A / (1 + np.exp(k * (H - H_c)))


def posterior_probability(H, A_L=0.40, H_c_L=4000, k_L=0.005,
                          A_S=0.20, H_c_S=1000, k_S=0.005,
                          P_large=0.01):
    """
    Compute Bayesian posterior P(Large | Hit, H).

    Parameters
    ----------
    H : float or np.ndarray
        Present-day drill depth (meters).
    A_L, H_c_L, k_L : float
        Logistic parameters for large systems.
    A_S, H_c_S, k_S : float
        Logistic parameters for small systems.
    P_large : float
        Prior probability of large system at targeted site.

    Returns
    -------
    float or np.ndarray
        Posterior probability P(Large | Hit, H).
    """
    P_small = 1 - P_large
    P_hit_L = logistic_survival(H, A_L, H_c_L, k_L)
    P_hit_S = logistic_survival(H, A_S, H_c_S, k_S)

    numerator = P_hit_L * P_large
    denominator = numerator + P_hit_S * P_small

    return numerator / denominator


def transition_zone_width(k):
    """
    Compute the 90%-to-10% transition zone width.

    Delta_Z ≈ 4.4 / k (derived from logistic function properties).

    Parameters
    ----------
    k : float
        Steepness parameter.

    Returns
    -------
    float
        Transition zone width in meters.
    """
    return 4.394 / k  # exact: 2 * ln(9) / k


def print_single_depth(H, k_S=0.005):
    """Print detailed results for a single depth."""
    params = DEFAULT_PARAMS

    P_hit_L = logistic_survival(H,
                                params["large_system"]["A_L"],
                                params["large_system"]["H_c_L"],
                                params["large_system"]["k_L"])
    P_hit_S = logistic_survival(H,
                                params["small_system"]["A_S"],
                                params["small_system"]["H_c_S"],
                                k_S)
    P_post = posterior_probability(H, k_S=k_S)

    print(f"\n{'='*60}")
    print(f"  BAYESIAN DEPTH-CONFIDENCE MODEL")
    print(f"  Depth H = {H} m  |  k_S = {k_S}")
    print(f"{'='*60}")
    print(f"  Prior P(Large)        = {params['prior']['P_large']}")
    print(f"  P(Hit | Large)        = {P_hit_L:.6f}")
    print(f"  P(Hit | Small)        = {P_hit_S:.6f}")
    print(f"  Likelihood ratio      = {P_hit_L/P_hit_S:.1f}x")
    print(f"  ----------------------------------------")
    print(f"  >>> P(Large | Hit)    = {P_post*100:.1f}%")
    print(f"{'='*60}")
    print(f"  Transition zone width = {transition_zone_width(k_S):.0f} m")
    print()


def print_sweep():
    """Print full depth sweep table with sensitivity analysis."""
    depths = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500]
    k_values = [0.003, 0.005, 0.008]

    print(f"\n{'='*72}")
    print(f"  DEPTH SWEEP: Posterior P(Large | Hit) [%]")
    print(f"  Prior = {DEFAULT_PARAMS['prior']['P_large']}, "
          f"A_L = {DEFAULT_PARAMS['large_system']['A_L']}, "
          f"A_S = {DEFAULT_PARAMS['small_system']['A_S']}")
    print(f"{'='*72}")

    header = f"  {'k_S':>14s}"
    for d in depths:
        header += f"  {d:>6d}m"
    print(header)
    print(f"  {'-'*68}")

    for k_S in k_values:
        label = f"k_S={k_S}"
        if k_S == 0.005:
            label += " *"
        row = f"  {label:>14s}"
        for d in depths:
            p = posterior_probability(d, k_S=k_S) * 100
            row += f"  {p:>6.1f}%"
        print(row)

    print(f"\n  * = central estimate")
    print(f"  Transition zone widths: "
          f"k=0.003 -> {transition_zone_width(0.003):.0f} m, "
          f"k=0.005 -> {transition_zone_width(0.005):.0f} m, "
          f"k=0.008 -> {transition_zone_width(0.008):.0f} m")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Bayesian Depth-Confidence Model for Deep Orogenic Gold Exploration"
    )
    parser.add_argument("--depth", type=float, default=2500,
                        help="Drill depth in meters (default: 2500)")
    parser.add_argument("--ks", type=float, default=0.005,
                        help="Small-system steepness k_S (default: 0.005)")
    parser.add_argument("--sweep", action="store_true",
                        help="Print full depth sweep table")
    parser.add_argument("--export-params", action="store_true",
                        help="Export parameters to data/parameters.json")

    args = parser.parse_args()

    if args.export_params:
        outpath = os.path.join(os.path.dirname(__file__), "..", "data", "parameters.json")
        os.makedirs(os.path.dirname(outpath), exist_ok=True)
        with open(outpath, "w") as f:
            json.dump(DEFAULT_PARAMS, f, indent=2)
        print(f"Parameters exported to {outpath}")
        return

    if args.sweep:
        print_sweep()
    else:
        print_single_depth(args.depth, args.ks)


if __name__ == "__main__":
    main()
