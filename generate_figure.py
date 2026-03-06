#!/usr/bin/env python3
"""
Generate Figure 1: Bayesian Depth-Confidence Model with Logistic Survival Functions

Author: Ruqing Chen, GUT Geoservice Inc., Montreal, Canada
Date:   March 2026

Outputs:
    figures/fig1_logistic_model.png  (300 DPI raster)
    paper/fig1_latex.pdf             (vector for LaTeX)

Usage:
    python generate_figure.py
"""

import os
import sys
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    print("Error: matplotlib is required. Install with: pip install matplotlib")
    sys.exit(1)

# Add parent directory for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from logistic_model import logistic_survival, posterior_probability


def generate_figure(output_dir_figures, output_dir_paper):
    """Generate the main Bayesian depth-confidence figure."""

    H = np.linspace(0, 4000, 500)

    # Parameters
    A_L, A_S = 0.4, 0.20
    H_c_L, H_c_S = 4000, 1000
    k_L, k_S = 0.005, 0.005

    # Compute curves
    P_hit_L = logistic_survival(H, A_L, H_c_L, k_L)
    P_hit_S = logistic_survival(H, A_S, H_c_S, k_S)
    P_post = np.array([posterior_probability(h, k_S=k_S) for h in H])

    # Sensitivity bands
    P_S_lo = logistic_survival(H, A_S, H_c_S, 0.008)
    P_S_hi = logistic_survival(H, A_S, H_c_S, 0.003)
    P_post_lo = np.array([posterior_probability(h, k_S=0.003) for h in H])
    P_post_hi = np.array([posterior_probability(h, k_S=0.008) for h in H])

    # ── Figure ──
    fig, ax1 = plt.subplots(figsize=(14, 8.5))
    c1, c2, c3 = '#d62728', '#1f77b4', '#2ca02c'

    # Left axis: survival functions
    ax1.set_xlabel('Present-Day Drill Depth $H$ (meters)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('$P(\\mathrm{Hit} \\mid \\mathrm{Type})$', fontsize=13)

    ax1.plot(H, P_hit_L, color=c3, linewidth=2.5,
             label=r'$P(\mathrm{Hit}|\mathrm{Large})$: $A_L$=0.4, $H_{c,L}$=4000 m')
    ax1.plot(H, P_hit_S, color=c1, linewidth=2.5, linestyle='--',
             label=r'$P(\mathrm{Hit}|\mathrm{Small})$: $A_S$=0.20, $H_{c,S}$=1000 m')
    ax1.fill_between(H, P_S_lo, P_S_hi, alpha=0.08, color=c1)
    ax1.tick_params(axis='y', labelsize=12)
    ax1.set_ylim(-0.02, 0.48)

    # H_c_S annotation
    ax1.axvline(x=1000, color=c1, linestyle=':', linewidth=1.2, alpha=0.4)
    ax1.annotate(r'$H_{c,S} = 1000$ m' + '\n(half-survival depth)',
                 xy=(1000, 0.10), xytext=(250, 0.03),
                 fontsize=10, color=c1, style='italic',
                 arrowprops=dict(arrowstyle='->', color=c1, lw=1.2))

    # Right axis: posterior
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'Posterior $P(\mathrm{Large} \mid \mathrm{Hit})$ (%)',
                   color=c2, fontsize=13, fontweight='bold')
    ax2.plot(H, P_post * 100, color=c2, linewidth=3.5,
             label=r'$P(\mathrm{Large}|\mathrm{Hit})$, central')
    ax2.fill_between(H, P_post_lo * 100, P_post_hi * 100, alpha=0.08, color=c2)
    ax2.tick_params(axis='y', labelcolor=c2, labelsize=12)
    ax2.set_ylim(-5, 105)

    # 2500 m marker
    target = 2500
    idx = (np.abs(H - target)).argmin()
    cert = P_post[idx] * 100
    ax2.axvline(x=target, color='#555', linestyle=':', linewidth=1.5, alpha=0.4)
    ax2.scatter([target], [cert], color='#FFD700', s=220,
                edgecolor='black', linewidths=1.8, zorder=10)
    ax2.annotate(f'$H = 2500$ m\n$P = {cert:.1f}\\%$',
                 xy=(target, cert), xytext=(3050, 68),
                 arrowprops=dict(facecolor='black', shrink=0.05, width=2, headwidth=9),
                 fontsize=12, fontweight='bold',
                 bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", ec="orange", lw=2))

    # Threshold markers
    for pct in [50, 90]:
        idx_p = np.argmin(np.abs(P_post * 100 - pct))
        d_p = H[idx_p]
        ax2.plot([d_p, d_p], [0, pct], 'b:', alpha=0.15, linewidth=0.8)
        ax2.text(d_p, -3, f'{pct}\\% @ {int(d_p)} m',
                 fontsize=8.5, color=c2, alpha=0.6, ha='center', va='top')

    # Formula box (lower right, clear area)
    ax1.text(3400, 0.13,
             'Logistic survival model:\n'
             r'$P = \frac{A}{1 + \exp\!\left(k(H - H_c)\right)}$' '\n'
             r'$H_c$ = critical closure depth' '\n'
             '(Sibson fault-valve)',
             fontsize=9.5, color='#333', ha='center', va='center',
             bbox=dict(boxstyle='round,pad=0.5', fc='white', ec='#999', alpha=1.0))

    # Legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc='upper right', fontsize=10, framealpha=0.95, edgecolor='gray',
               bbox_to_anchor=(0.98, 0.98))

    # Title
    plt.title(r'Bayesian Depth--Confidence Model with Logistic Survival Functions'
              '\n' r'(Prior: $P_0(\mathrm{Large}) = 0.01$ for targeted exploration site)',
              fontsize=14, fontweight='bold', pad=15)
    ax1.set_xlim(0, 4000)
    plt.tight_layout()

    # Save outputs
    os.makedirs(output_dir_figures, exist_ok=True)
    os.makedirs(output_dir_paper, exist_ok=True)

    png_path = os.path.join(output_dir_figures, "fig1_logistic_model.png")
    pdf_path = os.path.join(output_dir_paper, "fig1_latex.pdf")

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.close()

    print(f"  PNG: {png_path}")
    print(f"  PDF: {pdf_path}")


if __name__ == "__main__":
    repo_root = os.path.join(os.path.dirname(__file__), "..")
    print("Generating Figure 1...")
    generate_figure(
        output_dir_figures=os.path.join(repo_root, "figures"),
        output_dir_paper=os.path.join(repo_root, "paper")
    )
    print("Done.")
