"""Create final presentation figures for the V0399 analysis."""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis/project/contamination_tests/V0399_UMa"


def main() -> None:
    plot_dir = OUT / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    lc = pd.read_csv(
        ROOT
        / "analysis/project/calibration/flat_candidate_photometry_qc/flat_4s_masked_clipmean_common8_light_curve.csv"
    )
    aper = pd.read_csv(OUT / "aperture_annulus_grid_summary.csv")
    two = pd.read_csv(OUT / "two_star_contamination_model_summary.csv")
    inj = pd.read_csv(OUT / "injection_recovery_summary.csv")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    ax = axes[0, 0]
    good = lc[lc.good.astype(bool)].copy()
    for seq, g in good.groupby("sequence"):
        ax.plot(g.minutes_combined, 100 * (g.diff_flux_norm - 1), ".", ms=4, label=f"V0399 {seq}")
    ax.plot(good.minutes_combined, 100 * (good.check_diff_flux_norm - 1), ".", ms=3, alpha=0.35, label="DN/check")
    ax.axhline(0, color="0.45", lw=0.8)
    ax.set_xlabel("Minutes from first exposure")
    ax.set_ylabel("Differential flux residual [%]")
    ax.set_title("Preferred common-8 light curve")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)

    ax = axes[0, 1]
    by_radius = aper.groupby("radius").median(numeric_only=True).reset_index()
    ax.plot(by_radius.radius, by_radius.V0399_scatter_pct, "o-", label="V0399 scatter")
    ax.plot(by_radius.radius, by_radius.DN_scatter_pct, "o-", label="DN scatter")
    ax2 = ax.twinx()
    ax2.plot(by_radius.radius, by_radius.V0399_DN_corr_r, "s--", color="black", label="V0399-DN corr.")
    ax.set_xlabel("Aperture radius [px]")
    ax.set_ylabel("Robust scatter [%]")
    ax2.set_ylabel("Correlation r")
    ax.set_title("Aperture sensitivity")
    ax.grid(alpha=0.25)
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, fontsize=8, loc="best")

    ax = axes[1, 0]
    ordered = two.sort_values("bic").copy()
    ordered["bic_delta"] = ordered.bic - ordered.bic.min()
    labels = ordered.model.str.replace("v0399_", "", regex=False).str.replace("_", " ", regex=False)
    y = np.arange(len(ordered))
    ax.barh(y, ordered.bic_delta, color="tab:blue")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Delta BIC from best")
    ax.set_title("Two-star/common-mode model comparison")
    ax.grid(axis="x", alpha=0.25)

    ax = axes[1, 1]
    pivot = inj[inj.waveform == "sine"].pivot_table(
        index="control", columns="injected_amp_pct", values="fit_15p3_bic_delta", aggfunc="median"
    )
    im = ax.imshow(pivot.to_numpy(float), aspect="auto", cmap="coolwarm_r")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels([f"{c:.1f}%" for c in pivot.columns])
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_xlabel("Injected 15.3 min sine amplitude")
    ax.set_title("Injection recovery; negative favors signal")
    fig.colorbar(im, ax=ax, label="BIC(sine) - BIC(constant)")

    fig.suptitle("V0399 UMa robustness summary", fontsize=14)
    fig.savefig(plot_dir / "V0399_final_summary_panel.png", dpi=170)
    plt.close(fig)
    print(plot_dir / "V0399_final_summary_panel.png")


if __name__ == "__main__":
    main()
