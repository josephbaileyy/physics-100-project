"""Make a stacked AM CVn model-comparison plot.

The figure overlays four specific models on the same differential light curve:
single sine, fixed harmonic superhump template, beat-constrained precessing disk,
and the free double wave.  The output table reports the fitted parameters and
BIC/chi-square values used for comparing the models.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from amcvn.physical_model_suite import (
    ORBITAL_S,
    PER_DIR,
    PLOT_DIR,
    SUPERHUMP_S,
    design_matrix,
    fit_linear,
    load_light_curve,
    scan_beat_constrained_superhump,
    scan_period,
)


def amp_phase_from_beta(beta: np.ndarray, offset: int) -> tuple[float, float]:
    """Return semiamplitude in percent and phase in cycles for one sine term."""
    amp_pct = 100 * float(np.hypot(beta[offset], beta[offset + 1]))
    phase_cycles = float((math.atan2(beta[offset + 1], beta[offset]) / (2 * np.pi)) % 1)
    return amp_pct, phase_cycles


def summarize_model(
    model_id: str,
    display_name: str,
    interpretation: str,
    components: list[tuple[float, int]],
    fit: dict,
    free_double_bic: float,
    free_double_effective_bic: float,
    error_scale: float,
    nonlinear_parameters: int,
    precession_period_hr: float | None = None,
) -> dict:
    beta = fit["beta"]
    fundamental_amp_pct = np.nan
    fundamental_phase = np.nan
    harmonic_amp_pct = np.nan
    harmonic_phase = np.nan
    if components:
        fundamental_amp_pct, fundamental_phase = amp_phase_from_beta(beta, 1)
    if len(components) > 1:
        harmonic_amp_pct, harmonic_phase = amp_phase_from_beta(beta, 3)

    scaled_chi2 = float(fit["chi2"]) / (error_scale**2)
    effective_k = int(fit["k"]) + int(nonlinear_parameters)
    effective_bic = float(fit["chi2"]) + effective_k * math.log(int(fit["N"]))
    return {
        "model": model_id,
        "display_name": display_name,
        "physical_interpretation": interpretation,
        "period_s": float(60 * components[0][0]) if components else np.nan,
        "orbital_period_s": ORBITAL_S if model_id == "beat_constrained_precessing_disk" else np.nan,
        "precession_period_hr": precession_period_hr if precession_period_hr is not None else np.nan,
        "components": ";".join(f"{60 * period_min:.4f}s_h{harmonic}" for period_min, harmonic in components),
        "mean_differential_flux": float(beta[0]),
        "fundamental_amp_pct": fundamental_amp_pct,
        "fundamental_phase_cycles": fundamental_phase,
        "harmonic_amp_pct": harmonic_amp_pct,
        "harmonic_phase_cycles": harmonic_phase,
        "N": int(fit["N"]),
        "k": int(fit["k"]),
        "nonlinear_parameters": int(nonlinear_parameters),
        "effective_k": effective_k,
        "dof": int(fit["dof"]),
        "chi2": float(fit["chi2"]),
        "reduced_chi2": float(fit["reduced_chi2"]),
        "chi2_p_value_formal_errors": float(fit["chi2_p_value"]),
        "scaled_chi2": scaled_chi2,
        "scaled_reduced_chi2": scaled_chi2 / float(fit["dof"]),
        "scaled_chi2_p_value": float(stats.chi2.sf(scaled_chi2, int(fit["dof"]))),
        "bic_chi2": float(fit["bic_chi2"]),
        "delta_bic_vs_free_double_wave": float(fit["bic_chi2"] - free_double_bic),
        "bic_chi2_effective": effective_bic,
        "delta_effective_bic_vs_free_double_wave": float(effective_bic - free_double_effective_bic),
        "rms_pct": float(fit["rms_pct"]),
    }


def make_stacked_plot(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma_plot: np.ndarray,
    model_specs: list[dict],
    stats: pd.DataFrame,
) -> Path:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PLOT_DIR / "AM_CVn_stacked_model_fits.png"

    grid = np.linspace(float(t_min.min()), float(t_min.max()), 1200)
    y_pct = 100 * (y - 1)
    yerr_pct = 100 * sigma_plot
    ylim_pad = 0.6
    ymin = float(np.nanmin(y_pct - yerr_pct) - ylim_pad)
    ymax = float(np.nanmax(y_pct + yerr_pct) + ylim_pad)

    fig, axes = plt.subplots(
        len(model_specs),
        1,
        figsize=(11.5, 12.5),
        sharex=True,
        constrained_layout=True,
    )
    colors = ["0.2", "tab:blue", "tab:green", "tab:orange", "tab:red"]

    for ax, spec, color in zip(axes, model_specs, colors):
        row = stats.loc[stats.model.eq(spec["model"])].iloc[0]
        model_grid = design_matrix(grid, spec["components"]) @ spec["fit"]["beta"]
        ax.errorbar(
            t_min,
            y_pct,
            yerr=yerr_pct,
            fmt=".",
            ms=4,
            color="0.25",
            ecolor="0.82",
            elinewidth=0.6,
            capsize=0,
            alpha=0.9,
        )
        ax.plot(grid, 100 * (model_grid - 1), color=color, lw=2.0)
        ax.set_ylim(ymin, ymax)
        ax.set_ylabel("Flux - 1 [%]")
        period_text = f"P={row.period_s:.2f} s, " if np.isfinite(row.period_s) else ""
        ax.set_title(
            f"{row.display_name}: {period_text}"
            f"effective BIC={row.bic_chi2_effective:.2f}, "
            f"Delta={row.delta_effective_bic_vs_free_double_wave:.2f}",
            fontsize=11,
        )
        amp_text = ""
        if np.isfinite(row.fundamental_amp_pct):
            amp_text = f"A1={row.fundamental_amp_pct:.2f}%"
        if np.isfinite(row.harmonic_amp_pct):
            amp_text += f", A2={row.harmonic_amp_pct:.2f}%"
        if np.isfinite(row.precession_period_hr):
            amp_text += f"\nPprec={row.precession_period_hr:.2f} hr"
        if amp_text:
            amp_text += "\n"
        amp_text += (
            f"chi2/dof={row.chi2:.1f}/{int(row.dof)}\n"
            f"RMS={row.rms_pct:.2f}%"
        )
        ax.text(
            0.985,
            0.92,
            amp_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "0.8", "alpha": 0.9},
        )

    axes[-1].set_xlabel("Minutes from first exposure")
    fig.suptitle("AM CVn light curve with model overlays", fontsize=14)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    return out_path


def main() -> None:
    PER_DIR.mkdir(parents=True, exist_ok=True)
    good, t_min, y, sigma = load_light_curve()

    single_period_min, _, single_fit = scan_period(t_min, y, sigma, 5, 12, 3000, (1,))
    double_period_min, _, double_fit = scan_period(t_min, y, sigma, 14, 22, 4000, (1, 2))
    constant_fit = fit_linear(t_min, y, sigma, [])
    superhump_period_min = SUPERHUMP_S / 60
    harmonic_fit = fit_linear(t_min, y, sigma, [(superhump_period_min, 1), (superhump_period_min, 2)])
    beat_curve, beat_fit = scan_beat_constrained_superhump(t_min, y, sigma)
    beat_best = beat_curve.loc[int(beat_curve.chi2.idxmin())]
    beat_period_min = float(beat_best.superhump_period_min)

    model_specs = [
        {
            "model": "constant_flux",
            "display_name": "Constant flux baseline",
            "interpretation": "no coherent AM CVn variability",
            "components": [],
            "fit": constant_fit,
            "nonlinear_parameters": 0,
            "precession_period_hr": None,
        },
        {
            "model": "single_sine_free",
            "display_name": "Best-fit single sine",
            "interpretation": "free sinusoidal photometric period, dominated by the 525 s harmonic",
            "components": [(single_period_min, 1)],
            "fit": single_fit,
            "nonlinear_parameters": 1,
            "precession_period_hr": None,
        },
        {
            "model": "harmonic_superhump_template",
            "display_name": "Fixed harmonic superhump template",
            "interpretation": "literature 1051.2 s positive superhump plus first harmonic",
            "components": [(superhump_period_min, 1), (superhump_period_min, 2)],
            "fit": harmonic_fit,
            "nonlinear_parameters": 0,
            "precession_period_hr": None,
        },
        {
            "model": "beat_constrained_precessing_disk",
            "display_name": "Best beat-constrained precessing disk",
            "interpretation": "positive superhump with f_sh = f_orb - f_prec",
            "components": [(beat_period_min, 1), (beat_period_min, 2)],
            "fit": beat_fit,
            "nonlinear_parameters": 1,
            "precession_period_hr": float(beat_best.precession_period_hr),
        },
        {
            "model": "double_wave_free",
            "display_name": "Best free double wave",
            "interpretation": "free fundamental plus first harmonic waveform",
            "components": [(double_period_min, 1), (double_period_min, 2)],
            "fit": double_fit,
            "nonlinear_parameters": 1,
            "precession_period_hr": None,
        },
    ]

    free_double_bic = float(double_fit["bic_chi2"])
    free_double_effective_bic = float(double_fit["chi2"] + (int(double_fit["k"]) + 1) * math.log(int(double_fit["N"])))
    error_scale = math.sqrt(float(double_fit["reduced_chi2"]))
    rows = [
        summarize_model(
            spec["model"],
            spec["display_name"],
            spec["interpretation"],
            spec["components"],
            spec["fit"],
            free_double_bic,
            free_double_effective_bic,
            error_scale,
            spec["nonlinear_parameters"],
            spec["precession_period_hr"],
        )
        for spec in model_specs
    ]
    stats = pd.DataFrame(rows)
    stats_path = PER_DIR / "physical_stacked_model_fit_parameters.csv"
    stats.to_csv(stats_path, index=False)

    plot_path = make_stacked_plot(t_min, y, sigma * error_scale, model_specs, stats)

    print(f"Wrote {plot_path}")
    print(f"Wrote {stats_path}")
    print(
        stats[
            [
                "display_name",
                "period_s",
                "precession_period_hr",
                "fundamental_amp_pct",
                "harmonic_amp_pct",
                "chi2",
                "reduced_chi2",
                "chi2_p_value_formal_errors",
                "scaled_chi2_p_value",
                "bic_chi2",
                "delta_bic_vs_free_double_wave",
                "bic_chi2_effective",
                "delta_effective_bic_vs_free_double_wave",
                "rms_pct",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
