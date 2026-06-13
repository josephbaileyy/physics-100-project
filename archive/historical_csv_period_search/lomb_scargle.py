"""Run a Lomb-Scargle periodogram on an ASAS-SN photometry CSV.

Usage:
    python lomb_scargle.py "data/ASASSN-V_J093023.98+612034.3.csv"
    python lomb_scargle.py "data/ASASSN-V_J123454.75+373746.7.csv"
"""
import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle


def load_photometry_timeseries(csv_path, quantity):
    """Load every valid measurement as a continuous photometry time series."""
    df = pd.read_csv(csv_path, low_memory=False)

    if {"jd", "diff_flux", "diff_flux_err"}.issubset(df.columns):
        time_col = "jd"
        value_col = "diff_flux"
        error_col = "diff_flux_err"
    elif {"hjd", "mag", "mag_err"}.issubset(df.columns):
        time_col = "hjd"
        if quantity == "mag":
            value_col = "mag"
            error_col = "mag_err"
        elif quantity == "flux":
            value_col = "flux"
            error_col = "flux_err"
        else:
            raise ValueError(f"unsupported quantity: {quantity}")
    elif {"JD", "Magnitude"}.issubset(df.columns):
        time_col = "JD"
        value_col = "Magnitude"
        error_col = "Uncertainty" if "Uncertainty" in df.columns else None
    else:
        raise ValueError(
            f"{csv_path} must contain May 20 columns jd/diff_flux/diff_flux_err, "
            "ASAS-SN columns hjd/mag/mag_err, or AAVSO columns JD/Magnitude"
        )

    required = {time_col, value_col}
    if error_col:
        required.add(error_col)
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{csv_path} is missing required columns: {missing}")

    t = pd.to_numeric(df[time_col], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[value_col], errors="coerce").to_numpy(dtype=float)
    if error_col:
        dy = pd.to_numeric(df[error_col], errors="coerce").to_numpy(dtype=float)
    else:
        dy = np.full(len(t), np.nan)

    good = np.isfinite(t) & np.isfinite(y)
    t, y, dy = t[good], y[good], dy[good]

    finite_errors = np.isfinite(dy) & (dy > 0)
    if finite_errors.any():
        fill_error = np.nanmedian(dy[finite_errors])
    else:
        scatter = np.nanmedian(np.abs(y - np.nanmedian(y)))
        fill_error = scatter if np.isfinite(scatter) and scatter > 0 else 0.03
    dy = np.where(np.isfinite(dy) & (dy > 0), dy, fill_error)

    # Keep the full observing timeline for Lomb-Scargle; sorting changes only
    # the row order, not the timestamps or their multi-period spacing.
    order = np.argsort(t)
    return t[order], y[order], dy[order], value_col, time_col


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv", type=Path)
    p.add_argument("--min-period", type=float, default=0.02,
                   help="minimum period in days (default 0.02 d ≈ 29 min)")
    p.add_argument("--max-period", type=float, default=1.0,
                   help="maximum period in days (default 1 d)")
    p.add_argument("--samples-per-peak", type=int, default=10)
    p.add_argument("--quantity", choices=("mag", "flux"), default="mag",
                   help="ASAS-SN photometry column to analyze (default: mag)")
    p.add_argument("--time-min", type=float,
                   help="keep only observations at or after this time")
    p.add_argument("--time-max", type=float,
                   help="keep only observations at or before this time")
    p.add_argument("--detrend", choices=("none", "linear"), default="none",
                   help="remove a trend before Lomb-Scargle (default: none)")
    p.add_argument("--fold-period", type=float,
                   help="also write a folded plot at this period in days")
    p.add_argument("--outdir", type=Path, default=Path("results"))
    args = p.parse_args()

    t_abs, y, dy, value_col, time_col = load_photometry_timeseries(
        args.csv, args.quantity
    )
    keep = np.ones(len(t_abs), dtype=bool)
    if args.time_min is not None:
        keep &= t_abs >= args.time_min
    if args.time_max is not None:
        keep &= t_abs <= args.time_max
    t_abs, y, dy = t_abs[keep], y[keep], dy[keep]

    if len(t_abs) < 3:
        raise ValueError("Need at least 3 valid observations for Lomb-Scargle")

    # Subtracting a constant improves numerical conditioning while preserving
    # all gaps between nights/seasons. This is not phase folding or binning.
    t = t_abs - t_abs.min()
    y_for_ls = y.copy()
    if args.detrend == "linear":
        if len(t) < 2:
            raise ValueError("Need at least 2 observations for linear detrending")
        trend = np.polyval(np.polyfit(t, y_for_ls, deg=1), t)
        y_for_ls = y_for_ls - trend + np.nanmedian(y_for_ls)

    print(f"Loaded {len(t)} good points spanning {t.max() - t.min():.1f} d "
          f"(median cadence {np.median(np.diff(np.sort(t))):.3f} d)")
    print(f"Lomb-Scargle input uses the full continuous {time_col} timeline; "
          "phase folding happens only after the best period is found.")
    if args.detrend != "none":
        print(f"Applied {args.detrend} detrending before Lomb-Scargle.")

    # Convert period bounds to frequency bounds for autopower.
    f_min = 1.0 / args.max_period
    f_max = 1.0 / args.min_period

    ls = LombScargle(t, y_for_ls, dy)
    frequency, power = ls.autopower(
        minimum_frequency=f_min,
        maximum_frequency=f_max,
        samples_per_peak=args.samples_per_peak,
    )

    best_idx = np.argmax(power)
    best_freq = frequency[best_idx]
    best_period = 1.0 / best_freq
    best_power = power[best_idx]

    fap = ls.false_alarm_probability(best_power)

    print(f"Best period: {best_period:.6f} d  ({best_period * 24 * 60:.2f} min)")
    print(f"Best frequency: {best_freq:.4f} cycles/day")
    print(f"Peak power: {best_power:.4f}")
    print(f"False-alarm probability: {fap:.3e}")

    args.outdir.mkdir(exist_ok=True)
    stem = args.csv.stem.replace(" ", "_")

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.errorbar(t_abs, y, yerr=dy, fmt=".", ms=3, lw=0.5, alpha=0.75, color="C0")
    if value_col == "mag":
        ax.invert_yaxis()  # brighter up
    ax.set_xlabel(time_col)
    ax.set_ylabel(value_col)
    ax.set_title(f"{args.csv.stem} input time series")
    fig.tight_layout()
    timeseries_path = args.outdir / f"{stem}_timeseries.png"
    fig.savefig(timeseries_path, dpi=150)
    print(f"Wrote {timeseries_path}")

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(1.0 / frequency, power, lw=0.7)
    ax.axvline(best_period, color="C3", lw=1, ls="--",
               label=f"P = {best_period:.5f} d")
    ax.set_xscale("log")
    ax.set_xlabel("Period (days)")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title(f"{args.csv.stem} — Lomb-Scargle periodogram")
    ax.legend()
    fig.tight_layout()
    pgram_path = args.outdir / f"{stem}_periodogram.png"
    fig.savefig(pgram_path, dpi=150)
    print(f"Wrote {pgram_path}")

    phase = ((t - t.min()) * best_freq) % 1.0
    order = np.argsort(phase)

    fig, ax = plt.subplots(figsize=(9, 4))
    # Show two cycles for readability.
    ax.errorbar(phase[order], y_for_ls[order], yerr=dy[order],
                fmt=".", ms=3, lw=0.5, alpha=0.7, color="C0")
    ax.errorbar(phase[order] + 1, y_for_ls[order], yerr=dy[order],
                fmt=".", ms=3, lw=0.5, alpha=0.7, color="C0")
    if value_col == "mag":
        ax.invert_yaxis()  # brighter up
    ax.set_xlabel("Phase")
    ax.set_ylabel(value_col)
    ax.set_title(f"{args.csv.stem} folded at P = {best_period:.5f} d")
    fig.tight_layout()
    fold_path = args.outdir / f"{stem}_folded.png"
    fig.savefig(fold_path, dpi=150)
    print(f"Wrote {fold_path}")

    if args.fold_period:
        fold_frequency = 1.0 / args.fold_period
        phase = (t * fold_frequency) % 1.0
        order = np.argsort(phase)

        fig, ax = plt.subplots(figsize=(9, 4))
        ax.errorbar(phase[order], y_for_ls[order], yerr=dy[order],
                    fmt=".", ms=3, lw=0.5, alpha=0.7, color="C0")
        ax.errorbar(phase[order] + 1, y_for_ls[order], yerr=dy[order],
                    fmt=".", ms=3, lw=0.5, alpha=0.7, color="C0")
        if value_col == "mag":
            ax.invert_yaxis()  # brighter up
        ax.set_xlabel("Phase")
        ax.set_ylabel(value_col)
        ax.set_title(f"{args.csv.stem} folded at P = {args.fold_period:.6f} d")
        fig.tight_layout()
        fixed_fold_path = args.outdir / f"{stem}_folded_{args.fold_period:.6f}d.png"
        fig.savefig(fixed_fold_path, dpi=150)
        print(f"Wrote {fixed_fold_path}")


if __name__ == "__main__":
    sys.exit(main())
