"""Run a Lomb-Scargle periodogram on an ASAS-SN photometry CSV.

Usage:
    python lomb_scargle.py "data/ASASSN-V J093023.98+612034.3.csv"
"""
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv", type=Path)
    p.add_argument("--min-period", type=float, default=0.02,
                   help="minimum period in days (default 0.02 d ≈ 29 min)")
    p.add_argument("--max-period", type=float, default=1.0,
                   help="maximum period in days (default 1 d)")
    p.add_argument("--samples-per-peak", type=int, default=10)
    p.add_argument("--outdir", type=Path, default=Path("results"))
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    # ASAS-SN columns: hjd, camera, mag, mag_err, flux, flux_err
    t = df["hjd"].to_numpy()
    y = df["mag"].to_numpy()
    dy = df["mag_err"].to_numpy()

    good = np.isfinite(t) & np.isfinite(y) & np.isfinite(dy) & (dy > 0)
    t, y, dy = t[good], y[good], dy[good]
    print(f"Loaded {len(t)} good points spanning {t.max() - t.min():.1f} d "
          f"(median cadence {np.median(np.diff(np.sort(t))):.3f} d)")

    # Convert period bounds to frequency bounds for autopower.
    f_min = 1.0 / args.max_period
    f_max = 1.0 / args.min_period

    ls = LombScargle(t, y, dy)
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
    ax.errorbar(phase[order], y[order], yerr=dy[order],
                fmt=".", ms=3, lw=0.5, alpha=0.7, color="C0")
    ax.errorbar(phase[order] + 1, y[order], yerr=dy[order],
                fmt=".", ms=3, lw=0.5, alpha=0.7, color="C0")
    ax.invert_yaxis()  # brighter up
    ax.set_xlabel("Phase")
    ax.set_ylabel("mag")
    ax.set_title(f"{args.csv.stem} folded at P = {best_period:.5f} d")
    fig.tight_layout()
    fold_path = args.outdir / f"{stem}_folded.png"
    fig.savefig(fold_path, dpi=150)
    print(f"Wrote {fold_path}")


if __name__ == "__main__":
    sys.exit(main())
