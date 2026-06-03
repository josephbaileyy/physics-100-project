"""AM CVn AAVSO-sequence differential photometry and period analysis."""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy import stats
import astropy.units as u

from amcvn.pipeline import centroid_peak_in_box, detect_sources, load_config, robust_sigma, source_table


ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "observations/targets.json"
CHART = ROOT / "analysis/am_cvn/aavso/X42421BZ.json"
OUT = ROOT / "analysis/am_cvn/photometry/sequence_X42421BZ"
PER_OUT = ROOT / "analysis/am_cvn/periodograms/sequence_X42421BZ"
PLOT_OUT = ROOT / "analysis/am_cvn/plots"
FLAT = "flat_4s_masked_clipmean"
SATURATION = 58000.0
APERTURE_R = 6.0
ANN_IN = 12.0
ANN_OUT = 20.0
GOOD_REFERENCE_DIR = ROOT / "analysis/am_cvn/qc/good_reference"
GOOD_REFERENCE_SOURCES = GOOD_REFERENCE_DIR / "good_reference_sources.csv"
GOOD_REFERENCE_INFO = GOOD_REFERENCE_DIR / "good_reference_info.csv"
ANCHOR_RANKS = {"90": 1, "96": 2, "112": 5, "125": 9, "target": 16, "143": 17, "144": 18}
COMPARISON_LABELS = ["90", "96", "112", "125", "143", "144"]


def ensure_dirs() -> None:
    for path in [OUT, PER_OUT, PLOT_OUT, ROOT / "analysis/am_cvn/aavso", GOOD_REFERENCE_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def aperture_flux(data: np.ndarray, x: float, y: float, r: float, rin: float, rout: float) -> dict:
    height, width = data.shape
    pad = int(math.ceil(rout)) + 2
    if x < pad or y < pad or x > width - pad or y > height - pad:
        return {"flux": np.nan, "flux_err": np.nan, "peak": np.nan, "flag": "edge"}
    x0 = int(round(x)) - pad
    x1 = int(round(x)) + pad + 1
    y0 = int(round(y)) - pad
    y1 = int(round(y)) + pad + 1
    cut = data[y0:y1, x0:x1].astype(float)
    yy, xx = np.indices(cut.shape)
    xx = xx + x0
    yy = yy + y0
    rr2 = (xx - x) ** 2 + (yy - y) ** 2
    aper = rr2 <= r**2
    ann = (rr2 >= rin**2) & (rr2 <= rout**2)
    if ann.sum() < 20:
        return {"flux": np.nan, "flux_err": np.nan, "peak": np.nan, "flag": "bad_annulus"}
    bkg = np.nanmedian(cut[ann])
    bkg_std = np.nanstd(cut[ann])
    values = cut[aper]
    flux = float(np.nansum(values - bkg))
    err = float(math.sqrt(max(abs(flux), 0.0) + aper.sum() * bkg_std**2))
    peak = float(np.nanmax(values))
    flag = "ok" if np.isfinite(flux) and flux > 0 else "nonpositive_flux"
    return {"flux": flux, "flux_err": err, "peak": peak, "flag": flag}


def sky_offsets(chart: dict) -> pd.DataFrame:
    center = SkyCoord(chart["ra"], chart["dec"], unit=(u.hourangle, u.deg))
    rows = [{"label": "target", "role": "target", "ra": chart["ra"], "dec": chart["dec"], "v_mag": np.nan}]
    for row in chart["photometry"]:
        vmag = next((b["mag"] for b in row["bands"] if b["band"] == "V"), np.nan)
        rows.append(
            {
                "label": str(row["label"]),
                "role": "comparison",
                "ra": row["ra"],
                "dec": row["dec"],
                "v_mag": vmag,
                "auid": row.get("auid", ""),
            }
        )
    out = []
    for row in rows:
        coord = SkyCoord(row["ra"], row["dec"], unit=(u.hourangle, u.deg))
        dra = (coord.ra.deg - center.ra.deg) * math.cos(math.radians(center.dec.deg)) * 3600
        ddec = (coord.dec.deg - center.dec.deg) * 3600
        out.append({**row, "dra_arcsec": dra, "ddec_arcsec": ddec})
    return pd.DataFrame(out)


def reference_positions(seq: pd.DataFrame) -> pd.DataFrame:
    sources = good_reference_sources()
    anchors = []
    for label, rank in ANCHOR_RANKS.items():
        sky = seq[seq.label == label].iloc[0]
        src = sources.iloc[rank - 1]
        anchors.append((sky.dra_arcsec, sky.ddec_arcsec, src.x, src.y))
    anchors = np.asarray(anchors, dtype=float)
    a = np.column_stack([anchors[:, 0], anchors[:, 1], np.ones(len(anchors))])
    coef_x = np.linalg.lstsq(a, anchors[:, 2], rcond=None)[0]
    coef_y = np.linalg.lstsq(a, anchors[:, 3], rcond=None)[0]
    rows = []
    for row in seq.itertuples(index=False):
        vec = np.array([row.dra_arcsec, row.ddec_arcsec, 1.0])
        x = float(vec @ coef_x)
        y = float(vec @ coef_y)
        rows.append({**row._asdict(), "ref_x": x, "ref_y": y})
    ref = pd.DataFrame(rows)
    # Only use the sequence stars that are detected reliably in the reference.
    ref["used_for_photometry"] = ref.label.isin(["target", *COMPARISON_LABELS])
    ref["used_as_comparison"] = ref.label.isin(COMPARISON_LABELS)
    return ref


def good_reference_sources() -> pd.DataFrame:
    if GOOD_REFERENCE_SOURCES.exists() and GOOD_REFERENCE_INFO.exists():
        return pd.read_csv(GOOD_REFERENCE_SOURCES)
    manifest = pd.read_csv(
        ROOT / "analysis/am_cvn/reduced/AM_CVn" / f"reduction_manifest_{FLAT}.csv"
    ).sort_values("jd")
    manifest["minutes"] = (manifest.jd - manifest.jd.min()) * 24 * 60
    candidates = manifest[(manifest.minutes.between(45, 75)) & (manifest.max_adu_raw.between(30000, 56000))]
    if candidates.empty:
        candidates = manifest[manifest.max_adu_raw.between(30000, 56000)]
    row = candidates.iloc[len(candidates) // 2]
    path = ROOT / "analysis/am_cvn/reduced/AM_CVn" / row.reduced_file
    data = fits.getdata(path).astype(float)
    sources = source_table(detect_sources(data, threshold_sigma=5.0), limit=120)
    sources.to_csv(GOOD_REFERENCE_SOURCES, index=False)
    pd.DataFrame([row]).to_csv(GOOD_REFERENCE_INFO, index=False)
    return sources


def raw_path_for(cfg: dict, raw_file: str) -> Path:
    return Path(cfg["light_dir"]) / raw_file


def frame_photometry(cfg: dict, ref: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    manifest = pd.read_csv(
        ROOT / "analysis/am_cvn/reduced/AM_CVn" / f"reduction_manifest_{FLAT}.csv"
    ).sort_values("jd")
    reduced_dir = ROOT / "analysis/am_cvn/reduced/AM_CVn"
    objects = ref[ref.used_for_photometry].copy()
    rows = []
    for idx, frame in enumerate(manifest.itertuples(index=False), start=1):
        cal = fits.getdata(reduced_dir / frame.reduced_file).astype(float)
        raw = fits.getdata(raw_path_for(cfg, frame.raw_file)).astype(float)
        shifts = []
        centroids = {}
        for obj in objects[objects.used_as_comparison].itertuples(index=False):
            ok, cx, cy = centroid_peak_in_box(cal, float(obj.ref_x), float(obj.ref_y), half_size=35, aperture=10)
            if ok:
                centroids[obj.label] = (cx, cy)
                shifts.append((cx - obj.ref_x, cy - obj.ref_y))
        if len(shifts) >= 2:
            shift = np.asarray(shifts, dtype=float)
            dx = float(np.nanmedian(shift[:, 0]))
            dy = float(np.nanmedian(shift[:, 1]))
            scatter = float(np.nanmedian(np.hypot(shift[:, 0] - dx, shift[:, 1] - dy)))
            registered = scatter < 8 and abs(dx) < 80 and abs(dy) < 80
        else:
            dx = dy = scatter = np.nan
            registered = False
        frame_artifact = float(frame.max_adu_raw) < 10000
        object_results = []
        for obj in objects.itertuples(index=False):
            if obj.label in centroids:
                x, y = centroids[obj.label]
            else:
                x, y = obj.ref_x + (dx if registered else 0.0), obj.ref_y + (dy if registered else 0.0)
            phot = aperture_flux(cal, x, y, APERTURE_R, ANN_IN, ANN_OUT)
            raw_phot = aperture_flux(raw, x, y, APERTURE_R, ANN_IN, ANN_OUT)
            saturated = bool(np.isfinite(raw_phot["peak"]) and raw_phot["peak"] >= SATURATION)
            flag = phot["flag"]
            if frame_artifact:
                flag = "frame_artifact"
            elif not registered:
                flag = "registration_failed"
            elif saturated:
                flag = "saturated"
            flux_rate = phot["flux"] / frame.exptime if np.isfinite(phot["flux"]) else np.nan
            flux_rate_err = phot["flux_err"] / frame.exptime if np.isfinite(phot["flux_err"]) else np.nan
            result = {
                "file": frame.reduced_file,
                "raw_file": frame.raw_file,
                "jd": frame.jd,
                "exptime": frame.exptime,
                "object_id": obj.label,
                "role": obj.role,
                "v_mag": obj.v_mag,
                "x": x,
                "y": y,
                "registered": registered,
                "dx": dx,
                "dy": dy,
                "shift_scatter": scatter,
                "n_anchor_centroids": len(shifts),
                "flux": phot["flux"],
                "flux_err": phot["flux_err"],
                "flux_rate": flux_rate,
                "flux_rate_err": flux_rate_err,
                "raw_peak": raw_phot["peak"],
                "saturated": saturated,
                "frame_artifact": frame_artifact,
                "flag": flag,
            }
            rows.append(result)
            object_results.append(result)
    phot = pd.DataFrame(rows)
    medians = (
        phot[(phot.flag == "ok") & np.isfinite(phot.flux_rate) & (phot.flux_rate > 0)]
        .groupby("object_id")["flux_rate"]
        .median()
        .to_dict()
    )
    diff_rows = []
    for frame, g in phot.groupby("file", sort=False):
        target = g[g.object_id == "target"].iloc[0]
        comps = g[
            g.object_id.isin(COMPARISON_LABELS)
            & (g.flag == "ok")
            & np.isfinite(g.flux_rate)
            & (g.flux_rate > 0)
        ].copy()
        comps["relative_flux"] = [row.flux_rate / medians.get(row.object_id, np.nan) for row in comps.itertuples()]
        comps = comps[np.isfinite(comps.relative_flux) & (comps.relative_flux > 0)]
        common_mode = float(np.nanmedian(comps.relative_flux)) if len(comps) >= 2 else np.nan
        target_median = medians.get("target", np.nan)
        target_relative = target.flux_rate / target_median if np.isfinite(target_median) and target_median > 0 else np.nan
        good = (
            target.flag == "ok"
            and np.isfinite(target_relative)
            and target_relative > 0
            and np.isfinite(common_mode)
            and common_mode > 0
            and not bool(target.frame_artifact)
        )
        diff_rows.append(
            {
                "file": frame,
                "jd": target.jd,
                "exptime": target.exptime,
                "diff_flux": target_relative / common_mode if good else np.nan,
                "diff_flux_err": target.flux_rate_err / target_median / common_mode if good else np.nan,
                "target_flux_rate": target.flux_rate,
                "comparison_ensemble_flux_rate": common_mode,
                "n_comparisons": len(comps),
                "registered": bool(target.registered),
                "frame_artifact": bool(target.frame_artifact),
                "good": good,
                "reject_reason": "ok" if good else ("frame_artifact" if bool(target.frame_artifact) else target.flag),
            }
        )
    diff = pd.DataFrame(diff_rows).sort_values("jd")
    diff["minutes"] = (diff.jd - diff.jd.min()) * 24 * 60
    ok = diff.good.astype(bool) & np.isfinite(diff.diff_flux)
    diff["diff_flux_norm"] = np.nan
    diff.loc[ok, "diff_flux_norm"] = diff.loc[ok, "diff_flux"] / np.nanmedian(diff.loc[ok, "diff_flux"])
    phot.to_csv(OUT / "aperture_photometry.csv", index=False)
    diff.to_csv(OUT / "differential_light_curve.csv", index=False)
    return phot, diff


def comparison_controls(phot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    series_rows = []
    medians = (
        phot[(phot.flag == "ok") & np.isfinite(phot.flux_rate) & (phot.flux_rate > 0)]
        .groupby("object_id")["flux_rate"]
        .median()
        .to_dict()
    )
    for target in COMPARISON_LABELS:
        points = []
        for (file, jd), g in phot.groupby(["file", "jd"]):
            targ = g[(g.object_id == target) & (g.flag == "ok")]
            comps = g[
                g.object_id.isin([c for c in COMPARISON_LABELS if c != target])
                & (g.flag == "ok")
                & (g.flux_rate > 0)
            ].copy()
            if len(targ) == 1 and len(comps) >= 2:
                comps["relative_flux"] = [
                    row.flux_rate / medians.get(row.object_id, np.nan) for row in comps.itertuples()
                ]
                common_mode = float(np.nanmedian(comps.relative_flux))
                targ_median = medians.get(target, np.nan)
                points.append(
                    {
                        "object_id": target,
                        "file": file,
                        "jd": jd,
                        "norm_flux": (targ.flux_rate.iloc[0] / targ_median) / common_mode,
                    }
                )
        df = pd.DataFrame(points).sort_values("jd")
        if df.empty:
            continue
        df["norm_flux"] = df.norm_flux / np.nanmedian(df.norm_flux)
        df["minutes"] = (df.jd - phot.jd.min()) * 24 * 60
        series_rows.extend(df.to_dict("records"))
        sigma = 1.4826 * np.nanmedian(np.abs(df.norm_flux - np.nanmedian(df.norm_flux)))
        chi2 = float(np.sum(((df.norm_flux - 1) / max(sigma, 1e-9)) ** 2))
        dof = max(len(df) - 1, 1)
        rows.append(
            {
                "object_id": target,
                "N": len(df),
                "robust_scatter_pct": 100 * sigma,
                "reduced_chi2_constant_empirical": chi2 / dof,
                "constant_p_empirical": float(stats.chi2.sf(chi2, dof)),
            }
        )
    out = pd.DataFrame(rows).sort_values("robust_scatter_pct")
    out.to_csv(OUT / "comparison_control_summary.csv", index=False)
    series = pd.DataFrame(series_rows)
    series.to_csv(OUT / "comparison_control_light_curves.csv", index=False)
    plot_comparison_controls(series, out)
    return out


def plot_comparison_controls(series: pd.DataFrame, summary: pd.DataFrame) -> None:
    if series.empty or summary.empty:
        return
    labels = summary.object_id.astype(str).tolist()
    ncols = 2
    nrows = int(math.ceil(len(labels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 8), sharex=True, sharey=True, constrained_layout=True)
    flat_axes = np.ravel(axes)
    for ax, label in zip(flat_axes, labels):
        sub = series[series.object_id.astype(str).eq(label)].sort_values("minutes")
        stats_row = summary[summary.object_id.astype(str).eq(label)].iloc[0]
        ax.plot(sub.minutes, 100 * (sub.norm_flux - 1), ".", ms=3.5, color="tab:blue")
        ax.axhline(0, color="0.25", lw=0.9)
        scatter = float(stats_row.robust_scatter_pct)
        ax.axhspan(-scatter, scatter, color="tab:blue", alpha=0.08, lw=0)
        ax.set_title(
            f"Comp {label}: scatter={scatter:.2f}%, "
            f"chi2_nu={float(stats_row.reduced_chi2_constant_empirical):.2f}",
            fontsize=10,
        )
        ax.grid(alpha=0.25)
    for ax in flat_axes[len(labels) :]:
        ax.axis("off")
    for ax in flat_axes[-ncols:]:
        ax.set_xlabel("Minutes from first exposure")
    for ax in flat_axes[::ncols]:
        ax.set_ylabel("Pseudo-differential flux - 1 [%]")
    fig.suptitle("Comparison-star constant-flux controls", fontsize=14)
    fig.savefig(PLOT_OUT / "AM_CVn_comparison_star_controls.png", dpi=180)
    plt.close(fig)


def fit_sine(data: pd.DataFrame, period_min: float) -> dict:
    t_min = (data.jd.to_numpy(float) - data.jd.min()) * 24 * 60
    y = data.diff_flux_norm.to_numpy(float)
    phase = 2 * np.pi * t_min / period_min
    x0 = np.ones((len(y), 1))
    x1 = np.column_stack([np.ones(len(y)), np.sin(phase), np.cos(phase)])
    b0 = np.linalg.lstsq(x0, y, rcond=None)[0]
    b1 = np.linalg.lstsq(x1, y, rcond=None)[0]
    r0 = y - x0 @ b0
    r1 = y - x1 @ b1
    sse0 = float(np.sum(r0**2))
    sse1 = float(np.sum(r1**2))
    n = len(y)
    bic0 = n * math.log(max(sse0 / n, 1e-18)) + x0.shape[1] * math.log(n)
    bic1 = n * math.log(max(sse1 / n, 1e-18)) + x1.shape[1] * math.log(n)
    return {
        "period_min": period_min,
        "amp_pct": 100 * float(np.hypot(b1[1], b1[2])),
        "delta_sse": sse0 - sse1,
        "bic_constant": bic0,
        "bic_sine": bic1,
        "bic_delta": bic1 - bic0,
        "rms_sine_pct": 100 * float(np.sqrt(np.mean(r1**2))),
    }


def period_analysis(diff: pd.DataFrame) -> pd.DataFrame:
    good = diff[diff.good.astype(bool) & np.isfinite(diff.diff_flux_norm)].sort_values("jd").copy()
    t = good.jd.to_numpy(float) - good.jd.min()
    y = good.diff_flux_norm.to_numpy(float)
    y0 = y - np.nanmedian(y)
    dy = np.full_like(y0, max(robust_sigma(y0), 1e-4))
    ls = LombScargle(t, y0, dy)
    rows = []
    for label, pmin, pmax in [("short_5_to_90min", 5, 90), ("wide_3_to_180min", 3, 180)]:
        freq, power = ls.autopower(
            minimum_frequency=1 / (pmax / (24 * 60)),
            maximum_frequency=1 / (pmin / (24 * 60)),
            samples_per_peak=20,
        )
        idx = int(np.nanargmax(power))
        period_min = float((1 / freq[idx]) * 24 * 60)
        fap = float(ls.false_alarm_probability(power[idx]))
        pd.DataFrame({"period_min": (1 / freq) * 24 * 60, "frequency_cpd": freq, "power": power}).to_csv(
            PER_OUT / f"{label}_periodogram.csv", index=False
        )
        rows.append(
            {
                "run": label,
                "N": len(good),
                "best_period_min": period_min,
                "best_power": float(power[idx]),
                "false_alarm_probability": fap,
                **{f"best_sine_{k}": v for k, v in fit_sine(good, period_min).items()},
            }
        )
        plot_periodogram(freq, power, period_min, label)
    fixed_periods = [17.14, 17.2, 34.3, 13.0, 20.0]
    fixed = pd.DataFrame([{"run": "fixed", **fit_sine(good, p)} for p in fixed_periods])
    fixed.to_csv(PER_OUT / "fixed_period_sine_summary.csv", index=False)
    summary = pd.DataFrame(rows)
    summary.to_csv(PER_OUT / "period_results.csv", index=False)
    plot_light_curve(good)
    plot_fold(good, float(summary.iloc[0].best_period_min), "best_short")
    plot_fold(good, 17.14, "fixed_17p14min")
    return summary


def plot_light_curve(good: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    ax.plot(good.minutes, 100 * (good.diff_flux_norm - 1), ".", ms=4)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlabel("Minutes from first AM CVn exposure")
    ax.set_ylabel("Differential flux residual [%]")
    ax.set_title("AM CVn differential light curve; AAVSO X42421BZ")
    ax.grid(alpha=0.25)
    fig.savefig(PLOT_OUT / "AM_CVn_sequence_X42421BZ_light_curve.png", dpi=160)
    plt.close(fig)


def plot_periodogram(freq: np.ndarray, power: np.ndarray, best_period: float, label: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True)
    ax.plot((1 / freq) * 24 * 60, power, color="black", lw=0.8)
    ax.axvline(best_period, color="tab:red", ls="--", label=f"{best_period:.2f} min")
    ax.set_xscale("log")
    ax.set_xlabel("Period [min]")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title(f"AM CVn periodogram: {label}")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(PER_OUT / f"{label}_periodogram.png", dpi=160)
    plt.close(fig)


def plot_fold(good: pd.DataFrame, period_min: float, label: str) -> None:
    phase = ((good.jd - good.jd.min()) * 24 * 60 / period_min) % 1
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    ax.plot(phase, good.diff_flux_norm, ".", ms=4)
    ax.plot(phase + 1, good.diff_flux_norm, ".", ms=4, color="C0")
    bins = np.linspace(0, 1, 13)
    mids = 0.5 * (bins[:-1] + bins[1:])
    means = [np.nanmean(good.diff_flux_norm[(phase >= bins[i]) & (phase < bins[i + 1])]) for i in range(len(mids))]
    ax.plot(mids, means, "ko-", ms=4)
    ax.plot(mids + 1, means, "ko-", ms=4)
    ax.set_xlabel("Phase")
    ax.set_ylabel("Normalized differential flux")
    ax.set_title(f"AM CVn folded at {period_min:.2f} min")
    ax.grid(alpha=0.25)
    fig.savefig(PER_OUT / f"{label}_folded.png", dpi=160)
    plt.close(fig)


def plot_reference_overlay(ref: pd.DataFrame) -> None:
    if GOOD_REFERENCE_INFO.exists():
        info = pd.read_csv(GOOD_REFERENCE_INFO).iloc[0]
        ref_path = ROOT / "analysis/am_cvn/reduced/AM_CVn" / info.reduced_file
    else:
        good_reference_sources()
        info = pd.read_csv(GOOD_REFERENCE_INFO).iloc[0]
        ref_path = ROOT / "analysis/am_cvn/reduced/AM_CVn" / info.reduced_file
    data = fits.getdata(ref_path).astype(float)
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    lo, hi = np.nanpercentile(data, [2, 99.7])
    ax.imshow(data, origin="lower", cmap="gray", vmin=lo, vmax=hi)
    for row in ref[ref.used_for_photometry].itertuples(index=False):
        color = "tab:red" if row.label == "target" else "tab:cyan"
        circ = plt.Circle((row.ref_x, row.ref_y), 18, fill=False, ec=color, lw=1.2)
        ax.add_patch(circ)
        ax.text(row.ref_x + 12, row.ref_y + 12, row.label, color="yellow", fontsize=8)
    ax.set_title("AM CVn AAVSO X42421BZ sequence overlay")
    ax.set_xlabel("x [px]")
    ax.set_ylabel("y [px]")
    fig.savefig(OUT / "AM_CVn_sequence_reference_overlay.png", dpi=170)
    plt.close(fig)


def plot_qc(phot: pd.DataFrame) -> None:
    sat = (
        phot.groupby(["object_id", "flag"])
        .size()
        .reset_index(name="n")
        .pivot_table(index="object_id", columns="flag", values="n", fill_value=0)
    )
    sat.to_csv(OUT / "photometry_flag_summary.csv")
    fig, ax = plt.subplots(figsize=(8, 4), constrained_layout=True)
    for oid, g in phot[(phot.flag == "ok") & phot.object_id.isin(["target", *COMPARISON_LABELS])].groupby("object_id"):
        y = g.flux_rate.to_numpy(float)
        ax.plot((g.jd - phot.jd.min()) * 24 * 60, y / np.nanmedian(y), ".", ms=3, label=oid)
    ax.set_xlabel("Minutes")
    ax.set_ylabel("Flux rate / median")
    ax.set_title("AM CVn raw object flux-rate stability")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(OUT / "AM_CVn_object_flux_rate_qc.png", dpi=160)
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    cfg = load_config(CONFIG)
    chart = json.load(open(CHART))
    seq = sky_offsets(chart)
    ref = reference_positions(seq)
    seq.merge(ref[["label", "ref_x", "ref_y", "used_for_photometry"]], on="label").to_csv(
        ROOT / "analysis/am_cvn/aavso/X42421BZ_sequence_pixel_table.csv", index=False
    )
    plot_reference_overlay(ref)
    phot, diff = frame_photometry(cfg, ref)
    controls = comparison_controls(phot)
    plot_qc(phot)
    summary = period_analysis(diff)
    print("AM CVn sequence labels in X42421BZ:", ", ".join(seq.label.astype(str)))
    print("Comparison labels used:", ", ".join(COMPARISON_LABELS))
    print("Photometry rows:", len(phot), "good target frames:", int(diff.good.sum()))
    print(controls.to_string(index=False))
    print(summary.to_string(index=False))
    print(f"Wrote AM CVn products to {OUT} and {PER_OUT}")


if __name__ == "__main__":
    main()
