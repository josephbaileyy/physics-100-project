"""V0399 UMa contamination, blending, and detectability tests.

This script is intentionally analysis-oriented rather than a general pipeline
entrypoint. It uses the current best V0399 reduction products and writes a
focused suite of diagnostic CSVs/plots under analysis/project/contamination_tests.
"""
from __future__ import annotations

import json
import math
import os
from itertools import product
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle
from scipy import ndimage, optimize, stats

from project_pipeline import (
    detect_sources,
    detect_v0399_dn_pair,
    load_config,
    robust_sigma,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/targets.json"
OUT = ROOT / "analysis/project/contamination_tests/V0399_UMa"
REDUCED = ROOT / "analysis/project/reduced/V0399_UMa"
FLAT = "flat_4s_masked_clipmean"
COMMON8 = [
    "gaia_774271076054448768_src2",
    "gaia_774272862760870784_src1",
    "gaia_774287774887292672_src5",
    "gaia_774288805679445888_src7",
    "gaia_774289041901756032_src4",
    "gaia_774289488578357504_src6",
    "gaia_774290626745586944_src8",
    "gaia_774384016513181696_src3",
]
PERIODS_MIN = [4.68, 5.41, 7.65, 13.94, 14.5, 15.0, 15.3, 15.6, 16.0, 29.5]


def aperture_flux(
    data: np.ndarray,
    x: float,
    y: float,
    r: float,
    rin: float,
    rout: float,
    saturation_adu: float,
) -> dict:
    """Fast local-cutout aperture photometry for repeated grid tests."""
    height, width = data.shape
    pad = int(math.ceil(rout)) + 2
    if x < pad or y < pad or x > width - pad or y > height - pad:
        return {"flux": np.nan, "flux_err": np.nan, "saturated": False, "flag": "edge"}
    x0 = int(round(x)) - pad
    x1 = int(round(x)) + pad + 1
    y0 = int(round(y)) - pad
    y1 = int(round(y)) + pad + 1
    cut = data[y0:y1, x0:x1]
    yy, xx = np.indices(cut.shape)
    xx = xx + x0
    yy = yy + y0
    rr2 = (xx - x) ** 2 + (yy - y) ** 2
    aper = rr2 <= r**2
    ann = (rr2 >= rin**2) & (rr2 <= rout**2)
    if ann.sum() < 20:
        return {"flux": np.nan, "flux_err": np.nan, "saturated": False, "flag": "bad_annulus"}
    bkg = np.nanmedian(cut[ann])
    bkg_std = np.nanstd(cut[ann])
    aper_values = cut[aper]
    flux = float(np.nansum(aper_values - bkg))
    err = float(math.sqrt(max(abs(flux), 0.0) + aper.sum() * bkg_std**2))
    saturated = bool(np.nanmax(aper_values) >= saturation_adu)
    flag = "ok"
    if not np.isfinite(flux) or flux <= 0:
        flag = "nonpositive_flux"
    elif saturated:
        flag = "saturated"
    return {"flux": flux, "flux_err": err, "saturated": saturated, "flag": flag}


def ensure_dirs() -> None:
    for sub in ["plots", "difference_images", "tables"]:
        (OUT / sub).mkdir(parents=True, exist_ok=True)


def mad_pct(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    med = np.nanmedian(values)
    return float(100 * 1.4826 * np.nanmedian(np.abs(values / med - 1.0)))


def gaussian_aic_bic(sse: float, n: int, k: int) -> tuple[float, float]:
    """AIC/BIC for least-squares fits with an estimated Gaussian variance."""
    variance = max(float(sse) / max(int(n), 1), 1e-18)
    aic = n * math.log(variance) + 2 * k
    bic = n * math.log(variance) + k * math.log(n)
    return float(aic), float(bic)


def load_geometry() -> dict:
    cfg = load_config(CONFIG)
    objects = {o["id"]: o for o in cfg["targets"]["V0399_UMa"]["manual_objects"]}
    v8 = pd.read_csv(
        ROOT / "analysis/project/qc/V0399_UMa/sep_v8_reference_solve/v8_ref_gaia_sep_matches.csv"
    ).set_index("id")

    return {
        "V_8 generic": {
            "pattern": f"light_V_8s_2026_05_20-*_{FLAT}_cal.fits",
            "ref_target": v8.loc["target_V0399_direct"][["sep_x", "sep_y"]].to_numpy(float),
            "ref_dn": v8.loc["DN_UMa_direct"][["sep_x", "sep_y"]].to_numpy(float),
            "comps": {cid: v8.loc[cid][["sep_x", "sep_y"]].to_numpy(float) for cid in COMMON8},
            "hd": None,
        },
        "V0399 named": {
            "pattern": f"light_V0399_V_8s_2026_05_20-*_{FLAT}_cal.fits",
            "ref_target": np.array([objects["target"]["x"], objects["target"]["y"]], dtype=float),
            "ref_dn": np.array([objects["DN_UMa"]["x"], objects["DN_UMa"]["y"]], dtype=float),
            "comps": {cid: np.array([objects[cid]["x"], objects[cid]["y"]], dtype=float) for cid in COMMON8},
            "hd": np.array([objects["HD_103187"]["x"], objects["HD_103187"]["y"]], dtype=float),
        },
    }


def similarity_predict(xy: np.ndarray, ref_target: np.ndarray, ref_dn: np.ndarray, cur_target, cur_dn):
    ref_mid = (ref_target + ref_dn) / 2.0
    ref_vec = ref_target - ref_dn
    cur_target = np.asarray(cur_target, dtype=float)
    cur_dn = np.asarray(cur_dn, dtype=float)
    cur_mid = (cur_target + cur_dn) / 2.0
    cur_vec = cur_target - cur_dn
    theta = math.atan2(cur_vec[1], cur_vec[0]) - math.atan2(ref_vec[1], ref_vec[0])
    scale = np.hypot(*cur_vec) / np.hypot(*ref_vec)
    c, s = math.cos(theta), math.sin(theta)
    rot = np.array([[c, -s], [s, c]])
    return cur_mid + scale * (rot @ (np.asarray(xy, dtype=float) - ref_mid))


def frame_list(pattern: str) -> list[tuple[float, float, Path]]:
    frames = []
    for path in sorted(REDUCED.glob(pattern)):
        h = fits.getheader(path)
        frames.append((float(h["JD"]), float(h.get("EXPTIME", 8.0)), path))
    return sorted(frames)


def detect_positions(data: np.ndarray, geom: dict) -> dict[str, tuple[float, float]]:
    ok, pair = detect_v0399_dn_pair(data)
    if not ok:
        return {}
    sources = detect_sources(data, threshold_sigma=5.0)
    src_xy = np.array([[s.x, s.y] for s in sources], dtype=float) if sources else np.empty((0, 2))
    pos = {"V0399": pair["target"], "DN_UMa": pair["DN_UMa"]}
    for cid, refxy in geom["comps"].items():
        pred = similarity_predict(refxy, geom["ref_target"], geom["ref_dn"], pair["target"], pair["DN_UMa"])
        if src_xy.size:
            dist = np.hypot(src_xy[:, 0] - pred[0], src_xy[:, 1] - pred[1])
            idx = int(np.nanargmin(dist))
            if dist[idx] <= 70:
                pos[cid] = (sources[idx].x, sources[idx].y)
    if geom.get("hd") is not None:
        pred = similarity_predict(geom["hd"], geom["ref_target"], geom["ref_dn"], pair["target"], pair["DN_UMa"])
        if src_xy.size:
            dist = np.hypot(src_xy[:, 0] - pred[0], src_xy[:, 1] - pred[1])
            idx = int(np.nanargmin(dist))
            if dist[idx] <= 100:
                pos["HD_103187"] = (sources[idx].x, sources[idx].y)
    return pos


def differential_series(rows: pd.DataFrame, target: str) -> pd.DataFrame:
    out = []
    for (seq, file), g in rows.groupby(["sequence", "file"]):
        comps = g[(g.role == "comparison") & (g.flag == "ok") & (g.flux_rate > 0)]
        targ = g[(g.object_id == target) & (g.flag == "ok")]
        if len(targ) != 1 or len(comps) < 2:
            continue
        ens = float(np.nanmedian(comps.flux_rate))
        out.append(
            {
                "sequence": seq,
                "file": file,
                "jd": float(g.jd.iloc[0]),
                "target": target,
                "diff_flux": float(targ.flux_rate.iloc[0] / ens),
                "n_comps": int(len(comps)),
            }
        )
    df = pd.DataFrame(out)
    if df.empty:
        return df
    df = df.sort_values("jd")
    df["minutes_combined"] = (df.jd - df.jd.min()) * 24 * 60
    df["norm_flux"] = np.nan
    for seq, idx in df.groupby("sequence").groups.items():
        df.loc[idx, "norm_flux"] = df.loc[idx, "diff_flux"] / np.nanmedian(df.loc[idx, "diff_flux"])
    return df


def fit_sine(t_days: np.ndarray, y: np.ndarray, period_min: float, seq: np.ndarray | None = None) -> dict:
    if seq is None:
        seq = np.array(["all"] * len(y))
    seqs = sorted(set(seq))
    base = [(seq == s).astype(float) for s in seqs]
    phase = 2 * np.pi * (t_days * 24 * 60) / period_min
    x_const = np.column_stack(base)
    x_sine = np.column_stack(base + [np.sin(phase), np.cos(phase)])
    b0 = np.linalg.lstsq(x_const, y, rcond=None)[0]
    b1 = np.linalg.lstsq(x_sine, y, rcond=None)[0]
    r0 = y - x_const @ b0
    r1 = y - x_sine @ b1
    chi0 = float(np.sum(r0**2))
    chi1 = float(np.sum(r1**2))
    n = len(y)
    k0 = x_const.shape[1]
    k1 = x_sine.shape[1]
    _, bic0 = gaussian_aic_bic(chi0, n, k0)
    _, bic1 = gaussian_aic_bic(chi1, n, k1)
    dchi = chi0 - chi1
    return {
        "period_min": period_min,
        "amp_pct": float(100 * np.hypot(b1[-2], b1[-1])),
        "delta_chi2_unweighted": dchi,
        "nested_p_unweighted": float(stats.chi2.sf(max(dchi, 0), k1 - k0)),
        "bic_constant": bic0,
        "bic_sine": bic1,
        "rms_constant_pct": float(100 * np.sqrt(np.mean(r0**2))),
        "rms_sine_pct": float(100 * np.sqrt(np.mean(r1**2))),
    }


def lomb_scargle_summary(df: pd.DataFrame) -> dict:
    if len(df) < 10:
        return {"ls_period_min": np.nan, "ls_power": np.nan, "ls_fap": np.nan}
    t = df.jd.to_numpy(float) - df.jd.min()
    y = df.norm_flux.to_numpy(float)
    seq = df.sequence.to_numpy()
    y0 = y.copy()
    for s in sorted(set(seq)):
        m = seq == s
        y0[m] -= np.nanmean(y0[m])
    ls = LombScargle(t, y0, center_data=True, fit_mean=True)
    freq, power = ls.autopower(
        minimum_frequency=1 / (45 / (24 * 60)),
        maximum_frequency=1 / (3 / (24 * 60)),
        samples_per_peak=25,
    )
    idx = int(np.nanargmax(power))
    return {
        "ls_period_min": float((1 / freq[idx]) * 24 * 60),
        "ls_power": float(power[idx]),
        "ls_fap": float(ls.false_alarm_probability(power[idx])),
    }


def aperture_annulus_grid(geometry: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    radii = [3, 4, 5, 6, 7, 8, 10]
    annuli = [(8, 14), (10, 18), (14, 24), (18, 30)]
    rows = []
    summaries = []
    for sequence, geom in geometry.items():
        frames = frame_list(geom["pattern"])
        print(f"[aperture-grid] {sequence}: {len(frames)} frames", flush=True)
        for frame_idx, (jd, exptime, path) in enumerate(frames, start=1):
            if frame_idx % 25 == 0:
                print(f"[aperture-grid] {sequence}: {frame_idx}/{len(frames)}", flush=True)
            data = fits.getdata(path).astype(float)
            positions = detect_positions(data, geom)
            if not positions:
                continue
            for radius, (rin, rout) in product(radii, annuli):
                for oid, xy in positions.items():
                    role = "comparison" if oid in COMMON8 else ("target" if oid == "V0399" else "check")
                    phot = aperture_flux(data, xy[0], xy[1], radius, rin, rout, 58000.0)
                    rows.append(
                        {
                            "sequence": sequence,
                            "file": path.name,
                            "jd": jd,
                            "exptime": exptime,
                            "radius": radius,
                            "annulus_inner": rin,
                            "annulus_outer": rout,
                            "object_id": oid,
                            "role": role,
                            "x": xy[0],
                            "y": xy[1],
                            "flux_rate": phot["flux"] / exptime if np.isfinite(phot["flux"]) else np.nan,
                            "flag": phot["flag"],
                        }
                    )
    phot = pd.DataFrame(rows)
    phot.to_csv(OUT / "tables/aperture_annulus_grid_photometry.csv", index=False)
    for (radius, rin, rout), sub in phot.groupby(["radius", "annulus_inner", "annulus_outer"]):
        v = differential_series(sub, "V0399")
        dn = differential_series(sub, "DN_UMa")
        hd = differential_series(sub, "HD_103187")
        joined = v[["sequence", "file", "norm_flux", "jd"]].merge(
            dn[["sequence", "file", "norm_flux"]], on=["sequence", "file"], suffixes=("_v", "_dn")
        )
        corr = np.nan
        corr_p = np.nan
        if len(joined) > 5:
            corr, corr_p = stats.pearsonr(joined.norm_flux_v - 1, joined.norm_flux_dn - 1)
        sine15 = fit_sine(
            v.jd.to_numpy(float) - v.jd.min(),
            v.norm_flux.to_numpy(float),
            15.3,
            v.sequence.to_numpy(),
        ) if len(v) else {}
        summaries.append(
            {
                "radius": radius,
                "annulus_inner": rin,
                "annulus_outer": rout,
                "V0399_N": len(v),
                "DN_N": len(dn),
                "HD103187_N": len(hd),
                "V0399_scatter_pct": mad_pct(v.norm_flux) if len(v) else np.nan,
                "DN_scatter_pct": mad_pct(dn.norm_flux) if len(dn) else np.nan,
                "HD103187_scatter_pct": mad_pct(hd.norm_flux) if len(hd) else np.nan,
                "V0399_DN_corr_r": corr,
                "V0399_DN_corr_p": corr_p,
                "V0399_15p3_amp_pct": sine15.get("amp_pct", np.nan),
                "V0399_15p3_delta_chi2": sine15.get("delta_chi2_unweighted", np.nan),
                **lomb_scargle_summary(v),
            }
        )
    summary = pd.DataFrame(summaries).sort_values(["radius", "annulus_inner"])
    summary.to_csv(OUT / "aperture_annulus_grid_summary.csv", index=False)
    plot_aperture_summary(summary)
    return phot, summary


def plot_aperture_summary(summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for ann, g in summary.groupby(summary.annulus_inner.astype(str) + "-" + summary.annulus_outer.astype(str)):
        best = g.sort_values("V0399_scatter_pct")
        ax.plot(best.radius, best.V0399_scatter_pct, "o-", label=f"ann {ann}")
    ax.set_xlabel("Aperture radius [px]")
    ax.set_ylabel("V0399 robust scatter [%]")
    ax.set_title("Aperture/annulus sensitivity")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(OUT / "plots/aperture_radius_dependence.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for ann, g in summary.groupby(summary.annulus_inner.astype(str) + "-" + summary.annulus_outer.astype(str)):
        ax.plot(g.radius, g.V0399_DN_corr_r, "o-", label=f"ann {ann}")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlabel("Aperture radius [px]")
    ax.set_ylabel("corr(V0399, DN)")
    ax.set_title("V0399-DN correlation vs aperture")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(OUT / "plots/v0399_dn_correlation_vs_aperture.png", dpi=160)
    plt.close(fig)


def difference_image_photometry(geometry: dict) -> pd.DataFrame:
    rows = []
    examples = []
    for sequence, geom in geometry.items():
        frames = frame_list(geom["pattern"])
        print(f"[difference] {sequence}: {len(frames)} frames", flush=True)
        detected = []
        for frame_idx, (jd, exptime, path) in enumerate(frames, start=1):
            if frame_idx % 25 == 0:
                print(f"[difference] {sequence}: detect {frame_idx}/{len(frames)}", flush=True)
            data = fits.getdata(path).astype(float)
            pos = detect_positions(data, geom)
            if {"V0399", "DN_UMa"}.issubset(pos):
                detected.append((jd, exptime, path, data, pos))
        if len(detected) < 10:
            continue
        # Use the median target rate frames as reference candidates.
        rates = []
        for jd, exptime, path, data, pos in detected:
            phot = aperture_flux(data, pos["V0399"][0], pos["V0399"][1], 5, 10, 18, 58000)
            rates.append(phot["flux"] / exptime if np.isfinite(phot["flux"]) else np.nan)
        order = np.argsort(np.abs(np.asarray(rates) - np.nanmedian(rates)))[:15]
        ref_jd, ref_exp, ref_path, ref_data, ref_pos = detected[int(order[0])]
        ref_stack = []
        for idx in order:
            jd, exptime, path, data, pos = detected[int(idx)]
            dy = ref_pos["V0399"][1] - pos["V0399"][1]
            dx = ref_pos["V0399"][0] - pos["V0399"][0]
            shifted = ndimage.shift(data, shift=(dy, dx), order=1, mode="nearest")
            ref_stack.append(shifted / np.nanmedian(shifted))
        ref_norm = np.nanmedian(np.stack(ref_stack), axis=0)
        ref_image = ref_norm * np.nanmedian(ref_data)
        for frame_idx, (jd, exptime, path, data, pos) in enumerate(detected, start=1):
            if frame_idx % 25 == 0:
                print(f"[difference] {sequence}: subtract {frame_idx}/{len(detected)}", flush=True)
            dy = ref_pos["V0399"][1] - pos["V0399"][1]
            dx = ref_pos["V0399"][0] - pos["V0399"][0]
            shifted = ndimage.shift(data, shift=(dy, dx), order=1, mode="nearest")
            scale = np.nanmedian(shifted) / np.nanmedian(ref_image)
            diff = shifted - ref_image * scale
            for radius in [4, 5, 6]:
                for oid in ["V0399", "DN_UMa"]:
                    xy = ref_pos["V0399"] if oid == "V0399" else ref_pos["DN_UMa"]
                    phot = aperture_flux(diff, xy[0], xy[1], radius, radius + 4, radius + 12, 1e12)
                    rows.append(
                        {
                            "sequence": sequence,
                            "file": path.name,
                            "jd": jd,
                            "radius": radius,
                            "object_id": oid,
                            "residual_flux": phot["flux"],
                            "flag": phot["flag"],
                        }
                    )
            if len(examples) < 4 and path != ref_path:
                examples.append((sequence, ref_image, shifted, diff, ref_pos, path.name))
    df = pd.DataFrame(rows).sort_values("jd")
    if not df.empty:
        df["minutes_combined"] = (df.jd - df.jd.min()) * 24 * 60
    df.to_csv(OUT / "difference_image_residual_light_curves.csv", index=False)
    plot_difference_examples(examples)
    plot_difference_summary(df)
    return df


def plot_difference_examples(examples) -> None:
    if not examples:
        return
    fig, axes = plt.subplots(len(examples), 3, figsize=(12, 3.5 * len(examples)), constrained_layout=True)
    if len(examples) == 1:
        axes = np.array([axes])
    for row, (sequence, ref_image, shifted, diff, ref_pos, name) in enumerate(examples):
        for col, (img, title) in enumerate([(ref_image, "reference"), (shifted, "registered"), (diff, "difference")]):
            ax = axes[row, col]
            lo, hi = np.nanpercentile(img, [1, 99.7]) if col < 2 else np.nanpercentile(img, [0.5, 99.5])
            ax.imshow(img, origin="lower", cmap="gray", vmin=lo, vmax=hi)
            for oid, color in [("V0399", "red"), ("DN_UMa", "cyan")]:
                xy = ref_pos["V0399"] if oid == "V0399" else ref_pos["DN_UMa"]
                circ = plt.Circle((xy[0], xy[1]), 18, fill=False, ec=color, lw=1)
                ax.add_patch(circ)
            ax.set_xlim(ref_pos["DN_UMa"][0] - 130, ref_pos["V0399"][0] + 130)
            ax.set_ylim(min(ref_pos["V0399"][1], ref_pos["DN_UMa"][1]) - 130, max(ref_pos["V0399"][1], ref_pos["DN_UMa"][1]) + 130)
            ax.set_title(f"{sequence} {title}\n{name}" if col == 0 else title)
            ax.set_axis_off()
    fig.savefig(OUT / "plots/difference_image_contact_sheet.png", dpi=160)
    plt.close(fig)


def plot_difference_summary(df: pd.DataFrame) -> None:
    if df.empty:
        return
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for (oid, radius), g in df.groupby(["object_id", "radius"]):
        if radius == 5:
            y = g.residual_flux - np.nanmedian(g.residual_flux)
            ax.plot(g.minutes_combined, y, ".", label=oid)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlabel("Minutes")
    ax.set_ylabel("Difference residual flux, median-subtracted")
    ax.set_title("Difference-image residual light curves (r=5)")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(OUT / "plots/difference_image_residual_light_curves.png", dpi=160)
    plt.close(fig)


def load_default_series() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    grid = pd.read_csv(OUT / "tables/aperture_annulus_grid_photometry.csv")
    default = grid[(grid.radius == 5) & (grid.annulus_inner == 10) & (grid.annulus_outer == 18)]
    return (
        differential_series(default, "V0399"),
        differential_series(default, "DN_UMa"),
        differential_series(default, "HD_103187"),
    )


def two_star_model(v: pd.DataFrame, dn: pd.DataFrame) -> pd.DataFrame:
    joined = v[["sequence", "file", "jd", "norm_flux"]].merge(
        dn[["sequence", "file", "norm_flux"]], on=["sequence", "file"], suffixes=("_v", "_dn")
    )
    joined["v_resid"] = joined.norm_flux_v - 1
    joined["dn_resid"] = joined.norm_flux_dn - 1
    y = joined.v_resid.to_numpy(float)
    xdn = joined.dn_resid.to_numpy(float)
    seq = joined.sequence.to_numpy()
    seqs = sorted(set(seq))
    base = [(seq == s).astype(float) for s in seqs]
    models = []

    def fit(name, cols):
        x = np.column_stack(base + cols)
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
        pred = x @ beta
        resid = y - pred
        n, k = len(y), len(beta)
        sse = float(np.sum(resid**2))
        aic, bic = gaussian_aic_bic(sse, n, k)
        # Leave-one-out residual from hat matrix.
        h = np.sum((x @ np.linalg.pinv(x.T @ x)) * x, axis=1)
        loo = resid / np.clip(1 - h, 1e-6, None)
        models.append(
            {
                "model": name,
                "N": n,
                "k": k,
                "rms_pct": 100 * float(np.sqrt(np.mean(resid**2))),
                "loocv_rms_pct": 100 * float(np.sqrt(np.mean(loo**2))),
                "aic": aic,
                "bic": bic,
                "dn_leak_beta": beta[-1] if cols else np.nan,
            }
        )

    fit("constant_v0399", [])
    fit("v0399_mixed_with_dn", [xdn])
    for period in [15.3, 5.41]:
        phase = 2 * np.pi * ((joined.jd.to_numpy(float) - joined.jd.min()) * 24 * 60) / period
        fit(f"v0399_sine_{period:.2f}min", [np.sin(phase), np.cos(phase)])
        fit(f"v0399_sine_{period:.2f}min_plus_dn", [np.sin(phase), np.cos(phase), xdn])
    out = pd.DataFrame(models).sort_values("bic")
    out.to_csv(OUT / "two_star_contamination_model_summary.csv", index=False)
    return out


def model_basis(t_min: np.ndarray, period: float, kind: str) -> np.ndarray:
    phase = (t_min / period) % 1.0
    if kind == "sine":
        return np.column_stack([np.sin(2 * np.pi * phase), np.cos(2 * np.pi * phase)])
    if kind == "two_harmonic":
        return np.column_stack([
            np.sin(2 * np.pi * phase),
            np.cos(2 * np.pi * phase),
            np.sin(4 * np.pi * phase),
            np.cos(4 * np.pi * phase),
        ])
    if kind == "sawtooth":
        saw = 2 * phase - 1
        return saw[:, None]
    if kind == "gaussian_dip":
        widths = [0.06, 0.10, 0.16]
        return np.column_stack([-np.exp(-0.5 * np.minimum(np.abs(phase), 1 - np.abs(phase)) ** 2 / w**2) for w in widths])
    raise ValueError(kind)


def fit_template(series: pd.DataFrame, period: float, kind: str) -> dict:
    t_min = (series.jd.to_numpy(float) - series.jd.min()) * 24 * 60
    y = series.norm_flux.to_numpy(float)
    seq = series.sequence.to_numpy()
    seqs = sorted(set(seq))
    base = [(seq == s).astype(float) for s in seqs]
    x0 = np.column_stack(base)
    xt = np.column_stack(base + [c for c in model_basis(t_min, period, kind).T])
    b0 = np.linalg.lstsq(x0, y, rcond=None)[0]
    bt = np.linalg.lstsq(xt, y, rcond=None)[0]
    r0 = y - x0 @ b0
    rt = y - xt @ bt
    n, k0, kt = len(y), x0.shape[1], xt.shape[1]
    sse0 = float(np.sum(r0**2))
    sset = float(np.sum(rt**2))
    _, bic0 = gaussian_aic_bic(sse0, n, k0)
    _, bict = gaussian_aic_bic(sset, n, kt)
    return {
        "period_min": period,
        "model": kind,
        "N": n,
        "k": kt,
        "delta_chi2_unweighted": sse0 - sset,
        "nested_p_unweighted": float(stats.chi2.sf(max(sse0 - sset, 0), kt - k0)),
        "bic_constant": bic0,
        "bic_model": bict,
        "rms_model_pct": 100 * float(np.sqrt(np.mean(rt**2))),
        "amplitude_proxy_pct": 100 * float(np.nanmax(xt @ bt) - np.nanmin(xt @ bt)),
    }


def fixed_period_tests(v: pd.DataFrame, dn: pd.DataFrame, hd: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target_name, series in [("V0399", v), ("DN_UMa", dn), ("HD_103187", hd)]:
        if len(series) < 10:
            continue
        for period, kind in product(PERIODS_MIN, ["sine", "two_harmonic", "gaussian_dip", "sawtooth"]):
            rows.append({"target": target_name, **fit_template(series, period, kind)})
    out = pd.DataFrame(rows).sort_values(["target", "bic_model"])
    out.to_csv(OUT / "fixed_period_model_summary.csv", index=False)
    plot_fixed_period(v, out)
    return out


def plot_fixed_period(v: pd.DataFrame, fixed: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5), constrained_layout=True)
    phase = (((v.jd - v.jd.min()) * 24 * 60) / 15.3) % 1
    ax.plot(phase, v.norm_flux, ".", label="V0399")
    ax.plot(phase + 1, v.norm_flux, ".", color="C0")
    bins = np.linspace(0, 1, 11)
    mids = 0.5 * (bins[:-1] + bins[1:])
    means = [np.nanmean(v.norm_flux[(phase >= bins[i]) & (phase < bins[i + 1])]) for i in range(len(mids))]
    ax.plot(mids, means, "ko-", label="phase-bin mean")
    ax.plot(mids + 1, means, "ko-")
    ax.set_xlabel("Phase at 15.3 min")
    ax.set_ylabel("Normalized flux")
    ax.set_title("Fixed VSX-period fold")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(OUT / "plots/fixed_15p3min_folded.png", dpi=160)
    plt.close(fig)


def non_sinusoidal_tests(v: pd.DataFrame, dn: pd.DataFrame, hd: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for target_name, series in [("V0399", v), ("DN_UMa", dn), ("HD_103187", hd)]:
        if len(series) < 10:
            continue
        for period in [5.41, 7.65, 15.3, 29.5]:
            for kind in ["sawtooth", "gaussian_dip", "two_harmonic"]:
                rows.append({"target": target_name, **fit_template(series, period, kind)})
        # Time-domain flexible models.
        t = (series.jd.to_numpy(float) - series.jd.min()) * 24 * 60
        y = series.norm_flux.to_numpy(float)
        seq = series.sequence.to_numpy()
        seqs = sorted(set(seq))
        base = [(seq == s).astype(float) for s in seqs]
        x0 = np.column_stack(base)
        b0 = np.linalg.lstsq(x0, y, rcond=None)[0]
        r0 = y - x0 @ b0
        sse0 = float(np.sum(r0**2))
        for k in [4, 6, 8]:
            knots = np.linspace(t.min(), t.max(), k)
            basis = np.column_stack([np.maximum(0, t - knot) for knot in knots[1:-1]])
            x = np.column_stack(base + [t, *basis.T])
            beta = np.linalg.lstsq(x, y, rcond=None)[0]
            resid = y - x @ beta
            sse = float(np.sum(resid**2))
            _, bic0 = gaussian_aic_bic(sse0, len(y), x0.shape[1])
            _, bic = gaussian_aic_bic(sse, len(y), x.shape[1])
            rows.append(
                {
                    "target": target_name,
                    "period_min": np.nan,
                    "model": f"piecewise_linear_spline_k{k}",
                    "N": len(y),
                    "k": x.shape[1],
                    "delta_chi2_unweighted": sse0 - sse,
                    "nested_p_unweighted": float(stats.chi2.sf(max(sse0 - sse, 0), x.shape[1] - x0.shape[1])),
                    "bic_constant": bic0,
                    "bic_model": bic,
                    "rms_model_pct": 100 * float(np.sqrt(np.mean(resid**2))),
                    "amplitude_proxy_pct": 100 * float(np.nanmax(x @ beta) - np.nanmin(x @ beta)),
                }
            )
        # DRW/flickering proxy: AR(1) residual likelihood, compared descriptively.
        resid = r0[np.argsort(t)]
        rho = np.corrcoef(resid[:-1], resid[1:])[0, 1] if len(resid) > 3 else np.nan
        rows.append(
            {
                "target": target_name,
                "period_min": np.nan,
                "model": "damped_random_walk_proxy_ar1",
                "N": len(y),
                "k": 3,
                "delta_chi2_unweighted": np.nan,
                "nested_p_unweighted": np.nan,
                "bic_constant": gaussian_aic_bic(sse0, len(y), x0.shape[1])[1],
                "bic_model": np.nan,
                "rms_model_pct": 100 * float(np.sqrt(np.mean(r0**2))),
                "amplitude_proxy_pct": np.nan,
                "ar1_rho": rho,
            }
        )
    out = pd.DataFrame(rows).sort_values(["target", "bic_model"], na_position="last")
    out.to_csv(OUT / "non_sinusoidal_model_summary.csv", index=False)
    plot_model_bic(out)
    return out


def plot_model_bic(non: pd.DataFrame) -> None:
    v = non[non.target == "V0399"].dropna(subset=["bic_model"]).sort_values("bic_model").head(12)
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    ax.barh(np.arange(len(v)), v.bic_model - v.bic_constant)
    ax.set_yticks(np.arange(len(v)))
    ax.set_yticklabels([f"{r.model} {r.period_min:.2f}" if np.isfinite(r.period_min) else r.model for r in v.itertuples()])
    ax.axvline(0, color="0.5", lw=0.8)
    ax.set_xlabel("BIC(model) - BIC(constant)")
    ax.set_title("V0399 non-sinusoidal model comparison")
    fig.savefig(OUT / "plots/model_comparison_bic.png", dpi=160)
    plt.close(fig)


def injection_recovery(v: pd.DataFrame, hd: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(20260523)
    rows = []
    controls = [("V0399_noise", v), ("HD_103187", hd)]
    # Add five artificial stars by resampling V0399 residuals around constant flux.
    for i in range(5):
        sim = v.copy()
        residual = sim.norm_flux.to_numpy(float) - 1.0
        sim["norm_flux"] = 1.0 + rng.choice(residual, size=len(sim), replace=True)
        controls.append((f"artificial_{i+1}", sim))
    for name, base in controls:
        if len(base) < 10:
            continue
        t_min = (base.jd.to_numpy(float) - base.jd.min()) * 24 * 60
        phase = (t_min / 15.3) % 1
        for waveform in ["sine", "gaussian_dip", "sawtooth"]:
            amps = [0.002, 0.005, 0.01, 0.02] if waveform == "sine" else [0.005, 0.01, 0.02]
            for amp in amps:
                y = base.norm_flux.to_numpy(float).copy()
                if waveform == "sine":
                    y += amp * np.sin(2 * np.pi * phase)
                elif waveform == "gaussian_dip":
                    dist = np.minimum(phase, 1 - phase)
                    y -= amp * np.exp(-0.5 * dist**2 / 0.08**2)
                else:
                    y += amp * (2 * phase - 1)
                inj = base.copy()
                inj["norm_flux"] = y
                sine = fit_sine(inj.jd.to_numpy(float) - inj.jd.min(), y, 15.3, inj.sequence.to_numpy())
                ls = lomb_scargle_summary(inj)
                rows.append(
                    {
                        "control": name,
                        "waveform": waveform,
                        "injected_amp_pct": amp * 100,
                        "fit_15p3_amp_pct": sine["amp_pct"],
                        "fit_15p3_delta_chi2": sine["delta_chi2_unweighted"],
                        "fit_15p3_bic_delta": sine["bic_sine"] - sine["bic_constant"],
                        "ls_best_period_min": ls["ls_period_min"],
                        "ls_fap": ls["ls_fap"],
                        "recovered_fixed_period": bool(sine["bic_sine"] < sine["bic_constant"]),
                        "recovered_ls_near_15p3": bool(abs(ls["ls_period_min"] - 15.3) < 1.0 and ls["ls_fap"] < 0.1),
                    }
                )
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "injection_recovery_summary.csv", index=False)
    plot_injection(out)
    return out


def plot_injection(df: pd.DataFrame) -> None:
    pivot = df[df.waveform == "sine"].pivot_table(
        index="control", columns="injected_amp_pct", values="fit_15p3_bic_delta", aggfunc="median"
    )
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    im = ax.imshow(pivot.to_numpy(float), aspect="auto", cmap="coolwarm_r")
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels([f"{c:.1f}%" for c in pivot.columns])
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    ax.set_xlabel("Injected sine amplitude")
    ax.set_title("Injection recovery: BIC(sine 15.3) - BIC(constant)")
    fig.colorbar(im, ax=ax, label="negative favors injected-period model")
    fig.savefig(OUT / "plots/injection_recovery_heatmap.png", dpi=160)
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    geometry = load_geometry()
    print(f"Writing contamination tests to {OUT}", flush=True)
    aperture_annulus_grid(geometry)
    diff = difference_image_photometry(geometry)
    # Difference-image correlation table.
    corr_rows = []
    if not diff.empty:
        for radius, g in diff.groupby("radius"):
            piv = g.pivot_table(index=["sequence", "file"], columns="object_id", values="residual_flux")
            if {"V0399", "DN_UMa"}.issubset(piv.columns) and len(piv) > 5:
                r, p = stats.pearsonr(piv["V0399"], piv["DN_UMa"])
                corr_rows.append({"radius": radius, "N": len(piv), "residual_corr_r": r, "residual_corr_p": p})
    pd.DataFrame(corr_rows).to_csv(OUT / "difference_image_residual_correlation.csv", index=False)

    print("[models] loading default r=5, annulus=10-18 series", flush=True)
    v, dn, hd = load_default_series()
    print("[models] two-star contamination model", flush=True)
    two_star_model(v, dn)
    print("[models] fixed-period templates", flush=True)
    fixed_period_tests(v, dn, hd)
    print("[models] non-sinusoidal templates", flush=True)
    non_sinusoidal_tests(v, dn, hd)
    print("[models] injection recovery", flush=True)
    injection_recovery(v, hd)

    print(f"Wrote contamination tests to {OUT}")


if __name__ == "__main__":
    main()
