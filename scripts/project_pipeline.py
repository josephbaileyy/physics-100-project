"""Project FITS reduction, differential photometry, and Lomb-Scargle workflow.

The pipeline intentionally avoids online astrometric solving. Frames are
registered by source shifts against a reference frame, and target/comparison
apertures are supplied in reference-frame pixel coordinates.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sep
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy import signal
from scipy import ndimage


DEFAULT_CONFIG = Path("configs/targets.json")


@dataclass(frozen=True)
class Source:
    x: float
    y: float
    flux: float
    peak: float
    npix: int


def load_config(path: Path) -> dict:
    with path.open() as f:
        cfg = json.load(f)
    cfg["_config_path"] = str(path)
    cfg["_repo_root"] = str(Path.cwd())
    cfg["light_dir"] = str(Path(cfg["light_dir"]).expanduser())
    cfg["calibration_dir"] = str(Path(cfg["calibration_dir"]).expanduser())
    cfg["output_dir"] = str((Path.cwd() / cfg["output_dir"]).resolve())
    return cfg


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def native_data(path: Path) -> np.ndarray:
    data = fits.getdata(path).astype("f8")
    return np.ascontiguousarray(data)


def robust_sigma(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    med = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - med))
    return float(1.4826 * mad)


def header_time(header) -> float:
    if "JD" in header:
        return float(header["JD"])
    if "MJD-OBS" in header:
        return float(header["MJD-OBS"]) + 2400000.5
    raise ValueError("FITS header does not contain JD or MJD-OBS")


def exposure_from_path(path: Path) -> float:
    header = fits.getheader(path)
    if "EXPTIME" in header:
        return float(header["EXPTIME"])
    name = path.name
    for token in name.split("_"):
        if token.endswith("s"):
            try:
                return float(token[:-1])
            except ValueError:
                pass
    raise ValueError(f"Could not determine exposure for {path}")


def median_stack(paths: list[Path], subtract: np.ndarray | None = None) -> np.ndarray:
    if not paths:
        raise ValueError("No files supplied for median stack")
    planes = []
    for path in paths:
        data = native_data(path)
        if subtract is not None:
            data = data - subtract
        planes.append(data)
    return np.nanmedian(np.stack(planes), axis=0)


def masked_flat_stack(paths: list[Path], bias: np.ndarray, dark: np.ndarray) -> np.ndarray:
    """Median-combine flats after masking bright compact residuals per frame."""
    planes = []
    for path in paths:
        data = native_data(path) - bias - dark
        smooth = ndimage.median_filter(data, size=25)
        resid = data - smooth
        sig = robust_sigma(resid)
        if not np.isfinite(sig) or sig <= 0:
            mask = np.zeros(data.shape, dtype=bool)
        else:
            mask = resid > (5.0 * sig)
            mask = ndimage.binary_dilation(mask, iterations=3)
        planes.append(np.where(mask, np.nan, data))
    flat = np.nanmedian(np.stack(planes), axis=0)
    fill = np.nanmedian(flat)
    return np.where(np.isfinite(flat), flat, fill)


def normalize_flat(flat: np.ndarray) -> np.ndarray:
    med = np.nanmedian(flat)
    if not np.isfinite(med) or med <= 0:
        raise ValueError("Flat median is not positive")
    norm = flat / med
    return np.where(np.isfinite(norm) & (norm > 0), norm, 1.0)


def save_image(path: Path, data: np.ndarray, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 6), dpi=140)
    lo, hi = np.nanpercentile(data, [2, 99.5])
    ax.imshow(data, origin="lower", cmap="gray", vmin=lo, vmax=hi)
    ax.set_title(title)
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def calibration_groups(cal_dir: Path) -> dict[str, list[Path]]:
    return {
        "bias_early": sorted(cal_dir.glob("bias_2026_05_20-000[1-5].fit")),
        "bias_late": sorted(cal_dir.glob("bias_2026_05_20-000[6-9].fit"))
        + sorted(cal_dir.glob("bias_2026_05_20-0010.fit")),
        "dark_0p3": sorted(cal_dir.glob("darks_0.3s_*.fit")),
        "dark_4": sorted(cal_dir.glob("darks_4s_*.fit")),
        "dark_6": sorted(cal_dir.glob("darks_6s_*.fit")),
        "dark_8": sorted(cal_dir.glob("darks_8s_*.fit")),
        "dark_35": sorted(cal_dir.glob("darks_35s_*.fit")),
        "dark_40": sorted(cal_dir.glob("darks_40s_*.fit")),
        "dark_45": sorted(cal_dir.glob("darks_45s_*.fit")),
        "flat_0p3": sorted(cal_dir.glob("flats_V_0.3s_*.fit")),
        "flat_4": sorted(cal_dir.glob("flat_V_4s_*.fit")),
        "flat_6": sorted(cal_dir.glob("flat_V_6s_*.fit")),
    }


def command_inventory(args) -> None:
    cfg = load_config(args.config)
    light_dir = Path(cfg["light_dir"])
    cal_dir = Path(cfg["calibration_dir"])
    out = ensure_dir(Path(cfg["output_dir"]) / "inventory")

    rows = []
    for path in sorted(light_dir.glob("*.fit")) + sorted(cal_dir.glob("*.fit")):
        header = fits.getheader(path)
        rows.append(
            {
                "file": path.name,
                "directory": "lights" if path.parent == light_dir else "calibration",
                "object": header.get("OBJECT", ""),
                "filter": header.get("FILTER", ""),
                "exptime": header.get("EXPTIME", np.nan),
                "jd": header.get("JD", np.nan),
                "date_obs": header.get("DATE-OBS", ""),
                "xbinning": header.get("XBINNING", np.nan),
                "ybinning": header.get("YBINNING", np.nan),
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(out / "fits_inventory.csv", index=False)
    summary = (
        df.groupby(["directory", "object", "filter", "exptime"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["directory", "object", "exptime"])
    )
    summary.to_csv(out / "fits_inventory_summary.csv", index=False)
    print(summary.to_string(index=False))
    print(f"Wrote {out / 'fits_inventory.csv'}")


def command_calibrate(args) -> None:
    cfg = load_config(args.config)
    cal_dir = Path(cfg["calibration_dir"])
    out = ensure_dir(Path(cfg["output_dir"]) / "calibration")
    groups = calibration_groups(cal_dir)

    master_bias_early = median_stack(groups["bias_early"])
    master_bias_late = median_stack(groups["bias_late"])
    fits.writeto(out / "master_bias_early.fits", master_bias_early.astype("f4"), overwrite=True)
    fits.writeto(out / "master_bias_late.fits", master_bias_late.astype("f4"), overwrite=True)

    darks: dict[str, np.ndarray] = {}
    for label in ["0p3", "4", "6", "8", "35", "40", "45"]:
        dark = median_stack(groups[f"dark_{label}"], subtract=master_bias_late)
        darks[label] = dark
        fits.writeto(out / f"master_dark_{label}s.fits", dark.astype("f4"), overwrite=True)

    candidates = {
        "flat_0p3_late_bias": normalize_flat(
            median_stack(groups["flat_0p3"], subtract=master_bias_late + darks["0p3"])
        ),
        "flat_0p3_early_bias": normalize_flat(
            median_stack(groups["flat_0p3"], subtract=master_bias_early + darks["0p3"])
        ),
        "flat_4s_masked": normalize_flat(
            masked_flat_stack(groups["flat_4"], master_bias_late, darks["4"])
        ),
        "flat_6s_masked": normalize_flat(
            masked_flat_stack(groups["flat_6"], master_bias_late, darks["6"])
        ),
        "no_flat": np.ones_like(master_bias_late, dtype=float),
    }

    qc_rows = []
    for name, flat in candidates.items():
        fits.writeto(out / f"{name}.fits", flat.astype("f4"), overwrite=True)
        save_image(out / f"{name}.png", flat, name)
        high = np.mean(np.abs(flat - 1.0) > 0.10)
        qc_rows.append(
            {
                "candidate": name,
                "median": float(np.nanmedian(flat)),
                "robust_sigma": robust_sigma(flat),
                "p01": float(np.nanpercentile(flat, 1)),
                "p99": float(np.nanpercentile(flat, 99)),
                "frac_abs_deviation_gt_10pct": float(high),
            }
        )
    qc = pd.DataFrame(qc_rows).sort_values(
        ["frac_abs_deviation_gt_10pct", "robust_sigma"], ascending=True
    )
    preferred = qc[qc["candidate"] != "no_flat"].iloc[0]["candidate"]
    qc["preferred_by_metric"] = qc["candidate"] == preferred
    qc.to_csv(out / "flat_qc.csv", index=False)
    (out / "preferred_flat.txt").write_text(str(preferred) + "\n")
    print(qc.to_string(index=False))
    print(f"Preferred flat candidate by metric: {preferred}")


def load_master_dark(cal_out: Path, exptime: float) -> np.ndarray:
    label = str(int(round(exptime))) if abs(exptime - round(exptime)) < 1e-6 else str(exptime).replace(".", "p")
    path = cal_out / f"master_dark_{label}s.fits"
    if not path.exists():
        raise FileNotFoundError(f"Missing matched dark for {exptime}s: {path}")
    return native_data(path)


def selected_flat(cal_out: Path, flat_name: str | None) -> tuple[str, np.ndarray]:
    if flat_name is None:
        pref = (cal_out / "preferred_flat.txt").read_text().strip()
    else:
        pref = flat_name
    return pref, native_data(cal_out / f"{pref}.fits")


def target_light_files(cfg: dict, target: str) -> list[Path]:
    target_cfg = cfg["targets"][target]
    light_dir = Path(cfg["light_dir"])
    globs = target_cfg.get("light_globs", [target_cfg.get("light_glob")])
    files: list[Path] = []
    for pattern in globs:
        if pattern:
            files.extend(light_dir.glob(pattern))
    return sorted(set(files))


def command_reduce(args) -> None:
    cfg = load_config(args.config)
    cal_out = Path(cfg["output_dir"]) / "calibration"
    if not (cal_out / "master_bias_late.fits").exists():
        command_calibrate(args)

    targets = list(cfg["targets"]) if args.target == "all" else [args.target]
    bias = native_data(cal_out / "master_bias_late.fits")
    flat_name, flat = selected_flat(cal_out, args.flat)
    for target in targets:
        files = target_light_files(cfg, target)
        out = ensure_dir(Path(cfg["output_dir"]) / "reduced" / target)
        rows = []
        for path in files:
            data, header = fits.getdata(path, header=True)
            exptime = float(header.get("EXPTIME", exposure_from_path(path)))
            dark = load_master_dark(cal_out, exptime)
            cal = (data.astype("f8") - bias - dark) / flat
            out_path = out / path.name.replace(".fit", f"_{flat_name}_cal.fits")
            header["HISTORY"] = f"May20 pipeline calibration with {flat_name}"
            fits.writeto(out_path, cal.astype("f4"), header, overwrite=True)
            rows.append(
                {
                    "raw_file": path.name,
                    "reduced_file": out_path.name,
                    "flat": flat_name,
                    "exptime": exptime,
                    "jd": header_time(header),
                    "max_adu_raw": float(np.nanmax(data)),
                    "max_adu_cal": float(np.nanmax(cal)),
                }
            )
        manifest = pd.DataFrame(rows)
        manifest.to_csv(out / f"reduction_manifest_{flat_name}.csv", index=False)
        manifest.to_csv(out / "reduction_manifest.csv", index=False)
        print(f"{target}: wrote {len(rows)} reduced frames to {out}")


def background_subtracted(data: np.ndarray) -> tuple[np.ndarray, float]:
    background = ndimage.median_filter(data, size=41)
    sub = data - background
    sig = robust_sigma(sub)
    return sub, sig


def detect_sources(data: np.ndarray, threshold_sigma: float = 6.0) -> list[Source]:
    data = np.ascontiguousarray(data.astype("f8"))
    mask = ~np.isfinite(data)
    bkg = sep.Background(data, mask=mask, bw=128, bh=128)
    sub = data - bkg.back()
    extracted = sep.extract(
        sub,
        thresh=threshold_sigma,
        err=bkg.globalrms,
        mask=mask,
        deblend_cont=0.002,
    )
    sources: list[Source] = []
    for obj in extracted:
        if not (
            obj["a"] > 0.7
            and obj["b"] > 0.7
            and obj["a"] < 12
            and obj["b"] < 12
            and obj["flag"] == 0
        ):
            continue
        x = float(obj["x"])
        y = float(obj["y"])
        ix = int(round(x))
        iy = int(round(y))
        if iy < 0 or iy >= data.shape[0] or ix < 0 or ix >= data.shape[1]:
            peak = np.nan
        else:
            peak = float(data[iy, ix])
        sources.append(Source(x=x, y=y, flux=float(obj["flux"]), peak=peak, npix=int(obj["npix"])))
    sources.sort(key=lambda s: s.flux, reverse=True)
    return sources[:120]


def reduced_manifest(cfg: dict, target: str, flat: str | None = None) -> pd.DataFrame:
    reduced_dir = Path(cfg["output_dir"]) / "reduced" / target
    if flat:
        path = reduced_dir / f"reduction_manifest_{flat}.csv"
    else:
        path = reduced_dir / "reduction_manifest.csv"
    if not path.exists():
        raise FileNotFoundError(f"Run reduce first: {path}")
    return pd.read_csv(path)


def reference_frame(cfg: dict, target: str, flat: str | None = None) -> tuple[Path, np.ndarray, pd.Series]:
    manifest = reduced_manifest(cfg, target, flat)
    idx = manifest["max_adu_raw"].astype(float).idxmin()
    row = manifest.loc[idx]
    path = Path(cfg["output_dir"]) / "reduced" / target / row["reduced_file"]
    return path, native_data(path), row


def source_table(sources: list[Source], limit: int = 200) -> pd.DataFrame:
    rows = [
        {"rank": i + 1, "x": s.x, "y": s.y, "flux": s.flux, "peak": s.peak, "npix": s.npix}
        for i, s in enumerate(sources[:limit])
    ]
    return pd.DataFrame(rows)


def command_propose_stars(args) -> None:
    cfg = load_config(args.config)
    if args.target == "all":
        for target in cfg["targets"]:
            ns = argparse.Namespace(**vars(args))
            ns.target = target
            command_propose_stars(ns)
        return

    ref_path, ref_data, _ = reference_frame(cfg, args.target, getattr(args, "flat", None))
    sources = detect_sources(ref_data)
    out = ensure_dir(Path(cfg["output_dir"]) / "qc" / args.target)
    table = source_table(sources)
    table.to_csv(out / "proposed_sources.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 7), dpi=150)
    lo, hi = np.nanpercentile(ref_data, [2, 99.5])
    ax.imshow(ref_data, origin="lower", cmap="gray", vmin=lo, vmax=hi)
    for row in table.head(40).itertuples(index=False):
        ax.scatter(row.x, row.y, s=32, facecolors="none", edgecolors="tab:cyan", lw=0.8)
        ax.text(row.x + 6, row.y + 6, str(row.rank), color="yellow", fontsize=6)
    ax.set_title(f"{args.target} proposed sources on {ref_path.name}")
    ax.set_xlabel("x [pix]")
    ax.set_ylabel("y [pix]")
    fig.tight_layout()
    fig.savefig(out / "proposed_sources_overlay.png")
    plt.close(fig)
    print(f"Wrote {out / 'proposed_sources.csv'}")


def register_shift(ref_sources: list[Source], sources: list[Source]) -> tuple[bool, float, float, int]:
    ref = ref_sources[:80]
    cur = sources[:80]
    if len(ref) < 5 or len(cur) < 5:
        return False, np.nan, np.nan, 0
    shifts = []
    for a in cur:
        for b in ref:
            shifts.append((a.x - b.x, a.y - b.y))
    arr = np.asarray(shifts)
    rounded = np.round(arr / 1.0).astype(int)
    uniq, counts = np.unique(rounded, axis=0, return_counts=True)
    best = uniq[np.argmax(counts)].astype(float)
    dx, dy = best
    shifted = np.array([(s.x - dx, s.y - dy) for s in cur])
    matches = 0
    for r in ref:
        dist = np.hypot(shifted[:, 0] - r.x, shifted[:, 1] - r.y)
        if np.nanmin(dist) <= 4.0:
            matches += 1
    ok = matches >= 5 and abs(dx) <= 300 and abs(dy) <= 300
    return ok, float(dx), float(dy), int(matches)


def valid_manual_objects(cfg: dict, target: str) -> list[dict]:
    objects = cfg["targets"][target].get("manual_objects", [])
    valid = []
    for obj in objects:
        if obj.get("x") is None or obj.get("y") is None:
            continue
        valid.append(obj)
    return valid


def aperture_flux(
    data: np.ndarray,
    x: float,
    y: float,
    r: float,
    rin: float,
    rout: float,
    saturation_adu: float,
) -> dict:
    height, width = data.shape
    if x < rout or y < rout or x > width - rout or y > height - rout:
        return {"flux": np.nan, "flux_err": np.nan, "saturated": False, "flag": "edge"}
    yy, xx = np.indices(data.shape)
    rr2 = (xx - x) ** 2 + (yy - y) ** 2
    aper = rr2 <= r**2
    ann = (rr2 >= rin**2) & (rr2 <= rout**2)
    if ann.sum() < 20:
        return {"flux": np.nan, "flux_err": np.nan, "saturated": False, "flag": "bad_annulus"}
    bkg = np.nanmedian(data[ann])
    bkg_std = np.nanstd(data[ann])
    aper_values = data[aper]
    flux = float(np.nansum(aper_values - bkg))
    err = float(math.sqrt(max(abs(flux), 0.0) + aper.sum() * bkg_std**2))
    saturated = bool(np.nanmax(aper_values) >= saturation_adu)
    flag = "ok"
    if not np.isfinite(flux) or flux <= 0:
        flag = "nonpositive_flux"
    elif saturated:
        flag = "saturated"
    return {"flux": flux, "flux_err": err, "saturated": saturated, "flag": flag}


def centroid_near(data: np.ndarray, x: float, y: float, half_size: int = 80) -> tuple[bool, float, float]:
    height, width = data.shape
    x0 = max(0, int(round(x)) - half_size)
    x1 = min(width, int(round(x)) + half_size + 1)
    y0 = max(0, int(round(y)) - half_size)
    y1 = min(height, int(round(y)) + half_size + 1)
    if x1 - x0 < 5 or y1 - y0 < 5:
        return False, np.nan, np.nan
    cut = data[y0:y1, x0:x1].astype("f8")
    bkg = np.nanpercentile(cut, 25)
    weights = cut - bkg
    weights = np.where(weights > 0, weights, 0.0)
    total = np.nansum(weights)
    if not np.isfinite(total) or total <= 0:
        return False, np.nan, np.nan
    yy, xx = np.indices(cut.shape)
    cx = float(np.nansum(xx * weights) / total + x0)
    cy = float(np.nansum(yy * weights) / total + y0)
    if abs(cx - x) > half_size * 0.8 or abs(cy - y) > half_size * 0.8:
        return False, np.nan, np.nan
    return True, cx, cy


def local_registration_shift(data: np.ndarray, objects: list[dict]) -> tuple[bool, float, float, int]:
    anchors = [obj for obj in objects if obj["role"] in {"comparison", "check"}]
    shifts = []
    for obj in anchors:
        ok, cx, cy = centroid_near(data, float(obj["x"]), float(obj["y"]))
        if ok:
            shifts.append((cx - float(obj["x"]), cy - float(obj["y"])))
    if len(shifts) < 2:
        return False, np.nan, np.nan, len(shifts)
    shifts = np.asarray(shifts, dtype=float)
    dx = float(np.nanmedian(shifts[:, 0]))
    dy = float(np.nanmedian(shifts[:, 1]))
    scatter = np.nanmedian(np.hypot(shifts[:, 0] - dx, shifts[:, 1] - dy))
    ok = np.isfinite(dx) and np.isfinite(dy) and scatter < 8.0 and abs(dx) < 50 and abs(dy) < 50
    return bool(ok), dx, dy, len(shifts)


def highpass_for_registration(data: np.ndarray, factor: int = 4) -> np.ndarray:
    small = data[::factor, ::factor].astype("f8")
    small = small - ndimage.median_filter(small, size=21)
    lo, hi = np.nanpercentile(small, [5, 99.7])
    small = np.clip(small, lo, hi)
    small = small - np.nanmedian(small)
    norm = np.nanstd(small)
    return small / norm if np.isfinite(norm) and norm > 0 else small


def cross_correlation_shift(ref_data: np.ndarray, data: np.ndarray, factor: int = 4) -> tuple[bool, float, float, int]:
    ref = highpass_for_registration(ref_data, factor)
    cur = highpass_for_registration(data, factor)
    corr = signal.fftconvolve(cur, ref[::-1, ::-1], mode="same")
    cy, cx = np.unravel_index(np.nanargmax(corr), corr.shape)
    dy_small = cy - corr.shape[0] // 2
    dx_small = cx - corr.shape[1] // 2
    # fftconvolve(cur, reversed(ref)) returns the offset of ref relative to cur;
    # apertures defined in the reference frame need the opposite translation.
    dx = float(-dx_small * factor)
    dy = float(-dy_small * factor)
    ok = abs(dx) < 300 and abs(dy) < 300
    return ok, dx, dy, 1


def centroid_peak_in_box(
    data: np.ndarray,
    x: float,
    y: float,
    half_size: int,
    aperture: int = 12,
) -> tuple[bool, float, float]:
    height, width = data.shape
    x0 = max(0, int(round(x)) - half_size)
    x1 = min(width, int(round(x)) + half_size + 1)
    y0 = max(0, int(round(y)) - half_size)
    y1 = min(height, int(round(y)) + half_size + 1)
    if x1 - x0 < 5 or y1 - y0 < 5:
        return False, np.nan, np.nan
    cut = data[y0:y1, x0:x1].astype("f8")
    smoothed = ndimage.gaussian_filter(cut, 1.4)
    py, px = np.unravel_index(np.nanargmax(smoothed), smoothed.shape)
    peak_x = x0 + px
    peak_y = y0 + py

    ax0 = max(0, peak_x - aperture)
    ax1 = min(width, peak_x + aperture + 1)
    ay0 = max(0, peak_y - aperture)
    ay1 = min(height, peak_y + aperture + 1)
    sub = data[ay0:ay1, ax0:ax1].astype("f8")
    bkg = np.nanpercentile(sub, 20)
    weights = np.clip(sub - bkg, 0, None)
    total = np.nansum(weights)
    if not np.isfinite(total) or total <= 0:
        return True, float(peak_x), float(peak_y)
    yy, xx = np.indices(sub.shape)
    cx = float(np.nansum((xx + ax0) * weights) / total)
    cy = float(np.nansum((yy + ay0) * weights) / total)
    return True, cx, cy


def detect_v0399_dn_pair(data: np.ndarray) -> tuple[bool, dict[str, tuple[float, float]]]:
    """Find the bright V0399/DN pair directly and assign by geometry."""
    x0, x1, y0, y1 = 1760, 2045, 570, 760
    crop = data[y0:y1, x0:x1].astype("f8")
    smoothed = ndimage.gaussian_filter(crop, 1.5)
    work = smoothed.copy()
    peaks = []
    for _ in range(10):
        py, px = np.unravel_index(np.nanargmax(work), work.shape)
        x = x0 + px
        y = y0 + py
        if all(np.hypot(x - old_x, y - old_y) > 35 for old_x, old_y in peaks):
            ok, cx, cy = centroid_peak_in_box(data, x, y, half_size=16, aperture=10)
            if ok:
                peaks.append((cx, cy))
            if len(peaks) >= 2:
                break
        yy, xx = np.ogrid[: work.shape[0], : work.shape[1]]
        mask = (xx - px) ** 2 + (yy - py) ** 2 < 35**2
        work[mask] = np.nanmin(work)
    if len(peaks) < 2:
        return False, {}
    left, right = sorted(peaks[:2], key=lambda xy: xy[0])
    if right[1] < left[1]:
        target = right
        dn = left
    else:
        target = left
        dn = right
    return True, {"target": target, "DN_UMa": dn}


def direct_positions(data: np.ndarray, objects: list[dict]) -> tuple[bool, dict[str, tuple[float, float]], int]:
    ok_pair, positions = detect_v0399_dn_pair(data)
    if not ok_pair:
        return False, positions, 0
    anchors = 0
    bright_search = {
        "HD_103405": 140,
        "TYC_3452_787_1": 140,
        "HD_103187": 140,
    }
    for obj in objects:
        obj_id = obj["id"]
        if obj_id in positions:
            continue
        half_size = bright_search.get(obj_id)
        if half_size is None:
            continue
        ok, x, y = centroid_peak_in_box(data, float(obj["x"]), float(obj["y"]), half_size=half_size)
        if ok:
            positions[obj_id] = (x, y)
            anchors += 1
    return True, positions, anchors


def command_photometry(args) -> None:
    cfg = load_config(args.config)
    if args.target == "all":
        for target in cfg["targets"]:
            ns = argparse.Namespace(**vars(args))
            ns.target = target
            command_photometry(ns)
        return

    objects = valid_manual_objects(cfg, args.target)
    if not any(obj["role"] == "target" for obj in objects):
        raise ValueError(
            f"No target x/y positions set for {args.target}. Run propose-stars and edit {args.config}."
        )
    if sum(obj["role"] == "comparison" for obj in objects) < cfg.get("min_comparisons", 2):
        raise ValueError(f"Need at least {cfg.get('min_comparisons', 2)} comparison stars with x/y")

    flat = getattr(args, "flat", None)
    manifest = reduced_manifest(cfg, args.target, flat)
    ref_path, ref_data, _ = reference_frame(cfg, args.target, flat)
    use_registration = not getattr(args, "no_register", False)
    ref_sources = detect_sources(ref_data) if use_registration and getattr(args, "full_frame_register", False) else []
    out_name = args.target if not flat else f"{args.target}_{flat}"
    out = ensure_dir(Path(cfg["output_dir"]) / "photometry" / out_name)
    reduced_dir = Path(cfg["output_dir"]) / "reduced" / args.target

    r = float(cfg.get("aperture_radius", 5.0))
    rin = float(cfg.get("annulus_inner", 10.0))
    rout = float(cfg.get("annulus_outer", 18.0))
    saturation = float(cfg.get("saturation_adu", 58000))
    min_comps = int(cfg.get("min_comparisons", 2))

    rows = []
    diff_rows = []
    for frame in manifest.sort_values("jd").itertuples(index=False):
        path = reduced_dir / frame.reduced_file
        data, header = fits.getdata(path, header=True)
        data = data.astype("f8")
        direct_pos: dict[str, tuple[float, float]] = {}
        if getattr(args, "direct_v0399", False):
            ok_reg, direct_pos, matches = direct_positions(data, objects)
            dx, dy = 0.0, 0.0
        elif use_registration:
            if getattr(args, "full_frame_register", False):
                sources = detect_sources(data)
                ok_reg, dx, dy, matches = register_shift(ref_sources, sources)
            elif getattr(args, "xcorr_register", False):
                ok_reg, dx, dy, matches = cross_correlation_shift(ref_data, data)
            else:
                ok_reg, dx, dy, matches = local_registration_shift(data, objects)
        else:
            ok_reg, dx, dy, matches = True, 0.0, 0.0, 0
        exptime = float(frame.exptime)
        jd = float(frame.jd)
        object_results = []
        for obj in objects:
            if obj["id"] in direct_pos:
                x, y = direct_pos[obj["id"]]
            else:
                x = float(obj["x"]) + (dx if ok_reg else 0.0)
                y = float(obj["y"]) + (dy if ok_reg else 0.0)
            phot = aperture_flux(data, x, y, r, rin, rout, saturation)
            flux_rate = phot["flux"] / exptime if np.isfinite(phot["flux"]) else np.nan
            err_rate = phot["flux_err"] / exptime if np.isfinite(phot["flux_err"]) else np.nan
            result = {
                "file": frame.reduced_file,
                "jd": jd,
                "exptime": exptime,
                "object_id": obj["id"],
                "role": obj["role"],
                "x": x,
                "y": y,
                "registered": ok_reg,
                "dx": dx,
                "dy": dy,
                "matches": matches,
                "flux": phot["flux"],
                "flux_err": phot["flux_err"],
                "flux_rate": flux_rate,
                "flux_rate_err": err_rate,
                "saturated": phot["saturated"],
                "flag": phot["flag"] if ok_reg else "registration_failed",
            }
            object_results.append(result)
            rows.append(result)

        comps = [
            res
            for res in object_results
            if res["role"] == "comparison"
            and res["flag"] == "ok"
            and np.isfinite(res["flux_rate"])
            and res["flux_rate"] > 0
        ]
        comp_rates = np.array([res["flux_rate"] for res in comps], dtype=float)
        ensemble = float(np.nanmedian(comp_rates)) if comp_rates.size >= min_comps else np.nan
        target_res = next(res for res in object_results if res["role"] == "target")
        target_ok = target_res["flag"] == "ok" and np.isfinite(target_res["flux_rate"])
        good = bool(ok_reg and target_ok and np.isfinite(ensemble) and ensemble > 0)
        diff_flux = target_res["flux_rate"] / ensemble if good else np.nan
        diff_err = target_res["flux_rate_err"] / ensemble if good else np.nan
        reason = "ok"
        if not ok_reg:
            reason = "registration_failed"
        elif not target_ok:
            reason = target_res["flag"]
        elif comp_rates.size < min_comps:
            reason = "too_few_comparisons"
        diff_rows.append(
            {
                "file": frame.reduced_file,
                "jd": jd,
                "exptime": exptime,
                "diff_flux": diff_flux,
                "diff_flux_err": diff_err,
                "target_flux_rate": target_res["flux_rate"],
                "comparison_ensemble_flux_rate": ensemble,
                "n_comparisons": int(comp_rates.size),
                "registered": ok_reg,
                "good": good,
                "reject_reason": reason,
            }
        )

    phot = pd.DataFrame(rows)
    diff = pd.DataFrame(diff_rows)
    phot.to_csv(out / "aperture_photometry.csv", index=False)
    diff.to_csv(out / "differential_light_curve.csv", index=False)

    good = diff[diff["good"]]
    fig, ax = plt.subplots(figsize=(9, 4), dpi=150)
    ax.errorbar(good["jd"], good["diff_flux"], yerr=good["diff_flux_err"], fmt=".", ms=4)
    ax.set_xlabel("JD")
    ax.set_ylabel("differential flux")
    ax.set_title(f"{args.target} differential light curve")
    fig.tight_layout()
    fig.savefig(out / "differential_light_curve.png")
    plt.close(fig)

    rms_rows = []
    for obj_id, group in phot[phot["flag"] == "ok"].groupby("object_id"):
        values = group["flux_rate"].astype(float).to_numpy()
        rms_rows.append(
            {
                "object_id": obj_id,
                "n_ok": len(group),
                "median_flux_rate": float(np.nanmedian(values)),
                "relative_robust_rms": robust_sigma(values / np.nanmedian(values)),
            }
        )
    pd.DataFrame(rms_rows).to_csv(out / "photometry_rms_qc.csv", index=False)
    print(f"{args.target}: wrote {out / 'differential_light_curve.csv'}")


def run_periodogram(csv_path: Path, out: Path, target: str, min_period: float, max_period: float, suffix: str) -> dict:
    df = pd.read_csv(csv_path)
    good = df["good"].astype(bool) & np.isfinite(df["jd"]) & np.isfinite(df["diff_flux"])
    data = df[good].sort_values("jd")
    if len(data) < 5:
        raise ValueError(f"Need at least 5 good points for {target}; got {len(data)}")
    t_abs = data["jd"].to_numpy(dtype=float)
    t = t_abs - t_abs.min()
    y = data["diff_flux"].to_numpy(dtype=float)
    dy = data["diff_flux_err"].to_numpy(dtype=float)
    finite_dy = np.isfinite(dy) & (dy > 0)
    fill = np.nanmedian(dy[finite_dy]) if finite_dy.any() else robust_sigma(y)
    dy = np.where(finite_dy, dy, fill if np.isfinite(fill) and fill > 0 else 0.01)

    frequency, power = LombScargle(t, y - np.nanmedian(y), dy).autopower(
        minimum_frequency=1.0 / max_period,
        maximum_frequency=1.0 / min_period,
        samples_per_peak=10,
    )
    best_idx = int(np.nanargmax(power))
    best_freq = float(frequency[best_idx])
    best_period = float(1.0 / best_freq)
    best_power = float(power[best_idx])
    ls = LombScargle(t, y - np.nanmedian(y), dy)
    fap = float(ls.false_alarm_probability(best_power))

    stem = f"{target}_{suffix}"
    pd.DataFrame({"period_days": 1.0 / frequency, "frequency_cpd": frequency, "power": power}).to_csv(
        out / f"{stem}_periodogram.csv", index=False
    )

    fig, ax = plt.subplots(figsize=(9, 4), dpi=150)
    ax.errorbar(t_abs, y, yerr=dy, fmt=".", ms=4)
    ax.set_xlabel("JD")
    ax.set_ylabel("differential flux")
    ax.set_title(f"{target} Lomb-Scargle input")
    fig.tight_layout()
    fig.savefig(out / f"{stem}_timeseries.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4), dpi=150)
    ax.plot((1.0 / frequency) * 24 * 60, power, lw=0.8, color="black")
    ax.axvline(best_period * 24 * 60, color="tab:red", ls="--", label=f"{best_period * 24 * 60:.2f} min")
    ax.set_xscale("log")
    ax.set_xlabel("Period [minutes]")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title(f"{target} periodogram")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / f"{stem}_periodogram.png")
    plt.close(fig)

    phase = (t / best_period) % 1.0
    order = np.argsort(phase)
    fig, ax = plt.subplots(figsize=(9, 4), dpi=150)
    ax.errorbar(phase[order], y[order], yerr=dy[order], fmt=".", ms=4, alpha=0.75)
    ax.errorbar(phase[order] + 1, y[order], yerr=dy[order], fmt=".", ms=4, alpha=0.75)
    ax.set_xlabel("Phase")
    ax.set_ylabel("differential flux")
    ax.set_title(f"{target} folded at P = {best_period * 24 * 60:.2f} min")
    fig.tight_layout()
    fig.savefig(out / f"{stem}_folded.png")
    plt.close(fig)

    return {
        "target": target,
        "run": suffix,
        "n_good": int(len(data)),
        "min_period_days": min_period,
        "max_period_days": max_period,
        "best_period_days": best_period,
        "best_period_minutes": best_period * 24 * 60,
        "best_frequency_cpd": best_freq,
        "best_power": best_power,
        "false_alarm_probability": fap,
    }


def command_periodogram(args) -> None:
    cfg = load_config(args.config)
    if args.target == "all":
        for target in cfg["targets"]:
            ns = argparse.Namespace(**vars(args))
            ns.target = target
            command_periodogram(ns)
        return

    target_cfg = cfg["targets"][args.target]
    flat = getattr(args, "flat", None)
    out_name = args.target if not flat else f"{args.target}_{flat}"
    csv_path = Path(cfg["output_dir"]) / "photometry" / out_name / "differential_light_curve.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Run photometry first: {csv_path}")
    out = ensure_dir(Path(cfg["output_dir"]) / "periodograms" / out_name)
    min_days = float(target_cfg.get("period_min_minutes", 5.0)) / (24 * 60)
    max_days = float(target_cfg.get("period_max_minutes", 180.0)) / (24 * 60)
    rows = [
        run_periodogram(csv_path, out, args.target, min_days, max_days, "default"),
        run_periodogram(csv_path, out, args.target, min_days, max(max_days, 360.0 / (24 * 60)), "wide"),
    ]
    results = pd.DataFrame(rows)
    results.to_csv(out / "period_results.csv", index=False)
    print(results.to_string(index=False))


def command_run_all(args) -> None:
    command_inventory(args)
    command_calibrate(args)
    ns = argparse.Namespace(**vars(args))
    ns.target = "all"
    ns.flat = None
    command_reduce(ns)
    command_propose_stars(ns)
    try:
        command_photometry(ns)
        command_periodogram(ns)
    except ValueError as exc:
        print(f"Stopped before final photometry/periodogram: {exc}", file=sys.stderr)
        print("Fill manual x/y positions in configs/targets.json, then rerun photometry and periodogram.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("inventory")
    sub.add_parser("calibrate")

    reduce_p = sub.add_parser("reduce")
    reduce_p.add_argument("--target", default="all", choices=["all", "V0399_UMa", "AM_CVn"])
    reduce_p.add_argument("--flat", help="Flat candidate name; defaults to calibration/preferred_flat.txt")

    prop_p = sub.add_parser("propose-stars")
    prop_p.add_argument("--target", default="all", choices=["all", "V0399_UMa", "AM_CVn"])
    prop_p.add_argument("--flat", help="Use a flat-specific reduction manifest")

    phot_p = sub.add_parser("photometry")
    phot_p.add_argument("--target", default="all", choices=["all", "V0399_UMa", "AM_CVn"])
    phot_p.add_argument("--flat", help="Use a flat-specific reduction manifest and output suffix")
    phot_p.add_argument("--no-register", action="store_true", help="Use fixed reference-frame apertures without per-frame source registration")
    phot_p.add_argument("--full-frame-register", action="store_true", help="Use slower full-frame source matching instead of local centroid shifts")
    phot_p.add_argument("--xcorr-register", action="store_true", help="Use image cross-correlation translation registration")
    phot_p.add_argument("--direct-v0399", action="store_true", help="Directly locate the V0399/DN pair and bright V0399 comparison stars per frame")

    per_p = sub.add_parser("periodogram")
    per_p.add_argument("--target", default="all", choices=["all", "V0399_UMa", "AM_CVn"])
    per_p.add_argument("--flat", help="Read a flat-specific photometry output")

    sub.add_parser("run-all")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "inventory":
            command_inventory(args)
        elif args.command == "calibrate":
            command_calibrate(args)
        elif args.command == "reduce":
            command_reduce(args)
        elif args.command == "propose-stars":
            command_propose_stars(args)
        elif args.command == "photometry":
            command_photometry(args)
        elif args.command == "periodogram":
            command_periodogram(args)
        elif args.command == "run-all":
            command_run_all(args)
        else:
            parser.error(f"Unknown command {args.command}")
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
