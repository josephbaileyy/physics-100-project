# ph100-delta-scutis

Lomb-Scargle period analysis of δ Scuti candidates from ASAS-SN photometry, for PH 100.

## Layout

- `data/` — ASAS-SN light curves (HJD, mag, mag_err, …)
- `lomb_scargle.py` — runs `astropy.timeseries.LombScargle` on a CSV, writes a periodogram and phase-folded plot to `results/`
- `results/` — output PNGs
- `PracticalLombScargle/` — VanderPlas (2018) reference paper as a submodule

## Usage

```bash
pip install -r requirements.txt
python lomb_scargle.py "data/ASASSN-V_J093023.98+612034.3.csv"
```

## Result so far

`ASASSN-V J093023.98+612034.3` → P ≈ 0.0615 d (88.6 min), FAP ≈ 4×10⁻⁷². Folded light curve shows an asymmetric ~0.5 mag pulsation consistent with a δ Scuti.
