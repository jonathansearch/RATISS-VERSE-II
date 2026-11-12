#!/usr/bin/env python3
"""
Verse II build pipeline.

Recomputes everything from the published data and the shared physics core,
writes results/*.json, renders figures/*.png, and seals the artefacts with
SHA-256 in results/MANIFEST.sha256.

Nothing is hand-typed into results/: every number is recomputed here.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2 import __version__                     # noqa: E402
from v2 import background as bg                # noqa: E402
from v2 import constants as C                  # noqa: E402
from v2 import data as dat                     # noqa: E402
from v2 import degeneracy as dgn               # noqa: E402
from v2 import fit as ft                       # noqa: E402
from v2 import fr as frmod                     # noqa: E402
from v2 import growth as gw                    # noqa: E402
from v2 import mg                              # noqa: E402
from v2 import plots as pl                     # noqa: E402
from v2 import signatures as sgn               # noqa: E402

RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

Z_GRID = np.linspace(0.0, 1.6, 33)
K_REF = 0.1
FR0_SCAN = [3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7]


def _dump(name: str, obj) -> Path:
    p = RESULTS / name
    p.write_text(json.dumps(obj, indent=2, sort_keys=False) + "\n")
    print(f"  wrote results/{name}")
    return p


def model_curve(E, mu, z_grid=Z_GRID, k=K_REF):
    """f sigma_8(z) with sigma_8 fixed to the Planck value."""
    rows = gw.f_sigma8(E, mu, k, redshifts=tuple(z_grid))
    return [r["f_sigma8"] for r in rows]


def build_data() -> dict:
    print("[1/5] published growth-rate data (eBOSS DR16)")
    pts = dat.eboss_dr16()
    return {
        "dataset": "eBOSS DR16 RSD-only f sigma_8 (Table 3)",
        "reference": "eBOSS Collaboration 2021, PRD 103, 083533, arXiv:2007.08991",
        "n_points": len(pts),
        "points": [{"tracer": p.tracer, "z": p.z, "f_sigma8": p.fs8,
                    "sigma": p.sigma, "ref": p.ref} for p in pts],
    }


def build_growth_curves() -> dict:
    print("[2/5] model growth curves")
    pts = dat.eboss_dr16()
    z = [p.z for p in pts]
    out = {
        "k_ref_h_mpc": K_REF,
        "lcdm": {
            "z": list(Z_GRID),
            "f_sigma8": model_curve(bg.E_lcdm, 1.0),
            "at_data_z": gw.f_sigma8(bg.E_lcdm, 1.0, K_REF, redshifts=tuple(z)),
        },
        "fr_models": {},
    }
    for fR0 in (1e-5, 1e-6):
        mu = frmod.mu_fR_callable(fR0)
        out["fr_models"][f"{fR0:.0e}"] = {
            "z": list(Z_GRID),
            "f_sigma8": model_curve(bg.E_lcdm, mu),
            "at_data_z": gw.f_sigma8(bg.E_lcdm, mu, K_REF, redshifts=tuple(z)),
        }
    return out


def build_chi2() -> dict:
    print("[3/5] chi^2 against eBOSS DR16")
    pts = dat.eboss_dr16()
    fixed = ft.scan_fr0(pts, FR0_SCAN, k=K_REF, free_amplitude=False)
    free = ft.scan_fr0(pts, FR0_SCAN, k=K_REF, free_amplitude=True)
    # the strongest bound: smallest |f_R0| excluded at 2 sigma
    d2 = [r for r in fixed["rows"] if r["delta_chi2_vs_lcdm"] < 4.0]
    bound = max((r["f_R0"] for r in d2), default=None)
    return {
        "dataset": "eBOSS DR16",
        "k_h_mpc": K_REF,
        "sigma8_note": ("sigma_8 is an input, not a prediction (Verse I, N6); "
                        "the free-amplitude fit marginalises over it"),
        "sigma8_fixed": fixed,
        "amplitude_free": free,
        "bound_2sigma_f_R0": bound,
    }


def build_degeneracy() -> dict:
    print("[4/7] the degeneracy map")
    fixed = dgn.map_theories(free_amplitude=False)
    free = dgn.map_theories(free_amplitude=True)
    sep = dgn.separability()
    return {"sigma8_fixed": fixed, "amplitude_free": free, "separability": sep}


def build_family_curves() -> dict:
    print("[3b/7] family growth curves")
    z = list(Z_GRID)
    out = {}
    for th in mg.family():
        rows = gw.f_sigma8(th.E, th.mu, K_REF, redshifts=tuple(z))
        out[th.name] = [r["f_sigma8"] for r in rows]
    return {"z": z, "k_h_mpc": K_REF, "curves": out}


def build_signatures() -> dict:
    print("[4b/8] fine signatures")
    z_list = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    curves = {th.name: sgn.gamma_curve(th, z_list) for th in mg.family()}
    disc = sgn.discriminability()
    env = sgn.environment_table()
    return {
        "z_list": z_list,
        "gamma_z_curves": curves,
        "discriminability": disc,
        "environment": env,
    }


def build_figures(data_res, growth_res, chi2_res, fam_res, dg_res,
                  sig_res) -> None:
    print("[5/8] figures")
    pts = dat.eboss_dr16()
    zz = list(Z_GRID)

    lcdm = growth_res["lcdm"]["f_sigma8"]
    fr5 = growth_res["fr_models"]["1e-05"]["f_sigma8"]
    fr6 = growth_res["fr_models"]["1e-06"]["f_sigma8"]

    models = [
        (r"$\Lambda$CDM", zz, lcdm, "k"),
        (r"$f(R)$, $|f_{R0}|=10^{-5}$", zz, fr5, pl.C_BLUE),
        (r"$f(R)$, $|f_{R0}|=10^{-6}$", zz, fr6, pl.C_GREEN),
    ]
    pl.fig01_data_fs8(pts, models, FIGURES / "fig01_data_fs8.png")
    pl.fig02_chi2_scan(chi2_res["sigma8_fixed"],
                       FIGURES / "fig02_chi2_scan.png")
    pl.fig03_fs8_residuals(pts, models, FIGURES / "fig03_fs8_residuals.png")

    # the family figure and the degeneracy map
    pal = ["k", pl.C_BLUE, "#1f8fff", pl.C_GREEN, pl.C_ORANGE,
           "#b8860b", pl.C_PURPLE, pl.C_GREY]
    curves = []
    for i, (name, yy) in enumerate(fam_res["curves"].items()):
        curves.append((name, yy, pal[i % len(pal)]))
    pl.fig05_family_fs8(zz, curves, FIGURES / "fig05_family_fs8.png")
    pl.fig04_theory_map(dg_res["sigma8_fixed"]["rows"],
                        dg_res["amplitude_free"]["rows"],
                        FIGURES / "fig04_theory_map.png")

    # step 3: fine signatures
    zl = sig_res["z_list"]
    gcurves = []
    for i, (name, yy) in enumerate(sig_res["gamma_z_curves"].items()):
        gcurves.append((name, yy, pal[i % len(pal)]))
    pl.fig06_gamma_z(zl, gcurves, FIGURES / "fig06_gamma_z.png")
    pl.fig07_discriminator(sig_res["discriminability"]["vectors"],
                           FIGURES / "fig07_discriminator.png")
    pl.fig08_environment([1e-6, 1e-5, 1e-4],
                         sig_res["environment"]["rows"],
                         FIGURES / "fig08_environment.png")

    for p in sorted(FIGURES.glob("*.png")):
        print(f"    {p.name}")


def seal() -> dict:
    print("[6/8] sealing artefacts")
    man = RESULTS / "MANIFEST.sha256"
    entries = []
    for p in sorted(list(RESULTS.glob("*.json")) + list(FIGURES.glob("*.png"))):
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        entries.append(f"{h}  {p.relative_to(ROOT).as_posix()}")
    man.write_text("\n".join(entries) + "\n")
    print(f"    {len(entries)} artefacts sealed in results/MANIFEST.sha256")
    return {"n_artefacts": len(entries), "entries": entries}


def main() -> int:
    t0 = time.time()
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    print("=" * 72)
    print(f"Verse II build  (v{__version__})")
    print("=" * 72)

    d = build_data()
    _dump("data.json", d)
    g = build_growth_curves()
    _dump("growth_curves.json", g)
    fam = build_family_curves()
    _dump("family_curves.json", fam)
    c = build_chi2()
    _dump("chi2.json", c)
    dg = build_degeneracy()
    _dump("degeneracy.json", dg)
    sig = build_signatures()
    _dump("signatures.json", sig)
    build_figures(d, g, c, fam, dg, sig)

    summary = {
        "verse_ii_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "lcdm_chi2_sigma8_fixed": c["sigma8_fixed"]["lcdm"]["chi2"],
        "lcdm_chi2_amplitude_free": c["amplitude_free"]["lcdm"]["chi2"],
        "best_fit_amplitude": c["amplitude_free"]["lcdm"]["amplitude"],
        "bound_2sigma_f_R0": c["bound_2sigma_f_R0"],
        "gamma_spread_rel": dg["separability"]["gamma_spread_rel"],
        "n_degenerate_pairs_fixed": len(dg["sigma8_fixed"]["degenerate_pairs"]),
        "n_degenerate_pairs_free": len(dg["amplitude_free"]["degenerate_pairs"]),
        "n_theories": len(dg["sigma8_fixed"]["rows"]),
        "best_signature": max(sig["discriminability"]["signature_power"], key=lambda k: sig["discriminability"]["signature_power"][k]["spread_rel"]),
    }
    _dump("build.json", summary)
    manifest = seal()
    summary["n_artefacts"] = manifest["n_artefacts"]

    print()
    print(f"Lambda-CDM chi^2 (sigma8 fixed) : {summary['lcdm_chi2_sigma8_fixed']:.2f}")
    print(f"Lambda-CDM chi^2 (amp. free)    : {summary['lcdm_chi2_amplitude_free']:.2f} "
          f"(best-fit amp = {summary['best_fit_amplitude']:.3f})")
    print(f"2-sigma bound on |f_R0|         : {summary['bound_2sigma_f_R0']}")
    print(f"theories compared               : {summary['n_theories']}")
    print(f"gamma spread                    : {100*summary['gamma_spread_rel']:.1f}%")
    print(f"degenerate pairs (fixed / free) : "
          f"{summary['n_degenerate_pairs_fixed']} / {summary['n_degenerate_pairs_free']}")
    print()
    print(f"build complete in {time.time() - t0:.2f} s")
    print(f"{manifest['n_artefacts']} artefacts sealed")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
