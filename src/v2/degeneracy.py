"""
The degeneracy map: can the data tell the theories apart?

Verse I showed that H0*t0 cannot separate the universes. Verse II has real
f sigma_8 data, so the question becomes quantitative: compute each theory's
chi^2 against eBOSS DR16 and its growth index gamma, and look for pairs that
the data cannot distinguish.

Two theories are *degenerate* when their chi^2 differ by less than the
resolution of the data. The map below makes that visible.
"""

from __future__ import annotations

import numpy as np

from . import background as bg
from . import constants as C
from . import fit as ft
from . import growth as gw
from . import mg


def growth_index_of(th, k: float = 0.1) -> tuple[float, float]:
    """(f(z=0), gamma) for a Theory."""
    return gw.growth_index(th.E, th.mu, k)


def map_theories(theories=None, points=None, k: float = 0.1,
                 free_amplitude: bool = False) -> dict:
    """chi^2, gamma and mu(a=1) for every theory in the family."""
    if theories is None:
        theories = mg.family()
    if points is None:
        from . import data as dat
        points = dat.eboss_dr16()

    rows = []
    for th in theories:
        f0, gam = growth_index_of(th, k)
        c = ft.chi2_universe(th.E, th.mu, points, k=k,
                             free_amplitude=free_amplitude)
        rows.append({
            "name": th.name,
            "gamma": gam,
            "gamma_analytic": th.gamma_ref,
            "f_z0": f0,
            "mu_today": float(th.mu(k, 1.0)),
            "scale_dependent": th.scale_dependent,
            "chi2": c["chi2"],
            "chi2_per_dof": c["chi2_per_dof"],
            "amplitude": c["amplitude"],
            "source": th.source,
        })

    # reference is Lambda-CDM (first row)
    chi2_lcdm = rows[0]["chi2"]
    for r in rows:
        r["delta_chi2_vs_lcdm"] = r["chi2"] - chi2_lcdm

    # degeneracy pairs: theories the data cannot separate at 1 sigma
    degenerate = []
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if abs(rows[i]["chi2"] - rows[j]["chi2"]) < 1.0:
                degenerate.append({
                    "a": rows[i]["name"], "b": rows[j]["name"],
                    "delta_chi2": abs(rows[i]["chi2"] - rows[j]["chi2"]),
                })
    return {"k_h_mpc": k, "free_amplitude": free_amplitude,
            "rows": rows, "degenerate_pairs": degenerate}


def separability(theories=None, k: float = 0.1) -> dict:
    """How far apart are the growth indices, i.e. the shapes?

    Reports the spread of gamma and the widest separation between any two
    theories: a large spread means the shape alone can distinguish them.
    """
    if theories is None:
        theories = mg.family()
    gammas = {}
    for th in theories:
        _, g = growth_index_of(th, k)
        gammas[th.name] = g
    vals = np.array(list(gammas.values()))
    return {
        "gamma": gammas,
        "gamma_min": float(vals.min()),
        "gamma_max": float(vals.max()),
        "gamma_spread_abs": float(vals.max() - vals.min()),
        "gamma_spread_rel": float((vals.max() - vals.min()) / vals.min()),
    }
