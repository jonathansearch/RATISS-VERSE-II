"""
Confront the RATISS universes with published growth-rate data.

For a universe defined by ``(E(a), mu(k,a))`` the growth solver predicts

    f sigma_8(z) = sigma_8 * f(z) * D(z)/D(0),     f = dlnD/dlna,

where ``sigma_8`` is an amplitude. The repository's own negative result N6
says ``sigma_8`` is *not* predicted, so the fit either fixes it to the Planck
value or marginalises over it analytically (one free amplitude per model).
Both are reported; the difference is the honest part of the exercise.

The figure of merit is chi^2 against the eBOSS DR16 compilation, with the
Gaussian errors quoted in the paper.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from . import background as bg
from . import constants as C
from . import growth as gw
from . import fr as frmod


def predict_fs8(E, mu=1.0, k: float = 0.0, redshifts=(0.0,)) -> np.ndarray:
    """Predicted f sigma_8(z) with sigma_8 fixed to the Planck value."""
    rows = gw.f_sigma8(E, mu, k, redshifts=tuple(redshifts))
    return np.array([r["f_sigma8"] for r in rows])


def predict_shape(E, mu=1.0, k: float = 0.0, redshifts=(0.0,)) -> np.ndarray:
    """The *shape* f(z) D(z)/D(0): what a free amplitude multiplies."""
    rows = gw.f_sigma8(E, mu, k, redshifts=tuple(redshifts))
    return np.array([r["f_sigma8"] for r in rows]) / C.SIGMA_8


def _chi2_terms(m: np.ndarray, d: np.ndarray, s: np.ndarray,
                free_amplitude: bool) -> tuple[float, float]:
    """Return (chi^2, best-fit amplitude). Amplitude 1 if not free."""
    if not free_amplitude:
        return float(np.sum(((m - d) / s) ** 2)), 1.0
    w = 1.0 / s ** 2
    amp = float(np.sum(w * m * d) / np.sum(w * m * m))
    return float(np.sum(((amp * m - d) / s) ** 2)), amp


def chi2_universe(E, mu, points, k: float = 0.0,
                  free_amplitude: bool = False) -> dict:
    """chi^2 of one universe against a list of GrowthPoint."""
    z = np.array([p.z for p in points])
    d = np.array([p.fs8 for p in points])
    s = np.array([p.sigma for p in points])
    m = predict_shape(E, mu, k, redshifts=tuple(z))
    chi2, amp = _chi2_terms(m, d, s, free_amplitude)
    dof = len(points) - (1 if free_amplitude else 0)
    return {
        "chi2": chi2,
        "dof": dof,
        "chi2_per_dof": chi2 / dof,
        "amplitude": amp,
        "n_points": len(points),
    }


# ---------------------------------------------------------------------------
# the f(R) chameleon scan -- the scientific deliverable of Verse II
# ---------------------------------------------------------------------------
def scan_fr0(points, f_r0_values=None, k: float = 0.1,
             free_amplitude: bool = False) -> dict:
    """chi^2 of Hu--Sawicki f(R) versus |f_R0|, plus the Lambda-CDM reference.

    The background is Lambda-CDM for every |f_R0| (designer background), so
    the only thing that changes along the scan is the growth -- exactly the
    discriminating power Verse I identified and could not yet quantify.
    """
    if f_r0_values is None:
        f_r0_values = [3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 1e-7]

    ref = chi2_universe(bg.E_lcdm, 1.0, points, k=k,
                        free_amplitude=free_amplitude)
    rows = []
    for fR0 in f_r0_values:
        mu = frmod.mu_fR_callable(fR0)
        c = chi2_universe(bg.E_lcdm, mu, points, k=k,
                          free_amplitude=free_amplitude)
        rows.append({
            "f_R0": fR0,
            "compton_wavelength_mpc_h": frmod.compton_wavelength_today(fR0),
            "mu_at_k": float(frmod.mu_fR(k, 1.0, fR0)),
            **c,
            "delta_chi2_vs_lcdm": c["chi2"] - ref["chi2"],
        })
    return {
        "k_h_mpc": k,
        "free_amplitude": free_amplitude,
        "lcdm": ref,
        "rows": rows,
    }


def best_fit_amplitude(points, E, mu=1.0, k: float = 0.0) -> float:
    """The single best-fit sigma_8 amplitude for a model."""
    z = np.array([p.z for p in points])
    d = np.array([p.fs8 for p in points])
    s = np.array([p.sigma for p in points])
    m = predict_shape(E, mu, k, redshifts=tuple(z))
    w = 1.0 / s ** 2
    return float(np.sum(w * m * d) / np.sum(w * m * m))
