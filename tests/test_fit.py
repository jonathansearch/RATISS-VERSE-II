"""
Physics and statistics tests for Verse II.

Each test pins a number or an identity that must not drift. The data itself
is checked against the published table; the fitter is checked against an
independent brute-force computation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2 import background as bg
from v2 import constants as C
from v2 import data as dat
from v2 import fit as ft
from v2 import fr as frmod
from v2 import growth as gw


# ---------------------------------------------------------------------------
# the data
# ---------------------------------------------------------------------------
def test_eboss_compilation_matches_published_table():
    pts = dat.eboss_dr16()
    assert len(pts) == 6
    expected = {0.15: 0.53, 0.38: 0.497, 0.51: 0.459,
                0.70: 0.473, 0.85: 0.315, 1.48: 0.462}
    for p in pts:
        assert abs(p.fs8 - expected[p.z]) < 1e-9, (p.z, p.fs8)
        assert p.sigma > 0


# ---------------------------------------------------------------------------
# the fitter
# ---------------------------------------------------------------------------
def test_chi2_matches_bruteforce():
    pts = dat.eboss_dr16()
    z = np.array([p.z for p in pts])
    d = np.array([p.fs8 for p in pts])
    s = np.array([p.sigma for p in pts])
    shape = ft.predict_shape(bg.E_lcdm, 1.0, 0.1, redshifts=tuple(z))
    # brute force over amplitudes
    amps = np.linspace(0.5, 1.2, 20001)
    chi2 = [np.sum(((a * shape - d) / s) ** 2) for a in amps]
    i = int(np.argmin(chi2))
    got = ft.chi2_universe(bg.E_lcdm, 1.0, pts, k=0.1, free_amplitude=True)
    assert abs(got["chi2"] - chi2[i]) < 1e-2
    assert abs(got["amplitude"] - amps[i]) < 1e-3


def test_lcdm_amplitude_free_is_a_good_fit_but_planck_sigma8_is_tension():
    """The sigma_8 tension is visible in the data, not assumed."""
    pts = dat.eboss_dr16()
    fixed = ft.chi2_universe(bg.E_lcdm, 1.0, pts, k=0.1, free_amplitude=False)
    free = ft.chi2_universe(bg.E_lcdm, 1.0, pts, k=0.1, free_amplitude=True)
    assert free["chi2"] < fixed["chi2"]
    assert free["chi2_per_dof"] < 2.0
    assert fixed["chi2_per_dof"] > 3.0
    # the best-fit amplitude is below the Planck value: the S8 tension
    assert free["amplitude"] < 0.9


# ---------------------------------------------------------------------------
# the f(R) constraint
# ---------------------------------------------------------------------------
def test_strong_fr_is_disfavoured_and_screened_fr_is_not():
    pts = dat.eboss_dr16()
    scan = ft.scan_fr0(pts, [3e-4, 1e-5, 1e-7], k=0.1, free_amplitude=False)
    rows = {r["f_R0"]: r for r in scan["rows"]}
    assert rows[3e-4]["delta_chi2_vs_lcdm"] > 9.0      # excluded at >3 sigma
    assert rows[1e-7]["delta_chi2_vs_lcdm"] < 1.0      # invisible
    d = [rows[k]["delta_chi2_vs_lcdm"] for k in (3e-4, 1e-5, 1e-7)]
    assert d[0] > d[1] > d[2]                          # monotone in screening


def test_fr_scan_uses_lcdm_background():
    """The whole point: the background is identical, only growth changes."""
    pts = dat.eboss_dr16()
    scan = ft.scan_fr0(pts, [1e-5], k=0.1, free_amplitude=True)
    assert abs(bg.hubble_time_times_h0(bg.E_fr_hs)[0]
               - bg.hubble_time_times_h0(bg.E_lcdm)[0]) < 1e-12
    assert scan["rows"][0]["delta_chi2_vs_lcdm"] > 0.5


def test_amplitude_free_scan_is_gentler_than_fixed():
    """Marginalising over sigma_8 weakens (correctly) the constraint."""
    pts = dat.eboss_dr16()
    fixed = ft.scan_fr0(pts, [1e-4], k=0.1, free_amplitude=False)["rows"][0]
    free = ft.scan_fr0(pts, [1e-4], k=0.1, free_amplitude=True)["rows"][0]
    assert fixed["delta_chi2_vs_lcdm"] > free["delta_chi2_vs_lcdm"]
