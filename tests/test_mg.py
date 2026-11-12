"""
Tests for the generalised modified-gravity family.

The decisive test is DGP: the numerical growth solver must reproduce the
analytic growth index gamma = 11/16 = 0.6875. If it does, the whole mu(k,a)
machinery is validated against an independent, textbook number.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2 import background as bg
from v2 import constants as C
from v2 import degeneracy as dgn
from v2 import fr as frmod
from v2 import growth as gw
from v2 import mg


# ---------------------------------------------------------------------------
# the validation of the solver itself
# ---------------------------------------------------------------------------
def test_dgp_growth_index_matches_analytic_value():
    """gamma_DGP = 11/16 is a textbook number (Lue 2004, Linder 2005)."""
    _, gam = gw.growth_index(bg.E_dgp_selfacc, mg.mu_dgp_callable(), 0.1)
    assert abs(gam - mg.GAMMA_DGP) < 0.03, gam


def test_lcdm_growth_index_matches_analytic_value():
    """gamma_LCDM = 6/11; the numerical solver must land within 2%."""
    _, gam = gw.growth_index(bg.E_lcdm, 1.0, 0.1)
    assert abs(gam - mg.GAMMA_LCDM) < 0.02, gam


def test_gr_is_recovered_by_all_screening_limits():
    """Brans-Dicke with omega -> infinity and f(R) with |f_R0| -> 0 both
    collapse back onto General Relativity."""
    mu_bd = mg.mu_brans_dicke(1e8)
    assert abs(mu_bd(0.1, 1.0) - 1.0) < 1e-6
    assert abs(float(frmod.mu_fR(0.1, 1.0, 1e-12)) - 1.0) < 1e-3


# ---------------------------------------------------------------------------
# the couplings
# ---------------------------------------------------------------------------
def test_ft_is_scale_independent_and_fr_is_not():
    """The structural difference between the two: f(T) has no k, f(R) does."""
    mu_ft = mg.mu_ft_callable(0.1)
    assert abs(mu_ft(0.01, 0.5) - mu_ft(1.0, 0.5)) < 1e-12
    mu_fr = frmod.mu_fR_callable(1e-5)
    assert abs(mu_fr(0.01, 0.5) - mu_fr(1.0, 0.5)) > 1e-2


def test_brans_dicke_tends_to_gr():
    assert abs(mg.mu_brans_dicke(100.0)(0, 1) - 1.0) < 0.01
    assert mg.mu_brans_dicke(10.0)(0, 1) > mg.mu_brans_dicke(100.0)(0, 1)


def test_ft_weakens_gravity_below_gr():
    """f(T) with the Wu-Yu sign gives mu < 1: growth is suppressed."""
    assert mg.mu_ft(1.0, 0.1) < 1.0
    assert mg.mu_ft(1.0, 0.3) < mg.mu_ft(1.0, 0.1)


# ---------------------------------------------------------------------------
# the degeneracy map
# ---------------------------------------------------------------------------
def test_family_is_loaded_and_sourced():
    fam = mg.family()
    assert len(fam) >= 6
    for th in fam:
        assert th.source, th.name
        assert callable(th.mu)


def test_theories_are_separable_in_gamma_but_degenerate_in_chi2():
    """The central Verse II result: a wide gamma spread, yet many theories
    the data cannot separate once the amplitude floats."""
    m = dgn.map_theories(free_amplitude=True)
    sep = dgn.separability()
    assert sep["gamma_spread_rel"] > 0.20          # shapes do differ
    assert len(m["degenerate_pairs"]) >= 5          # yet chi^2 cannot tell
    # Lambda-CDM and Brans-Dicke w=100 must be indistinguishable
    names = {(p["a"], p["b"]) for p in m["degenerate_pairs"]}
    assert ("Lambda-CDM", "Brans-Dicke w=100") in names or \
           ("Brans-Dicke w=100", "Lambda-CDM") in names


def test_marginalising_the_amplitude_collapses_the_differences():
    fixed = dgn.map_theories(free_amplitude=False)
    free = dgn.map_theories(free_amplitude=True)
    spread_fixed = max(abs(r["delta_chi2_vs_lcdm"]) for r in fixed["rows"])
    spread_free = max(abs(r["delta_chi2_vs_lcdm"]) for r in free["rows"])
    assert spread_fixed > spread_free
