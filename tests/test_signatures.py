"""
Tests for the fine signatures (Verse II step 3).

The scientific claim is that the time shape of gamma and the scale dependence
of f sigma_8 separate theories that a single-amplitude f sigma_8 fit cannot.
These tests pin the qualitative structure of that claim.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from v2 import mg
from v2 import signatures as sg


def _by_name(fam):
    return {th.name: th for th in fam}


# ---------------------------------------------------------------------------
# 1. gamma(z)
# ---------------------------------------------------------------------------
def test_gamma_slope_separates_lcdm_from_fr_and_dgp():
    fam = _by_name(mg.family())
    s_lcdm = sg.gamma_slope(fam["Lambda-CDM"])
    s_fr = sg.gamma_slope(fam["f(R) |fR0|=1e-5"])
    s_dgp = sg.gamma_slope(fam["DGP self-acc"])
    # f(R) falls much faster, DGP rises
    assert s_fr < s_lcdm < s_dgp
    assert s_fr < -0.1
    assert s_dgp > 0.0


def test_gamma_slope_is_the_best_discriminator():
    d = sg.discriminability()
    power = d["signature_power"]
    assert power["gamma_slope"]["spread_rel"] > power["gamma_0"]["spread_rel"]
    assert power["gamma_slope"]["spread_rel"] > 5.0


# ---------------------------------------------------------------------------
# 2. scale dependence
# ---------------------------------------------------------------------------
def test_only_fr_is_scale_dependent():
    fam = mg.family()
    for th in fam:
        v = sg.fs8_vs_k(th)
        if th.scale_dependent:
            assert v["k_variation_rel"] > 0.02, th.name
        else:
            assert v["k_variation_rel"] < 1e-6, th.name


def test_fr_scale_dependence_grows_with_fr0():
    """A stronger f(R) is more scale-dependent -- monotone, as expected."""
    v5 = sg.fs8_vs_k(mg.family()[1])["k_variation_rel"]      # 1e-5
    v6 = sg.fs8_vs_k(mg.family()[2])["k_variation_rel"]      # 1e-6
    assert v5 > v6


# ---------------------------------------------------------------------------
# 3. environment
# ---------------------------------------------------------------------------
def test_thin_shell_is_monotone_in_environment_depth():
    """Shallow environments are unscreened, deep ones are screened."""
    S_void = sg.thin_shell_factor(1e-5, sg.PHI_ENVIRONMENTS["void"])
    S_cluster = sg.thin_shell_factor(1e-5, sg.PHI_ENVIRONMENTS["cluster"])
    assert S_void > S_cluster
    assert S_void > 0.9
    assert S_cluster < 0.2


def test_environment_split_is_nonzero_only_for_chameleon():
    split = sg.environment_growth_split(1e-5)
    assert split["split_rel"] > 0.01
    assert split["void_minus_cluster"] > 0.0


def test_thin_shell_increases_with_fr0():
    S_small = sg.thin_shell_factor(1e-7, 1e-5)
    S_large = sg.thin_shell_factor(1e-4, 1e-5)
    assert S_large > S_small


# ---------------------------------------------------------------------------
# the discriminator vector
# ---------------------------------------------------------------------------
def test_discriminator_vector_separates_degenerate_pair():
    """Lambda-CDM and Brans-Dicke w=100 were degenerate in chi^2 (v0.2.0);
    they must still differ in the time shape of gamma."""
    fam = _by_name(mg.family())
    v_l = sg.discriminator_vector(fam["Lambda-CDM"])
    v_bd = sg.discriminator_vector(fam["Brans-Dicke w=100"])
    assert abs(v_l["gamma_slope"] - v_bd["gamma_slope"]) > 1e-3
