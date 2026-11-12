"""
Generalised modified-gravity couplings mu(k,a).

Verse I carried a single genuine modification (f(R)). Verse II needs the
whole family, because the scientific question is now: *can the data tell the
theories apart?* A family where every member predicts the same f sigma_8
would be a degeneracy worth reporting; a family that is separable is a
discovery channel.

Every coupling below is taken from the literature, not invented. The
formula and its source are given at each definition.

    mu(k,a) = G_eff / G   in the quasi-static, sub-horizon limit

so mu = 1 is General Relativity. The background E(a) is supplied separately:
for the designer-background models (f(R), f(T)) it is Lambda-CDM, so the
entire difference lives in the growth -- exactly the discriminator Verse I
identified.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from . import background as bg
from . import constants as C
from . import fr as frmod


# ---------------------------------------------------------------------------
# analytic growth-index benchmarks (used to validate the numerical solver)
# ---------------------------------------------------------------------------
GAMMA_LCDM = 6.0 / 11.0          # 0.5455  (Peebles; standard)
GAMMA_DGP = 11.0 / 16.0          # 0.6875  (Lue, Scoccimarro & Starkman 2004;
                                 #          Linder 2005; Wang 2008 review)


def _dlnH_dlna(E: Callable[[float], float], a: float, eps: float = 1e-4) -> float:
    """dlnH/dlna = Hdot/H^2, by a central difference in N = ln a."""
    N = np.log(a)
    return (np.log(E(np.exp(N + eps))) - np.log(E(np.exp(N - eps)))) / (2 * eps)


# ---------------------------------------------------------------------------
# DGP -- Vainshtein screening, scale-independent
# ---------------------------------------------------------------------------
def rcH(a: float, om: float = C.OMEGA_M) -> float:
    """r_c H(a) for the self-accelerating DGP branch.

    From H^2 - H/r_c = rho_m/3 one gets Omega_m = 1 - 1/(r_c H_0), hence
    r_c H_0 = 1/(1 - Omega_m) and r_c H(a) = E(a)/(1 - Omega_m).
    """
    return bg.E_dgp_selfacc(a, om) / (1.0 - om)


def beta_dgp(a: float, om: float = C.OMEGA_M) -> float:
    """DGP beta function.

        beta(a) = 1 - 2 r_c H (1 + Hdot/(3 H^2)).

    Source: Koyama & Maartens 2006, JCAP 0601, 016, Eq. (18); Wang 2008,
    arXiv:0710.3885. The convention is fixed so that the numerical growth
    index reproduces the analytic DGP value gamma = 11/16 = 0.6875
    (Lue, Scoccimarro & Starkman 2004; Linder 2005).
    """
    E = lambda x: bg.E_dgp_selfacc(x, om)
    return 1.0 - 2.0 * rcH(a, om) * (1.0 + _dlnH_dlna(E, a) / 3.0)


def mu_dgp(a: float, om: float = C.OMEGA_M) -> float:
    """DGP effective Newton constant, G_eff/G = 1 + 1/(3 beta).

    Source: Wang 2008, arXiv:0710.3885; Koyama & Maartens 2006, JCAP 0601, 016.
    """
    b = beta_dgp(a, om)
    return 1.0 + 1.0 / (3.0 * b) if b != 0 else 1.0


def mu_dgp_callable(om: float = C.OMEGA_M) -> Callable[[float, float], float]:
    return lambda k, a: mu_dgp(a, om)


# ---------------------------------------------------------------------------
# f(T) teleparallel gravity -- scale-independent
# ---------------------------------------------------------------------------
def f_T_powerlaw(a: float, n: float, om: float = C.OMEGA_M) -> float:
    """f_T for the power-law model f(T) = alpha (-T)^n.

        f_T = n/(1-2n) (1 - Omega_m0) (H0/H)^{2(1-n)}.

    Source: Wu & Yu 2011, "Growth factor in f(T) gravity",
    arXiv:1010.3512, Eq. (59). Designer background: H is Lambda-CDM.
    """
    E = bg.E_lcdm(a, om)
    return (n / (1.0 - 2.0 * n)) * (1.0 - om) * E ** (-2.0 * (1.0 - n))


def mu_ft(a: float, n: float, om: float = C.OMEGA_M) -> float:
    """f(T) effective Newton constant, G_eff/G = 1/(1 + f_T).

    Source: Wu & Yu 2011, arXiv:1010.3512, Eq. (50).
    """
    return 1.0 / (1.0 + f_T_powerlaw(a, n, om))


def mu_ft_callable(n: float, om: float = C.OMEGA_M) -> Callable[[float, float], float]:
    return lambda k, a: mu_ft(a, n, om)


# ---------------------------------------------------------------------------
# Brans-Dicke -- scale-independent, a clean one-parameter benchmark
# ---------------------------------------------------------------------------
def mu_brans_dicke(omega_BD: float) -> Callable[[float, float], float]:
    """G_eff/G = (4 + 2 omega_BD)/(3 + 2 omega_BD).

    Source: Brans & Dicke 1961; the quasi-static, massless limit (Wang 2008
    review, arXiv:0710.3885). omega_BD -> infinity recovers GR.
    """
    val = (4.0 + 2.0 * omega_BD) / (3.0 + 2.0 * omega_BD)
    return lambda k, a: val


# ---------------------------------------------------------------------------
# phenomenological mu(k,a) -- the model-independent way to report a bound
# ---------------------------------------------------------------------------
def mu_pheno_callable(mu0: float, k_ref: float = 0.1,
                      scale_free: bool = True) -> Callable[[float, float], float]:
    """A minimal two-parameter form, mu = 1 + mu0 (a/1)^s, and its
    scale-free variant where mu does not depend on k.

    This is the standard way to quote a *model-independent* constraint on
    deviations from GR: the data bound mu0, and the bound is then mapped onto
    named theories. Here mu is frozen to its present value for simplicity,
    which is the most conservative (largest-mu) choice at z > 0.
    """
    def mu(k, a):
        return 1.0 + mu0
    return mu


# ---------------------------------------------------------------------------
# the family registry
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Theory:
    name: str
    E: Callable[[float], float]
    mu: Callable[[float, float], float]
    scale_dependent: bool
    gamma_ref: float | None      # analytic growth index, if known
    source: str


def family(om: float = C.OMEGA_M) -> list[Theory]:
    """The theories Verse II compares, with sources."""
    return [
        Theory("Lambda-CDM", bg.E_lcdm, lambda k, a: 1.0, False,
               GAMMA_LCDM, "GR; Peebles 1980"),
        Theory("f(R) |fR0|=1e-5", bg.E_fr_hs,
               frmod.mu_fR_callable(1e-5), True, None,
               "Hu & Sawicki 2007, PRD 76, 064004"),
        Theory("f(R) |fR0|=1e-6", bg.E_fr_hs,
               frmod.mu_fR_callable(1e-6), True, None,
               "Hu & Sawicki 2007, PRD 76, 064004"),
        Theory("DGP self-acc", bg.E_dgp_selfacc, mu_dgp_callable(om), False,
               GAMMA_DGP, "Lue et al 2004, PRD 69, 124015"),
        Theory("f(T) n=0.1", bg.E_ft, mu_ft_callable(0.1, om), False, None,
               "Wu & Yu 2011, arXiv:1010.3512"),
        Theory("f(T) n=0.3", bg.E_ft, mu_ft_callable(0.3, om), False, None,
               "Wu & Yu 2011, arXiv:1010.3512"),
        Theory("Brans-Dicke w=100", bg.E_lcdm, mu_brans_dicke(100.0), False,
               None, "Brans & Dicke 1961"),
        Theory("Brans-Dicke w=10", bg.E_lcdm, mu_brans_dicke(10.0), False,
               None, "Brans & Dicke 1961"),
    ]
