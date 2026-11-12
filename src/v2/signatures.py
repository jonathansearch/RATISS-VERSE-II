"""
Fine signatures that might break the degeneracy found in degeneracy.py.

Verse II v0.2.0 showed that the *amplitude* of f sigma_8 separates theories
but its *shape* (one amplitude marginalised) does not: 14 of 28 pairs became
degenerate. This module looks for observables that are not a single-amplitude
f sigma_8, and asks whether any of them separates the family.

Three signatures are computed.

1. gamma(z): the growth index as a *function* of redshift. A theory with the
   same gamma(0) but a different slope is distinguishable by a wide-z survey.

2. f sigma_8(k): the scale dependence. f(R) has a scale-dependent mu, while
   DGP, f(T) and Brans-Dicke are scale-independent. Measuring the same growth
   at two wavenumbers is therefore a clean discriminator -- the one thing the
   single-scale fit cannot see.

3. environment: chameleon screening depends on the local Newtonian potential,
   so a void and a cluster probe different mu. The thin-shell factor is
   computed and turned into an environment-dependent growth rate.
"""

from __future__ import annotations

import numpy as np

from . import background as bg
from . import constants as C
from . import fr as frmod
from . import growth as gw


# ---------------------------------------------------------------------------
# 1. gamma(z)
# ---------------------------------------------------------------------------
def gamma_at(E, mu, z: float, k: float = 0.1) -> float:
    """gamma(z) = ln f / ln Omega_m(a) at one redshift."""
    N = -np.log1p(z)
    _, D, dD, om = gw.growth_factor(E, mu, k, N_fin=N)
    if D == 0 or dD / D <= 0 or not (0.0 < om < 1.0):
        return float("nan")
    f = dD / D
    return float(np.log(f) / np.log(om))


def gamma_curve(th, z_list, k: float = 0.1) -> list[float]:
    return [gamma_at(th.E, th.mu, z, k) for z in z_list]


def gamma_slope(th, k: float = 0.1) -> float:
    """d gamma / dz near z = 0, by a central difference. The *time* shape."""
    h = 0.05
    return (gamma_at(th.E, th.mu, h, k) - gamma_at(th.E, th.mu, -h + 1e-6, k)) / h


# ---------------------------------------------------------------------------
# 2. f sigma_8(k) -- the scale dependence
# ---------------------------------------------------------------------------
def fs8_at_k(th, k: float, z: float = 0.5) -> float:
    rows = gw.f_sigma8(th.E, th.mu, k, redshifts=(z,))
    return rows[0]["f_sigma8"]


def fs8_vs_k(th, k_list=(0.02, 0.05, 0.1, 0.2, 0.4), z: float = 0.5) -> dict:
    vals = [fs8_at_k(th, k, z) for k in k_list]
    vals = np.array(vals)
    return {
        "k": list(k_list),
        "fs8": vals.tolist(),
        "k_variation_rel": float((vals.max() - vals.min()) / vals.mean()),
    }


def scale_dependence_table(theories, k_list=(0.02, 0.05, 0.1, 0.2, 0.4),
                           z: float = 0.5) -> dict:
    """Which theories have a k-dependent f sigma_8, and by how much.

    This is the discriminator the single-scale fit of v0.2.0 could not use.
    """
    out = {}
    for th in theories:
        out[th.name] = {
            "scale_dependent": th.scale_dependent,
            **fs8_vs_k(th, k_list, z),
        }
    return {"z": z, "k_list": list(k_list), "theories": out}


# ---------------------------------------------------------------------------
# 3. environment-dependent screening
# ---------------------------------------------------------------------------
# Representative Newtonian potential depths |Phi_N| = GM/(R c^2):
#   void        ~ 1e-8   (large, empty, shallow)
#   field/group ~ 1e-6
#   galaxy      ~ 1e-5
#   cluster     ~ 1e-4   (deep, dense)
PHI_ENVIRONMENTS = {
    "void": 1e-8,
    "group": 1e-6,
    "galaxy": 1e-5,
    "cluster": 1e-4,
}


def thin_shell_factor(fR0: float, phi_N: float) -> float:
    """Chameleon screening of the fifth force in an environment of depth phi_N.

    The thin-shell condition (Khoury & Weltman 2004; Hu & Sawicki 2007) says
    an object is unscreened when the background field amplitude exceeds its
    potential depth. We use the standard smooth interpolation

        S = x/(1+x),    x = |f_R0| / ((3/2) phi_N),

    so S -> 1 in a shallow, unscreened environment (void) and S -> 0 in a
    deep, screened one (cluster). This is a controlled approximation, not a
    full chameleon solve; it captures the monotonic dependence on |f_R0| and
    phi_N that the real solution has.
    """
    x = abs(fR0) / (1.5 * phi_N)
    return float(x / (1.0 + x))


def mu_environment(a: float, fR0: float, environment: str, k: float = 0.1) -> float:
    """Environment-dependent effective Newton constant.

    The bare chameleon enhancement of fr.mu_fR is multiplied by the thin-shell
    factor, so the enhancement survives in voids and is erased in clusters.
    """
    phi_N = PHI_ENVIRONMENTS[environment]
    bare = float(frmod.mu_fR(k, a, fR0))
    return 1.0 + (bare - 1.0) * thin_shell_factor(fR0, phi_N)


def mu_environment_callable(fR0: float, environment: str, k: float = 0.1):
    return lambda kk, aa: mu_environment(aa, fR0, environment, kk)


def environment_table(fR0_values=(1e-6, 1e-5, 1e-4)) -> dict:
    """Screening factor and mu enhancement per environment and |f_R0|."""
    rows = []
    for fR0 in fR0_values:
        for env, phi in PHI_ENVIRONMENTS.items():
            S = thin_shell_factor(fR0, phi)
            rows.append({
                "f_R0": fR0,
                "environment": env,
                "phi_N": phi,
                "screening_factor": S,
                "mu_today": mu_environment(1.0, fR0, env),
                "enhancement_pct": 100.0 * (mu_environment(1.0, fR0, env) - 1.0),
            })
    return {"rows": rows}


def environment_growth_split(fR0: float, z: float = 0.5) -> dict:
    """How much does the growth rate differ between a void and a cluster?

    A non-zero split is a *screening* signature: it says the fifth force is
    environment-dependent, which no scale-independent theory reproduces.
    """
    out = {}
    for env in PHI_ENVIRONMENTS:
        mu = mu_environment_callable(fR0, env)
        rows = gw.f_sigma8(bg.E_lcdm, mu, 0.1, redshifts=(z,))
        out[env] = rows[0]["f_sigma8"]
    vals = list(out.values())
    out["void_minus_cluster"] = out["void"] - out["cluster"]
    out["split_rel"] = (max(vals) - min(vals)) / np.mean(vals)
    out["f_R0"] = fR0
    out["z"] = z
    return out


# ---------------------------------------------------------------------------
# 4. the discriminator vector -- which signature breaks which degeneracy?
# ---------------------------------------------------------------------------
def discriminator_vector(th, k: float = 0.1) -> dict:
    """A four-component fingerprint per theory, each component a *different*
    observable, so that theories degenerate in one can differ in another.

        (gamma(0), dgamma/dz, k-variation of f sigma_8, void-cluster split)
    """
    gam0 = gamma_at(th.E, th.mu, 0.0, k)
    slope = gamma_slope(th, k)
    kv = fs8_vs_k(th, (0.02, 0.1, 0.4), z=0.5)["k_variation_rel"]
    # the environment split only exists for a chameleon-like mu; for the
    # scale-independent theories it is zero by construction
    split = 0.0
    if th.scale_dependent:
        fR0 = _fR0_of(th.name)
        if fR0 is not None:
            split = environment_growth_split(fR0)["split_rel"]
    return {
        "name": th.name,
        "gamma_0": gam0,
        "gamma_slope": slope,
        "k_variation_rel": kv,
        "environment_split_rel": split,
    }


def _fR0_of(name: str) -> float | None:
    import re
    m = re.search(r"\|fR0\|=([0-9.eE+-]+)", name)
    return float(m.group(1)) if m else None


def discriminability(theories=None, k: float = 0.1) -> dict:
    """Which signature separates the family best?

    For each signature we normalise its spread across the family and report
    it. A signature with a large spread is the one to build a survey around.
    """
    if theories is None:
        from . import mg
        theories = mg.family()
    vecs = [discriminator_vector(th, k) for th in theories]
    import numpy as np
    out = {}
    for key in ("gamma_0", "gamma_slope", "k_variation_rel",
                "environment_split_rel"):
        vals = np.array([v[key] for v in vecs])
        out[key] = {
            "min": float(vals.min()),
            "max": float(vals.max()),
            "spread": float(vals.max() - vals.min()),
            "spread_rel": float((vals.max() - vals.min())
                                / (abs(vals.mean()) + 1e-30)),
        }
    return {"k_h_mpc": k, "vectors": vecs, "signature_power": out}
