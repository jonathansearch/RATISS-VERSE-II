"""
Figures for Verse II. Rendering only, no physics.
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


C_RED = "#c0392b"
C_BLUE = "#2c6fbb"
C_GREEN = "#27ae60"
C_ORANGE = "#e67e22"
C_GREY = "#7f8c8d"
C_PURPLE = "#8e44ad"


def _save(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig01_data_fs8(points, models, path):
    """The real measurements, and the models that face them.

    models : list of (label, z_array, fs8_array, color)
    """
    fig, ax = plt.subplots(figsize=(8.6, 5.2))

    z = np.array([p.z for p in points])
    y = np.array([p.fs8 for p in points])
    s = np.array([p.sigma for p in points])
    ax.errorbar(z, y, yerr=s, fmt="o", ms=8, lw=0, elinewidth=2,
                capsize=4, color="k", label="eBOSS DR16 (RSD, real data)",
                zorder=5)

    for label, zz, yy, col in models:
        ax.plot(zz, yy, lw=2.4, color=col, label=label)

    ax.set_xlabel(r"redshift  $z$")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title("Growth of structure: measurements vs the RATISS universes")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, path)


def fig02_chi2_scan(scan, path, label="f(R) Hu-Sawicki"):
    """Delta chi^2 versus |f_R0|: the constraint delivered by the data."""
    rows = scan["rows"]
    fR0 = np.array([r["f_R0"] for r in rows])
    dchi2 = np.array([r["delta_chi2_vs_lcdm"] for r in rows])

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.plot(fR0, dchi2, "o-", lw=2.4, ms=8, color=C_BLUE)
    ax.set_xscale("log")
    ax.axhline(0, color="k", lw=1.0, ls=":")
    for level, col, name in ((1.0, C_GREEN, r"$1\sigma$"),
                             (4.0, C_ORANGE, r"$2\sigma$"),
                             (9.0, C_RED, r"$3\sigma$")):
        ax.axhline(level, color=col, lw=1.2, ls="--")
        ax.text(fR0.min(), level + 0.3, name, color=col, fontsize=9)
    ax.set_xlabel(r"$|f_{R0}|$")
    ax.set_ylabel(r"$\Delta\chi^2$ vs $\Lambda$CDM")
    ax.set_title(f"{label}: what eBOSS DR16 says")
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    _save(fig, path)


def fig04_theory_map(rows_fixed, rows_free, path):
    """chi^2 versus growth index gamma, for the whole theory family.

    Left: sigma_8 fixed to Planck (amplitude is a prediction).
    Right: one amplitude marginalised (only the shape is tested).
    """
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), sharey=True)

    pal = [C_RED, C_BLUE, "#1f8fff", C_GREEN, C_ORANGE, "#b8860b",
           C_PURPLE, C_GREY]
    for ax, rows, title in ((axes[0], rows_fixed, r"$\sigma_8$ fixed to Planck"),
                            (axes[1], rows_free, r"one amplitude marginalised")):
        for i, r in enumerate(rows):
            col = pal[i % len(pal)]
            ax.scatter(r["gamma"], r["chi2"], s=130, color=col, zorder=5,
                       edgecolor="k", linewidth=0.6)
            ax.annotate(r["name"], (r["gamma"], r["chi2"]),
                        textcoords="offset points", xytext=(7, 5),
                        fontsize=8, color=col)
        ax.axhline(rows[0]["chi2"], color="k", lw=1.0, ls=":",
                   label=r"$\Lambda$CDM")
        ax.set_xlabel(r"growth index  $\gamma$")
        ax.set_title(title)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$\chi^2$ vs eBOSS DR16")
    axes[0].legend(loc="upper right", fontsize=9)
    fig.suptitle("The degeneracy map: which universes can the data separate?",
                 y=1.02)
    fig.tight_layout()
    _save(fig, path)


def fig05_family_fs8(z_grid, curves, path):
    """f sigma_8(z) for the whole family, against the measurements.

    curves : list of (label, y_array, color)
    """
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    for label, yy, col in curves:
        ax.plot(z_grid, yy, lw=2.2, color=col, label=label)
    ax.set_xlabel(r"redshift  $z$")
    ax.set_ylabel(r"$f\sigma_8(z)$")
    ax.set_title("The whole family of universes: same background, different growth")
    ax.legend(loc="upper right", fontsize=8.5, ncol=2)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, path)


def fig06_gamma_z(z_list, curves, path):
    """gamma(z) for the whole family: the strongest discriminator.

    curves : list of (label, y_array, color)
    """
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    for label, yy, col in curves:
        ax.plot(z_list, yy, lw=2.3, color=col, label=label)
    ax.set_xlabel(r"redshift  $z$")
    ax.set_ylabel(r"growth index  $\gamma(z)=\ln f/\ln\Omega_m$")
    ax.set_title(r"The time shape of $\gamma$: a discriminator the single-scale fit cannot use")
    ax.legend(loc="best", fontsize=8.5, ncol=2)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, path)


def fig07_discriminator(vectors, path):
    """Each theory as a point in the (gamma(0), dgamma/dz) plane, coloured by
    whether it has a scale-dependent coupling."""
    fig, ax = plt.subplots(figsize=(9.0, 6.0))
    pal = [C_RED, C_BLUE, "#1f8fff", C_GREEN, C_ORANGE, "#b8860b",
           C_PURPLE, C_GREY]
    for i, v in enumerate(vectors):
        col = pal[i % len(pal)]
        ax.scatter(v["gamma_0"], v["gamma_slope"], s=150, color=col,
                   zorder=5, edgecolor="k", linewidth=0.6)
        ax.annotate(v["name"], (v["gamma_0"], v["gamma_slope"]),
                    textcoords="offset points", xytext=(8, 5),
                    fontsize=8.5, color=col)
    ax.axhline(0, color="k", lw=1.0, ls=":")
    ax.set_xlabel(r"$\gamma(0)$")
    ax.set_ylabel(r"$d\gamma/dz$")
    ax.set_title(r"Discriminator plane: $\gamma(0)$ alone is degenerate, the slope is not")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, path)


def fig08_environment(fr0_values, env_rows, path):
    """The chameleon screening factor per environment and |f_R0|."""
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    envs = ["void", "group", "galaxy", "cluster"]
    pal = {1e-6: C_GREEN, 1e-5: C_BLUE, 1e-4: C_RED}
    import numpy as np
    for fR0 in fr0_values:
        ys = [next(r["screening_factor"] for r in env_rows
                   if r["f_R0"] == fR0 and r["environment"] == e) for e in envs]
        ax.plot(range(len(envs)), ys, "o-", lw=2.3, ms=8,
                color=pal.get(fR0, "k"), label=rf"$|f_{{R0}}|={fR0:.0e}$")
    ax.set_xticks(range(len(envs)))
    ax.set_xticklabels(envs)
    ax.set_ylabel(r"screening factor  $S$  (1 = unscreened)")
    ax.set_title("Chameleon screening: voids see the fifth force, clusters do not")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, path)


def fig03_fs8_residuals(points, models, path):
    """Residuals of each model against the measurements."""
    z = np.array([p.z for p in points])
    y = np.array([p.fs8 for p in points])
    s = np.array([p.sigma for p in points])

    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    ax.axhline(0, color="k", lw=1.0, ls=":")
    ax.axhspan(-1, 1, color=C_GREY, alpha=0.15, label=r"$\pm1\sigma$")
    for label, zz, yy, col in models:
        yy = np.asarray(yy)
        zz = np.asarray(zz)
        # interpolate model onto the data redshifts
        mz = np.interp(z, zz, yy)
        ax.plot(z, (mz - y) / s, "o-", lw=2.0, ms=6, color=col, label=label)
    ax.set_xlabel(r"redshift  $z$")
    ax.set_ylabel(r"$(f\sigma_8^{\rm model}-f\sigma_8^{\rm data})/\sigma$")
    ax.set_title("Residuals: where each universe fails")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    _save(fig, path)
