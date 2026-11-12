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
