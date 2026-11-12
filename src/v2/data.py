"""
Published growth-rate measurements used by Verse II.

Every number here is transcribed from a peer-reviewed paper, not simulated.
The provenance is recorded per dataset so that any result can be traced back
to its source.

Primary dataset
---------------
eBOSS DR16 consensus, RSD-only f sigma_8(z), Table 3 of

    eBOSS Collaboration 2021, "Completed SDSS-IV extended Baryon
    Oscillation Spectroscopic Survey: Cosmological implications from two
    decades of spectroscopic surveys", Phys. Rev. D 103, 083533,
    arXiv:2007.08991.

The six tracers and their effective redshifts are the standard compilation
used in the modified-gravity literature (MGS, BOSS galaxies, eBOSS LRG/ELG/
QSO). Values are the RSD-only f sigma_8 with 1-sigma Gaussian errors.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GrowthPoint:
    tracer: str
    z: float
    fs8: float
    sigma: float
    ref: str = "eBOSS DR16 (2021), PRD 103, 083533, arXiv:2007.08991"


# RSD-only f sigma_8(z), eBOSS DR16 consensus (Table 3).
EBOSS_DR16 = [
    GrowthPoint("MGS",        0.15, 0.53,  0.16),
    GrowthPoint("BOSS gal",   0.38, 0.497, 0.045),
    GrowthPoint("BOSS gal",   0.51, 0.459, 0.038),
    GrowthPoint("eBOSS LRG",  0.70, 0.473, 0.041),
    GrowthPoint("eBOSS ELG",  0.85, 0.315, 0.095),
    GrowthPoint("eBOSS QSO",  1.48, 0.462, 0.045),
]


def eboss_dr16() -> list[GrowthPoint]:
    """The eBOSS DR16 growth compilation (six points)."""
    return list(EBOSS_DR16)


def as_arrays(points: list[GrowthPoint]):
    """(z, fs8, sigma) as plain lists, for the fitter."""
    z = [p.z for p in points]
    y = [p.fs8 for p in points]
    s = [p.sigma for p in points]
    return z, y, s
