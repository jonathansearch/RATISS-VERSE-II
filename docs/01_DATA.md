# The data

Verse II uses one primary dataset: the **eBOSS DR16 consensus** growth-rate
measurement.

## Provenance

> eBOSS Collaboration 2021, *Completed SDSS-IV extended Baryon Oscillation
> Spectroscopic Survey: Cosmological implications from two decades of
> spectroscopic surveys*, Phys. Rev. D **103**, 083533, arXiv:2007.08991.

We use the **RSD-only** `fσ₈(z)` values (Table 3), the standard compilation
used throughout the modified-gravity literature:

| tracer | z_eff | fσ₈ | σ |
|---|---|---|---|
| MGS | 0.15 | 0.530 | 0.160 |
| BOSS galaxy | 0.38 | 0.497 | 0.045 |
| BOSS galaxy | 0.51 | 0.459 | 0.038 |
| eBOSS LRG | 0.70 | 0.473 | 0.041 |
| eBOSS ELG | 0.85 | 0.315 | 0.095 |
| eBOSS QSO | 1.48 | 0.462 | 0.045 |

These are transcribed by hand from the published table; a test
(`test_eboss_compilation_matches_published_table`) pins every value so that a
typo cannot enter silently.

## What `fσ₈` is

`fσ₈(z) = f(z) · σ₈ · D(z)/D(0)`, with `f = dlnD/dlna` the logarithmic growth
rate. It is the combination that redshift-space distortions measure directly:
it depends on how fast structure grows, not only on how the universe expands.

This is precisely why it discriminates. Verse I showed that `H₀t₀` is the same
for many of the registry universes; `fσ₈(z)` is not.

## Error model

The quoted `σ` are treated as 1σ Gaussian and uncorrelated, as in the
original RSD-only analysis. No covariance between tracers is used; this is
stated so the constraint is not over-sold.
