# RATISS Verse II — confronting the registry with real data

Verse I (`RATISS-VERSE-I`) built a registry of ten universes and proved one
negative result: the dimensionless age `H₀t₀` cannot separate them. Only the
growth of structure — the growth index `γ`, and `fσ₈(z)` — carries the
signal.

Verse II takes the next step. It takes **published** growth-rate measurements
and asks which universes survive.

## The data

`eBOSS DR16` consensus, RSD-only `fσ₈(z)`, six tracers (Table 3 of
eBOSS Collaboration 2021, PRD 103, 083533, arXiv:2007.08991):

| tracer | z | fσ₈ |
|---|---|---|
| MGS | 0.15 | 0.530 ± 0.160 |
| BOSS gal | 0.38 | 0.497 ± 0.045 |
| BOSS gal | 0.51 | 0.459 ± 0.038 |
| eBOSS LRG | 0.70 | 0.473 ± 0.041 |
| eBOSS ELG | 0.85 | 0.315 ± 0.095 |
| eBOSS QSO | 1.48 | 0.462 ± 0.045 |

Every value is transcribed from the paper and checked by a test.

## The results

**1. ΛCDM with the Planck `σ₈` is in tension with the data.**
Fixing `σ₈ = 0.8111` gives `χ²/dof = 4.47`. Letting one amplitude float gives
`χ²/dof = 1.17` with a best-fit amplitude of **0.833** — i.e. the data prefer
a growth amplitude about 17 % *below* Planck. This is the well-known
`σ₈`/`S₈` tension, reproduced here from the fit rather than assumed.

**2. The data constrain f(R).** Scanning the Hu–Sawicki chameleon coupling
(`μ(k,a) = 1 + (1/3) x/(1+x)`, `x = (k/aM(a))²`) at `k = 0.1 h/Mpc`:

| `|f_R0|` | Compton wavelength | `Δχ²` vs ΛCDM | verdict |
|---|---|---|---|
| `3e-4` | 42 Mpc/h | +31.9 | excluded |
| `1e-4` | 24 Mpc/h | +25.6 | excluded |
| `3e-5` | 13 Mpc/h | +16.1 | excluded |
| `1e-5` | 7.7 Mpc/h | +8.1 | disfavoured |
| `3e-6` | 4.2 Mpc/h | +3.0 | allowed |
| `1e-6` | 2.4 Mpc/h | +1.1 | allowed |
| `1e-7` | 0.8 Mpc/h | +0.1 | invisible |

The 2σ bound is `|f_R0| < 3e-6`: **f(R) must be strongly screened** to survive
the observed growth of structure. This is the first quantitative constraint
this line of work produces.

## What this is, and is not

- The background is ΛCDM for every `|f_R0|` (designer background), so the
  entire constraint comes from the growth sector — exactly the discriminator
  Verse I identified.
- `σ₈` is an input, not a prediction (Verse I, N6). The scan is reported both
  with `σ₈` fixed and with one amplitude marginalised, because the two answer
  different questions.
- `k = 0.1 h/Mpc` is a representative scale, inside the DESI/Euclid window.
- The physics core (`background`, `growth`, `fr`) is shared with Verse I
  verbatim, so the two repositories cannot drift apart.

## Reproduce

```
pip install -r requirements.txt
python scripts/build.py     # ~2 s, seals 6 artefacts with SHA-256
python scripts/verify.py    # re-verifies every hash
python -m pytest tests/ -q  # 6 tests
```

## Layout

```
src/v2/
├── constants.py     Planck 2018 parameters (shared with Verse I)
├── background.py    E(a) for every universe (shared with Verse I)
├── growth.py        the growth ODE and f sigma_8 (shared with Verse I)
├── fr.py            Hu-Sawicki chameleon coupling (shared with Verse I)
├── data.py          the published eBOSS DR16 compilation
├── fit.py           chi^2, the f(R) scan, marginalised amplitude
└── plots.py         rendering only, no physics
scripts/
├── build.py         recompute everything, seal with SHA-256
└── verify.py        re-verify every hash
```
