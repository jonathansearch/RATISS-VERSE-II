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

**3. The degeneracy map.** The family is extended with DGP (Vainshtein),
f(T) (teleparallel) and Brans-Dicke, every coupling taken from the
literature:

| theory | γ | μ(a=1) | Δχ² (σ₈ fixed) | Δχ² (amp. free) |
|---|---|---|---|---|
| ΛCDM | 0.554 | 1.000 | — | — |
| f(R) `1e-5` | 0.499 | 1.124 | +8.1 | +1.1 |
| f(R) `1e-6` | 0.546 | 1.019 | +1.1 | +0.1 |
| DGP self-acc | **0.668** | 0.727 | **−20.1** | −1.7 |
| f(T) n=0.1 | 0.579 | 0.921 | −2.1 | −0.3 |
| f(T) n=0.3 | 0.692 | 0.661 | −12.6 | −1.1 |
| Brans-Dicke ω=100 | 0.552 | 1.005 | +0.5 | +0.03 |
| Brans-Dicke ω=10 | 0.532 | 1.044 | +4.8 | +0.3 |

Two things stand out.

- **With `σ₈` fixed, the theories separate strongly.** DGP and f(T) *prefer*
  the data over ΛCDM (negative Δχ²) because they weaken gravity (`μ < 1`),
  which mimics the low `S₈` the data want. The growth index spans 38.7 %,
  from f(R) `γ = 0.499` to f(T) `γ = 0.692`.
- **Once the amplitude floats, the map collapses.** 14 of 28 pairs become
  indistinguishable (`Δχ² < 1`); ΛCDM is degenerate with Brans-Dicke ω=100
  at `Δχ² = 0.03`. The `fσ₈` *shape* alone cannot separate these theories —
  the discrimination lives in the amplitude, i.e. in `S₈`, not in `γ`.

That second point is the honest, important negative result: a family with a
wide `γ` spread is still degenerate in `fσ₈(z)` once a single amplitude is
marginalised. Breaking it needs a probe of the *shape* that is not `fσ₈`, or
an independent amplitude.

## Validation

The numerical growth solver is validated against two textbook numbers:

- ΛCDM: `γ = 0.554` vs analytic `6/11 = 0.545` (1.6 %);
- DGP: `γ = 0.668` vs analytic `11/16 = 0.6875` (2.8 %).

That DGP test is the decisive one: it exercises `μ(k,a) ≠ 1`, a
scale-independent coupling, and a non-ΛCDM background all at once.

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
├── mg.py            DGP, f(T), Brans-Dicke and the family registry
├── data.py          the published eBOSS DR16 compilation
├── fit.py           chi^2, the f(R) scan, marginalised amplitude
├── degeneracy.py    the theory map and the degeneracy pairs
└── plots.py         rendering only, no physics
scripts/
├── build.py         recompute everything, seal with SHA-256
└── verify.py        re-verify every hash
```
