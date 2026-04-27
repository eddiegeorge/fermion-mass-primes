# Fermion Mass Primes

Scripts and data supporting the paper:

**"Fermion Mass Hierarchies from Prime Stability Conditions"**

Eddie George (2026)

**Version:** 1.0
**DOI (this version):** [10.5281/zenodo.19828753](https://doi.org/10.5281/zenodo.19828753)
**DOI (all versions):** [10.5281/zenodo.19346807](https://doi.org/10.5281/zenodo.19346807)
**Supersedes:** v0.9 (March 2026, DOI: 10.5281/zenodo.19346808)
**License:** Creative Commons Attribution 4.0 International (CC-BY 4.0). Copyright © 2026 Eddie George.

## Contents

- `scripts/` — Python implementations of the four tests cited in the paper
- `logs/` — output logs of the runs cited in the paper
- `results/` — supporting analysis documents
- `README.md` — this file

## Running the scripts

Requirements: Python 3.10+ with `sympy` (for `nextprime`).

To reproduce the results cited in the paper, from the repository root:

```
cd scripts
python fermion_mass_primes_test11.py
python fermion_mass_primes_test12.py
python fermion_mass_primes_test14.py
python fermion_mass_primes_test15.py
```

Each script prints results to standard output. To capture the output in the `logs/` format used in this repo:

```
python fermion_mass_primes_test11.py | tee ../logs/test11_$(date +%Y%m%d_%H%M%S).log
```

Runtimes (approximate):

- test11, test14: sub-second
- test15: ~4 seconds
- test12: ~7–8 minutes for 1,000,000 trials (configurable in script; reduce trial count for faster smoke tests)

## Scripts

### fermion_mass_primes_test11.py

**Multi-scale single-prime substrate test (minimal grammar)**

Purpose:
Tests whether fermion mass ratios are best described as a single odd prime, optionally multiplied by a power of 2. This isolates the hypothesis that ratios are either bare prime modes or substrate-dressed versions of them.

Physical motivation:
If fermion masses have the form m = 2^a × (mode), then ratios take the form 2^(a₁ − a₂) × p.

* If a₁ = a₂, the substrate cancels → pure odd prime
* If a₁ ≠ a₂, the ratio carries a residual power of 2

Grammar:
2^a × p

* a ∈ {0, 1, 2, 3}
* p = single odd prime
* no composite products
* no prime powers
* no denominators

Scope:

* 6 within-family ratios only
  (mc/mu, mt/mc, ms/md, mb/ms, mmu/me, mtau/mmu)
* evaluated at 10 energy scales
  (PDG mixed + 9 Antusch scales)

Method:

* For each ratio at each scale:

  * find best-fit 2^a × p
  * compute residual
* no enumeration over prime sets; direct nearest-fit search

Output:

* best-fit (a, p) per ratio and scale
* residuals
* frequency of a = 0 vs a > 0
* stability of prime assignment across scales

Key questions:

* Do leptonic ratios resolve as pure primes (a = 0)?
* Do quark ratios require substrate dressing (a > 0)?
* Is a single-prime model sufficient without composite products?

### fermion_mass_primes_test12.py

**Hierarchy-matched single-prime substrate null test (multi-scale comparison)**

Purpose:
Tests whether the proximity of fermion mass ratios to single-prime substrate forms 2^a × p (test11) is statistically unusual, by comparing against hierarchy-matched random ratios.

Scope:

* 6 within-family ratios only
  (mc/mu, mt/mc, ms/md, mb/ms, mmu/me, mtau/mmu)
* Real data evaluated at all 10 energy scales
  (PDG mixed + 9 Antusch scales)
* Random trials evaluated once per trial (no scale dependence)

Null model:

* For each trial, generate six independent ratios:
  rᵢ ∼ LogUniform(8.4, 1179)
* The range [8.4, 1179] is an envelope spanning roughly half the minimum to twice the maximum of the observed ratios.

Scoring (test11 metric):

* For each ratio, find best-fit candidate of the form:
  2^a × p
  where:

  * a ∈ {0, 1, 2, 3}
  * p = single odd prime
* Compute residual:
  |rᵢ − 2^a × p| / rᵢ × 100
* Aggregate per trial: mean residual across 6 ratios

Procedure:

* For each real scale:

  * compute mean residual using test11 best fits
* For each of N trials (e.g. 1000 or 1,000,000):

  * generate random six-ratio set
  * compute mean residual using same 2^a × p fitting
* Compare real vs random distribution

Output:

* Per scale:

  * real mean residual
  * percentile within random distribution
* Random distribution summary:

  * mean, median, percentiles, min, max
* Distinct prime usage (real vs random)

Key questions:

* Are real fermion ratios closer to 2^a × p than hierarchy-matched random ratios?
* Does the signal persist across all energy scales?
* Do real ratios exhibit a stronger preference for low substrate powers (a = 0, 2) than random data?

Interpretation:

* If real residuals lie consistently in the lower tail of the random distribution, this supports a structured single-prime substrate model.
* If real and random overlap substantially, the single-prime model may not capture unique structure.

### fermion_mass_primes_test14.py

**Cross-family same-generation boundary test**

#### Purpose

Tests whether the test11 grammar,

R ≈ 2^a × p

describes same-generation cross-family mass ratios as well as it describes within-family ratios.

This is a boundary test, not a grammar extension.

#### Grammar

Identical to test11:

2^a × p

where:

* a ∈ {0, 1, 2, 3}
* p = single odd prime
* no composite products
* no prime powers
* no denominators

The smallest grammar candidate is:

3

so ratios below 3 are flagged as below-floor.

#### Scope

Nine same-generation cross-family ratios:

d/u, c/s, t/b, u/e, c/μ, t/τ, d/e, μ/s, b/τ

Evaluated at 10 energy scales:

* PDG mixed
* 9 Antusch common scales, identical to test11

Neutrino sector excluded.

#### Orientation

For each cell, compute the heavier/lighter ratio:

R = max(mᵢ, mⱼ) / min(mᵢ, mⱼ)

Record the orientation.

If orientation swaps between scales, flag the pair as a unity-crossing.

#### Method

For each ratio at each scale:

* compute R
* flag whether R < 3
* find best-fit 2^a × p
* compute residual:

|R − 2^a × p| / R × 100

#### Output

Per ratio and scale:

* R
* orientation
* below-floor flag
* best (a, p)
* residual

Per ratio:

* mean residual across scales
* number of below-floor cells
* prime stability: locked / walking
* substrate stability: whether a is stable or changes
* unity-crossing flag

Per scale:

* mean cross-family residual
* number of below-floor ratios

Comparison to test11:

* within-family mean residual per scale
* cross-family mean residual per scale
* overall within-family mean residual
* overall cross-family mean residual
* below-floor count for cross-family ratios

#### Key questions

* Do cross-family residuals look like within-family residuals?
* Are cross-family ratios often below the grammar floor?
* Do any cross-family ratios show locked or walking behavior like test11?
* Are apparent good fits clustered by ratio, or just isolated hits?
* Does any unity-crossing occur?

#### Interpretation

If cross-family residuals are worse, or many ratios fall below the grammar floor, the test supports the view that the 2^a × p grammar describes within-family excitation structure rather than cross-family parameter differences.

If cross-family ratios fit comparably well, especially with stable locked or walking behavior, then the family-independence interpretation needs rethinking.

### fermion_mass_primes_test15.py

**Sector prime-set examination**

Purpose:
Examines the common-scale prime vocabularies produced by the test11 grammar, split by charged-fermion sector. No new fits are computed. The goal is to ask whether sector prime sets are structurally informative or mostly artifacts of magnitude range and nearest-candidate fitting.

Source data:
- test11 results at common scales only; PDG mixed excluded
- sector-prime mapping from `test11_sector_audit.md`

Common-scale prime sets:
- Up quark: {37, 71, 73, 127, 139, 151, 271, 277, 281, 307, 503}
- Down quark: {5, 13, 47, 53}
- Charged lepton: {17, 211}

Scope:
- within-family ratios only
- charged fermions only
- unchanged test11 grammar: 2^a × p

---

**T0 — Occupancy table**

For each sector, report each visited prime with:

- count across common-scale cells
- ratio(s) where it appears
- substrate powers (a) used with it
- locked / walking status

Purpose: preserve frequency information before reducing primes to sets.

---

**T1 — Skip-pattern / window sparsity check**

For each sector:

- define the prime window from smallest to largest visited prime
- count all odd primes in that window
- count visited primes
- compute skip fraction:

  skip fraction = 1 − (#visited primes / #available primes in window)

Interpretation:
- low skip fraction = dense coverage
- high skip fraction = sparse/selective coverage

Purpose: ask whether a sector visits most available primes in its magnitude range or only a sparse subset.

---

**T2 — Magnitude-controlled disjointness null**

Tests whether the observed cross-sector disjointness is surprising.

Two nulls:

*A. Window-matched null.*
- For each sector, sample random prime subsets from that sector's own magnitude window
- Preserve each sector's cardinality
- Compute pairwise and total cross-sector overlaps
- Compare actual disjointness to the random distribution

*B. Pooled-window null.*
- Pool: all odd primes ≤ 503 (the largest observed prime); the smallest natural set that contains every observed value
- Sample all sector prime sets from this shared pool, preserving sector cardinalities
- Compute pairwise and total overlaps

Purpose: determine whether disjoint sector vocabularies require explanation, or follow naturally from magnitude separation and small set sizes.

---

**T3 — Algebraic structure scan**

For each sector's prime set, test membership in a fixed list of structurally distinguished prime families:

- modular residue classes:
  - 1 mod 3, 2 mod 3
  - 1 mod 4, 3 mod 4
  - 1 mod 6, 5 mod 6
  - 1 mod 8, 3 mod 8, 5 mod 8, 7 mod 8
  - 1 mod 12, 5 mod 12, 7 mod 12, 11 mod 12
- Eisenstein-prime structure: p ≡ 1 mod 3 or p = 3
- Gaussian-prime structure: p ≡ 3 mod 4
- Sophie Germain primes: 2p+1 also prime
- twin-prime participation: p−2 or p+2 prime
- Mersenne form: 2^n − 1
- Fermat form: 2^(2^n) + 1
- primorial ± 1
- n² + n + 1
- n² + 1
- 4n² + 1

For each property:
- report hit fraction in the actual sector set
- compare against magnitude-matched random prime subsets
- report percentile rank

No multiple-comparisons correction is applied. T3 is diagnostic: isolated hits are not evidence by themselves, but sector-wide patterns or patterns consistent with T0/T1/T2/T4/T5 are worth noting.

---

**T4 — Gap structure within sector**

For each sector with at least 3 primes:

- sort visited primes
- compute consecutive gaps within the visited set
- compare against random subsets of the same cardinality drawn from the same sector window
- report mean gap, gap variance, and percentile ranks versus random subsets

Purpose: ask whether a sector's visited primes are unusually clustered, spread out, or regularly spaced within its magnitude window.

Charged lepton sector has only two primes, so only its single gap is reported descriptively.

---

**T5 — Local candidate-density diagnostic**

For each test11 common-scale cell:

- record target ratio
- record best 2^a × p candidate and residual
- record nearest alternative grammar candidates above and below
- compute local grammar spacing around the target
- compare observed residual to random ratios of similar magnitude

Purpose: ask whether good fits occur in sparse regions of the 2^a × p candidate lattice, or mainly where nearby candidates are dense.

This does not test whether test11 skipped a closer candidate, because test11 always selects the nearest candidate by definition.

---

Output:

For each test, print:
- raw numbers
- compact summary statistics
- plain-language interpretation

Also report:
- common-scale sector prime sets
- common-scale substrate counts by sector
- repeated or locked signatures, especially:
  - lepton bare-prime locks
  - down-sector 2² × 5
  - up-sector walking behavior

---

Key questions:

1. Are sector prime vocabularies genuinely disjoint beyond magnitude effects?
2. Are the vocabularies sparse or dense within their own prime windows?
3. Which primes are repeatedly occupied versus visited once?
4. Do sector prime sets show simple algebraic structure?
5. Are good residuals meaningful, or expected from local candidate density?
6. Do the sector signatures support different internal parameterizations of the same underlying structure?

---

Interpretation:

Test15 is diagnostic, not a new model.

A stronger result would be:
- repeated or locked sector-specific primes
- sparse sector vocabularies
- disjointness not fully explained by magnitude windows
- simple algebraic structure that is sector-wide rather than isolated
- good fits occurring in sparse grammar regions

A weaker result would be:
- disjointness explained by magnitude separation
- apparent algebraic hits isolated to individual primes
- good residuals occurring mainly where grammar candidates are dense

No grammar expansion is introduced.
