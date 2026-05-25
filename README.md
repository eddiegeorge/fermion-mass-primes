# Fermion Mass Primes

Scripts and data supporting the paper:

**"Fermion Mass Hierarchies from Prime Stability Conditions"**

Eddie George (2026)

**Version:** 2.0
**DOI (this version):** [10.5281/zenodo.20377363](https://doi.org/10.5281/zenodo.20377363)
**DOI (all versions):** [10.5281/zenodo.19346807](https://doi.org/10.5281/zenodo.19346807)
**Supersedes:** v1.0 (April 2026, DOI: 10.5281/zenodo.19828753), v0.9 (March 2026, DOI: 10.5281/zenodo.19346808)
**License:** Creative Commons Attribution 4.0 International (CC-BY 4.0). Copyright © 2026 Eddie George.

## Contents

- `scripts/` — Python implementations of the nine tests cited in the paper
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
python fermion_mass_primes_test16.py
python fermion_mass_primes_test17.py
python fermion_mass_primes_test18.py
python fermion_mass_primes_test19.py
python fermion_mass_primes_test21.py
```

Each script prints results to standard output. To capture the output in the `logs/` format used in this repo:

```
python fermion_mass_primes_test11.py | tee ../logs/test11_$(date +%Y%m%d_%H%M%S).log
```

Runtimes (approximate):

- test11, test14, test16: sub-second
- test15, test17: ~4 seconds
- test19: ~1 second
- test12: ~7–8 minutes for 1,000,000 trials (configurable in script; reduce trial count for faster smoke tests)
- test18: ~10–15 minutes for 1,000,000 trials per scale
- test21: ~5 minutes for 1,000,000 trials per random null plus ~30 seconds for exact permutation enumeration

## Scripts

### fermion_mass_primes_test11.py

**Multi-scale single-prime substrate test (minimal grammar)**

Purpose:
Tests whether fermion mass ratios are best described as a single odd prime, optionally multiplied by a power of 2. This isolates the hypothesis that ratios are either bare prime modes or substrate-dressed versions of them.

Physical motivation:
If fermion masses have the form ( m = 2^a \times (\text{mode}) ), then ratios take the form ( 2^{(a_1 - a_2)} \times p ).

* If ( a_1 = a_2 ), the substrate cancels → pure odd prime
* If ( a_1 \neq a_2 ), the ratio carries a residual power of 2

Grammar:
( 2^a \times p )

* ( a \in {0,1,2,3} )
* ( p ) = single odd prime
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

  * find best-fit ( 2^a \times p )
  * compute residual
* no enumeration over prime sets; direct nearest-fit search

Output:

* best-fit ( (a, p) ) per ratio and scale
* residuals
* frequency of ( a = 0 ) vs ( a > 0 )
* stability of prime assignment across scales

Key questions:

* Do leptonic ratios resolve as pure primes (( a = 0 ))?
* Do quark ratios require substrate dressing (( a > 0 ))?
* Is a single-prime model sufficient without composite products?

### fermion_mass_primes_test12.py

**Hierarchy-matched single-prime substrate null test (multi-scale comparison)**

Purpose:
Tests whether the proximity of fermion mass ratios to single-prime substrate forms (2^a \times p) (test11) is statistically unusual, by comparing against hierarchy-matched random ratios.

Scope:

* 6 within-family ratios only
  (mc/mu, mt/mc, ms/md, mb/ms, mmu/me, mtau/mmu)
* Real data evaluated at all 10 energy scales
  (PDG mixed + 9 Antusch scales)
* Random trials evaluated once per trial (no scale dependence)

Null model (same generator as test07):

* For each trial, generate six independent ratios:
  (r_i \sim \mathrm{LogUniform}(8.4, 1179))

Scoring (test11 metric):

* For each ratio, find best-fit candidate of the form:
  (2^a \times p)
  where:

  * (a \in {0,1,2,3})
  * (p) = single odd prime
* Compute residual:
  [
  \left|\frac{r_i - (2^a p)}{r_i}\right| \times 100
  ]
* Aggregate per trial: mean residual across 6 ratios

Procedure:

* For each real scale:

  * compute mean residual using test11 best fits
* For each of N trials (e.g. 1000 or 1M):

  * generate random six-ratio set
  * compute mean residual using same (2^a \times p) fitting
* Compare real vs random distribution

Output:

* Per scale:

  * real mean residual
  * percentile within random distribution
* Random distribution summary:

  * mean, median, percentiles, min, max
* Distinct prime usage (real vs random)

Key questions:

* Are real fermion ratios closer to (2^a \times p) than hierarchy-matched random ratios?
* Does the signal persist across all energy scales?
* Do real ratios exhibit a stronger preference for low substrate powers ((a=0,2)) than random data?

Interpretation:

* If real residuals lie consistently in the lower tail of the random distribution, this supports a structured single-prime substrate model
* If real and random overlap substantially, the single-prime model may not capture unique structure

Note:
This test directly extends test10 by replacing the nearest-prime metric with the single-prime substrate model from test11.

### fermion_mass_primes_test14.py

**Cross-family same-generation boundary test**

##### Purpose

Tests whether the test11 grammar,

[
R \approx 2^a \times p
]

describes same-generation cross-family mass ratios as well as it describes within-family ratios.

This is a boundary test, not a grammar extension.

##### Grammar

Identical to test11:

[
2^a \times p
]

where:

* (a \in {0,1,2,3})
* (p) = single odd prime
* no composite products
* no prime powers
* no denominators

The smallest grammar candidate is:

[
3
]

so ratios below 3 are flagged as below-floor.

##### Scope

Nine same-generation cross-family ratios:

[
d/u,\ c/s,\ t/b,\ u/e,\ c/\mu,\ t/\tau,\ d/e,\ \mu/s,\ b/\tau
]

Evaluated at 10 energy scales:

* PDG mixed
* 9 Antusch common scales, identical to test11

Neutrino sector excluded.

##### Orientation

For each cell, compute the heavier/lighter ratio:

[
R = \frac{\max(m_i,m_j)}{\min(m_i,m_j)}
]

Record the orientation.

If orientation swaps between scales, flag the pair as a unity-crossing.

##### Method

For each ratio at each scale:

* compute (R)
* flag whether (R < 3)
* find best-fit (2^a \times p)
* compute residual:

[
\left|\frac{R - 2^a p}{R}\right| \times 100
]

##### Output

Per ratio and scale:

* (R)
* orientation
* below-floor flag
* best ((a,p))
* residual

Per ratio:

* mean residual across scales
* number of below-floor cells
* prime stability: locked / walking
* substrate stability: whether (a) is stable or changes
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

##### Key questions

* Do cross-family residuals look like within-family residuals?
* Are cross-family ratios often below the grammar floor?
* Do any cross-family ratios show locked or walking behaviour like test11?
* Are apparent good fits clustered by ratio, or just isolated hits?
* Does any unity-crossing occur?

##### Interpretation

If cross-family residuals are worse, or many ratios fall below the grammar floor, the test supports the view that the (2^a p) grammar describes within-family excitation structure rather than cross-family parameter differences.

If cross-family ratios fit comparably well, especially with stable locked or walking behaviour, then the family-independence interpretation needs rethinking.

### fermion_mass_primes_test15.py

**Sector prime-set examination**

Purpose:
Examines the common-scale prime vocabularies produced by the test11 grammar, split by charged-fermion sector. No new fits are computed. The goal is to ask whether sector prime sets are structurally informative or mostly artefacts of magnitude range and nearest-candidate fitting.

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
  - up-sector walking behaviour

---

Key questions:

1. Are sector prime vocabularies genuinely disjoint beyond magnitude effects?
2. Are the vocabularies sparse or dense within their own prime windows?
3. Which primes are repeatedly occupied versus visited once?
4. Do sector prime sets show simple algebraic structure?
5. Are good residuals meaningful, or expected from local candidate density?
6. Do the sector signatures support different internal parameterisations of the same underlying soliton hardware?

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

### fermion_mass_primes_test16.py

**Multi-scale chain composition (accounting check).**

What this test is. An accounting and propagation check on whether existing grammar-fitted edges from test11 and test14 compose coherently along a chosen tree path. It is not independent evidence for modal equivalence; it confirms only that the chosen edge fits propagate correctly under multiplicative composition.

The chains tested at each common scale:

- Lepton chain: m_t / m_e = (t/τ) × (τ/μ) × (μ/e) — one cross-family bridge plus two within-family locks
- Down chain: m_t / m_d = (t/b) × (b/s) × (s/d) — one cross-family bridge plus two within-family edges
- Up chain: m_t / m_u = (t/c) × (c/u) — two within-family edges only

The lepton and down chains are the more informative ones because they pass through cross-family bridges that v1.0 found to fit individually. The up chain is mostly arithmetic consistency on within-family fits and is included for completeness.

This is a stronger hypothesis than v1.0. v1.0 found that the grammar fails as a general cross-family description, with three of nine cross-family ratios fitting individually at sub-1% (t/τ, t/b, c/s) and six failing systematically. test16 admits the two cross-family bridges (t/τ, t/b) that fit individually as candidate modal edges. The test is therefore a stress check on a stronger modal-equivalence reading than v1.0 commits to, not a confirmation of the v1.0 claim.

Why this is largely an arithmetic identity. Observed ratios are products of observed component ratios: m_t/m_e = (m_t/m_τ)(m_τ/m_μ)(m_μ/m_e) by definition. Predictions are products of fitted component ratios. The composed residual is therefore mathematically determined by the signed component residuals: |∏(1+δᵢ) − 1| with δᵢ the signed errors of each component. This identity must hold; the test is informative chiefly through what it reveals about cancellation versus reinforcement of component errors.

Signed residuals, not absolute. Component errors can cancel or reinforce under composition. The test reports signed component residuals δᵢ, the multiplicative product ∏(1+δᵢ), and the directly computed composed residual at each scale. Cancellation patterns — for example, a chain whose composed residual is small only because component over- and under-predictions happen to cancel — are themselves informative even when the chain looks like it "passes". A chain that composes well through systematic cancellation tells a different structural story than one that composes well through uniformly small component errors.

PDG excluded. The lepton chain breaks at PDG because m_μ/m_e fits 2 × 103 rather than the locked p = 211 (v1.0 §3.1, attributed to mixed-convention conventions). Including PDG selectively — using it where it works and dropping where it doesn't — would constitute post-hoc data selection. The test operates on the nine common-scale rows only, where the convention is uniform.

Diagnostic criterion. For each chain, the directly computed composed residual must equal the residual implied by the signed component errors,

Δ_chain = ∏ᵢ(1 + δᵢ) − 1,

to numerical tolerance. Any discrepancy indicates a bookkeeping, orientation, sign-convention, or implementation error. A successful diagnostic is not evidence for modal equivalence; it only confirms that the chosen edge fits propagate correctly along the tree.

Inputs. test11 fit data at nine common scales (within-family edges); test14 fit data at nine common scales (cross-family bridges t/τ, t/b).

Outputs. Per-scale, per-chain table showing component signed residuals, multiplicative product ∏(1 + δᵢ), predicted versus observed composed ratio, directly computed composed residual, and the implied-versus-direct discrepancy (which should be zero to numerical tolerance). Summary of cancellation versus reinforcement patterns by chain.

Note on terminology. "Charged fermions" in this work refers to the nine charged-fermion mass eigenstates (three charged leptons plus six quarks). Neutrinos are outside the scope of v1.0 and v1.1; see v1.0 §2.3.

### fermion_mass_primes_test17.py

**Bridge substitution null (local enumeration).**

What this test is. A local null check on whether the cross-family bridges (t/τ, t/b) used in the chain hypothesis carry structural information beyond their numerical magnitude. Where test16 confirmed bookkeeping, test17 asks whether the specific bridge candidate selected by the grammar is unusually well-aligned with the value the chain endpoint actually requires, compared to other grammar candidates of similar magnitude.

Procedure. For each common scale and each chain that passes through a cross-family bridge:

1. Hold the within-family components fixed at their v1.0 grammar fits.
2. Compute the required bridge value:

   B_required = (m_t / m_X)_observed / Π(within-family fitted factors)

   where X is the gen-1 charged fermion at the chain endpoint (e for lepton chain, d for down chain). This is the value the bridge would need to take for the chain to predict the observed endpoint exactly, given the fitted within-family edges.

3. Take the observed bridge ratio B_obs from running mass data at that scale (e.g., m_t/m_τ for the lepton chain).
4. Enumerate every grammar candidate 2^a × p with a ∈ {0,1,2,3} and p odd prime that falls in the magnitude window [B_obs/2, 2·B_obs]. This is a finite set, typically dozens of candidates.
5. For each candidate in the pool, substitute it into the chain in place of the actual bridge and compute the absolute residual against the observed endpoint:

   residual(c) = |Π(within-family fits) × c − (m_t/m_X)_observed| / (m_t/m_X)_observed

6. The actual bridge candidate B_fit = 2^a × p selected by test14 is also in this pool. Rank B_fit within the pool by residual; report its rank and the pool size.

Why required versus observed bridge. Random substitutes get scored against B_required, not against B_obs. The within-family fits aren't perfect; they introduce small offsets at each chain step. B_required is the bridge value that would compose with those imperfect within-family fits to land exactly on the observed endpoint. Comparing every candidate (including B_fit) to B_required measures who composes best in chain context — not who is closest to the observed bridge in isolation.

Two chains in scope. Lepton chain (bridge t/τ) and down chain (bridge t/b). The up chain has no cross-family bridge; nothing to substitute.

Magnitude window of factor two. Pre-committed before running. Wide enough to admit a non-trivial pool (typically dozens of candidates), narrow enough to control for gross magnitude mismatch. The pool size is reported per scale because a low percentile means different things if the pool has 12 versus 80 candidates.

Exact enumeration, not Monte Carlo. The candidate pool inside [B_obs/2, 2·B_obs] is a finite, deterministic set. Every candidate is evaluated. No sampling, no percentile-resolution issues, no random seeds to record.

Percentile definition. Rank 1 is the best (lowest residual) candidate in the pool. Percentile is defined as

percentile = 100 × (rank − 1) / (N − 1)

where N is the pool size. The best candidate sits at percentile 0; the worst at percentile 100.

Tie handling. If two or more candidates have identical residuals to floating-point precision, average rank is assigned to the tied set. Pool members and their tied ranks are reported alongside B_fit's rank.

Pool membership assertion. Before ranking, the script asserts that B_fit (the bridge candidate selected by test14) is in the enumerated pool. If B_fit is not in the window, the null definition has failed: a hard diagnostic prints and the test halts. This catches edge cases where the test14 fit lies just outside [B_obs/2, 2·B_obs].

PDG excluded. Same reason as test16 — the lepton chain's μ/e component breaks at PDG.

Pre-committed reading rules (not pass/fail).

The output is a rank distribution. Phrased as evidence strength, not binary verdict:

- Strong bridge evidence at a scale: percentile ≤ 5
- Modest bridge evidence at a scale: 5 < percentile ≤ 25
- Weak bridge evidence at a scale: 25 < percentile ≤ 50
- No bridge evidence at a scale: percentile > 50 — the actual bridge is no better than typical magnitude-comparable candidates.

A single weak scale doesn't collapse the interpretation. The chain hypothesis is supported in proportion to how often and how strongly the actual bridge ranks across the nine common scales. Aggregate readings are reported as "strong at N/9 scales, modest at M/9, weak at K/9, none at J/9".

Scope of the claim. This is a local bridge-density null. It tests whether B_fit is unusually well-aligned with B_required among nearby grammar candidates. A strong rank supports the bridge edge being more than magnitude-compatible — the specific candidate selected appears to do structural work. It does not establish modal equivalence on its own.

Caution on interpretation. Because B_fit was selected by test14 as the nearest grammar candidate to B_obs, and B_required is usually close to B_obs, a strong rank partly reflects nearest-candidate behaviour at the bridge step. The test is strongest at scales where B_required is noticeably shifted from B_obs by within-family fit errors and B_fit still ranks well — those cases show the bridge composing favourably in chain context, not just sitting close to the observed bridge ratio. Per-scale reporting includes both B_obs and B_required so this shift is visible to the reader.

Inputs. test11 fits at nine common scales (within-family components, used to compute B_required). test14 fitted bridge values (the actual B_fit being ranked). Running mass data (for B_obs and observed endpoints).

Outputs. Per-scale, per-chain table showing: pool size, B_obs, B_required, B_fit, rank of B_fit in pool, percentile, residuals of B_fit and best-pool-candidate. Tied ranks reported where they occur. Aggregate summary across nine scales.

### fermion_mass_primes_test18.py

**Hierarchy-matched 9 charged-mass null.**

What this test is. A global null check on whether the joint structure of v1.0's grammar fits — within-family sub-percent fits, gen-3 cross-family bridges, and their compositional consistency — occurs more often in real charged-fermion data than in random hierarchies of comparable shape.

This is the strongest joint-structure test in the v1.1 programme. Where test16 confirmed bookkeeping and test17 tested individual bridge candidates locally, test18 asks whether the combination of structural features observed in the real data is statistically distinguishable from the same features in random hierarchies.

Scope: charged fermions only. The null is generated over the nine charged-fermion mass eigenstates (three charged leptons plus six quarks). Neutrinos are excluded entirely from the null — including them as hierarchy constraints would import out-of-domain assumptions about absolute neutrino masses that v1.0 §2.3 explicitly places outside the framework.

Per-scale evaluation. The test is evaluated independently at each common scale. It does not test cross-scale stability of locks; it tests whether the real per-scale co-occurrence pattern is unusual under the hierarchy-matched charged-mass null.

Random-system generator (pre-committed). Random charged-mass systems are generated as follows, independently per scale:

1. Use the real charged-fermion rank order at the given scale (e.g., at M_Z: m_t > m_b > m_c > m_τ > m_s > m_μ > m_d > m_u > m_e). Some inter-sector orderings change with scale (e.g., b/τ crosses around 10⁷ GeV in the data); the rank order used here is the empirical one at each scale.
2. Draw the total log span S = log(m_max / m_min) from a uniform distribution on [S_real − log 2, S_real + log 2], where S_real is the real-data log span at that scale.
3. Draw eight ordered log gaps from a uniform Dirichlet distribution on the 8-simplex, scaled to sum to S.
4. Construct nine ordered log masses by cumulative sum from m_min (held fixed at the real-data value to set the absolute scale, though only ratios enter the scoring).
5. Assign mass labels (u, d, s, c, b, t, e, μ, τ) according to the real rank order from step 1.

This gives plausibly-ordered random hierarchies without baking in any prime-grammar structure.

Scoring features (pre-committed).

For each system (real or random) at each scale, four features are scored:

1. Within-family sub-percent fit count: number of the six within-family adjacent-generation ratios fitting 2^a × p at residual ≤ 1%.
2. Selected cross-family fit count (low-residual cross-family count): number of the three selected cross-family ratios (t/τ, t/b, c/s) fitting 2^a × p at residual ≤ 3%. These are the three cross-family ratios that v1.0/test14 found to fit individually; c/s is included for completeness even though it does not appear in test16's chain construction.
3. Lepton chain absolute residual: absolute value of the composed residual of m_t/m_e under test16's lepton chain.
4. Down chain absolute residual: absolute value of the composed residual of m_t/m_d under test16's down chain.

Joint score. A system "matches or exceeds" the real data at a given scale if it has within-family fit count ≥ real, cross-family fit count ≥ real, and absolute lepton and down chain residuals ≤ real absolute residuals.

Note on feature independence. The four features are not statistically independent. Chain residuals depend on the within-family fits and bridge fits that compose them. The joint score therefore tests co-occurrence of the structural features, not their independent multiplication into a combined p-value. The reported fraction is "fraction of random systems matching or exceeding the real-data joint score", interpreted as co-occurrence frequency.

Number of trials. 1,000,000 random systems per scale. Matches test12's trial count and gives 1e-6 fraction resolution at the tail.

Implementation note. The grammar fit search runs per ratio per random system. Naive linear scan over the candidate set would be unnecessarily slow; the implementation pre-sorts candidates and uses bisect for nearest-candidate lookup, reducing per-fit cost from O(N) to O(log N). The candidate set is sized to cover the full range of ratios producible by random systems (up to exp(S_real + log 2) at the largest-span scale), not just the real-data ratio range, so that random systems are not artificially punished by candidate truncation. Below-floor protection: ratios below the grammar floor (R = 3) are not eligible for fit counting, and chain residuals with any below-floor component are reported as infinite (cannot match real data).

Zero-count reporting. If no random sample at a given scale matches or exceeds the real data, the reported fraction is < 1/N (i.e., < 1e-6 at one million trials), not 0%. This preserves the resolution limit of the test rather than asserting zero.

PDG excluded. Same reason as test16/17 — convention-mixing breaks the lepton chain at PDG.

Pre-committed reading rules.

Per-scale evidence bands:

- Strong: fewer than 0.1% of random hierarchies match or exceed the real data at that scale.
- Modest: fewer than 1%.
- Weak: fewer than 5%.
- None: 5% or more — the real-data joint structure is comparable to typical random plausible hierarchies.

Aggregate reading uses the N/9 framing established in test17: "strong at N/9 scales, modest at M/9, weak at K/9, none at J/9". Not combined into a single p-value, since random samples are generated independently per scale.

Why this is the strongest joint test.

The within-family lock claim alone admits a null (v1.0's test12 ran one). The bridge-fit claim alone is testable (test17). The chain-composition claim alone is testable (test16, weakly). Test18 asks whether the structural features observed in the real data co-occur more often than in random plausible hierarchies. A real underlying structure should produce sub-percent fits, bridges, and clean chain compositions together; random hierarchies should produce them independently with low joint frequency.

If the real data's joint pattern is reproduced by random hierarchies more than 5% of the time at multiple scales, the v1.1 modal-equivalence reading loses its empirical foundation. If it's reproduced rarely, that supports the structural reading.

Caveat on the null definition. The Dirichlet log-gap generator is one specific choice for "random plausible hierarchy". Alternatives (log-uniform on individual masses, geometric distributions on gaps, etc.) would produce different random pools and potentially different outcomes. The Dirichlet choice is symmetric, simple, and uninformative — it places no preference on any particular gap structure beyond ordering. Alternative null definitions can be reported as robustness checks if the headline result is borderline.

Inputs. Real running mass data at nine common scales (Antusch et al. tabulation). Within-family fits computed on each random system using the same algorithm as test11. Cross-family fits computed on each random system using the same algorithm as test14. Chain residuals computed using test16's chain construction.

Outputs. Per-scale table showing real-data scores on all four features, random-system match-or-exceed fractions (with < 1/N reporting where appropriate), and per-scale evidence band. Aggregate summary across nine common scales using N/9 framing. Histograms of random-system scores per feature for inspection.

### fermion_mass_primes_test19.py

**Held-out cycle closure.**

What this test is. A graph-closure test on the modal-equivalence reading of the chain hypothesis. Test16 confirmed bookkeeping along the tree. Test17 tested individual bridges locally. Test18 tested joint co-occurrence of structural features. Test19 asks the harder question: do held-out cross-family edges — those NOT used to construct the chain tree — close consistently when predicted by the tree's composition?

If the modal-equivalence reading holds (charged fermions as mode positions on a single underlying object), the tree predicts a value for every cross-family ratio, including ones that weren't used as tree bridges. Comparing tree-predicted values against observed values at every common scale tests whether the modal structure is internally consistent beyond the edges it was built from.

The tree's bridges (used). The chain construction in test16 uses two cross-family edges as tree bridges:

- t/τ (lepton chain)
- t/b (down chain)

The up chain uses two within-family edges (t/c, c/u) and no cross-family bridge.

Held-out cross-family edges (predicted via composition). Of the nine cross-family ratios in v1.0/test14, two are used as bridges (t/τ, t/b). The remaining seven are held out. For each, the tree predicts a value via composition through the within-family edges and the two bridges:

| Held-out ratio | Tree prediction |
|---|---|
| c/s | (t/b × b/s) / (t/c) |
| d/u | (t/c × c/u) / (t/b × b/s × s/d) |
| u/e | (t/τ × τ/μ × μ/e) / (t/c × c/u) |
| c/μ | (t/τ × τ/μ) / (t/c) |
| d/e | (t/τ × τ/μ × μ/e) / (t/b × b/s × s/d) |
| μ/s | (t/b × b/s) / (t/τ × τ/μ) |
| b/τ | (t/τ) / (t/b) |

Seven held-out edges × nine common scales = 63 cells of compositional prediction versus observation.

Procedure. For each cell:

1. Compute the tree-predicted value for the held-out ratio using the composition above, with all factors taken from their v1.0 grammar fits (test11 within-family edges, test14 cross-family bridges).
2. Compute the observed (directed) value of the held-out ratio from running mass data at that scale. Tree composition predicts directed labelled ratios (e.g., b/τ as m_b/m_τ, which crosses unity across scales), so compare directed-to-directed.
3. Compute the compositional residual: (predicted − observed) / observed (signed); report absolute value as percentage.
4. Compute the implied compositional residual from exact signed error propagation:

   Δ_implied = ∏ᵢ (1 + δᵢ)^{sᵢ} − 1

   with sᵢ = +1 for numerator factors and sᵢ = −1 for denominator factors, and δᵢ the signed residuals of the v1.0 grammar fits used in the composition.
5. Compute the identity check discrepancy:

   identity_discrepancy = signed_compositional_residual − Δ_implied

   This must be zero to numerical tolerance by construction. Any non-zero value indicates a bookkeeping, orientation, sign-convention, or implementation error — not physical inconsistency. (Reported as a discrepancy rather than a quotient, since |Δ_implied| can be near zero for cells with cancelling component errors.)
6. Compute the directly-fitted grammar residual for the same ratio (test14 result), if the ratio is above the grammar floor. For the test14 fit, report the heavier/lighter orientation and the fitted (a, p) value alongside the directed-comparison residual.

This gives a per-cell record: predicted, observed, compositional residual, implied compositional residual, identity check discrepancy, direct-fit residual (or "below floor" marker).

The 2×2 interpretive matrix per cell.

A cell's "direct grammar fit small" means direct-fit residual ≤ 3%, matching the test14/test18 cross-family convention.

| | direct grammar fit small (≤ 3%) | direct grammar fit large or below floor |
|---|---|---|
| compositional residual small | tree closes AND ratio has direct compact form | tree closes for a ratio that doesn't fit the direct grammar at this scale; the fitted tree predicts the directed ratio better than the local direct candidate |
| compositional residual large | challenges the chosen tree — ratio fits its own grammar but doesn't compose | tree open for a ratio with no compact direct form either |

The "direct large, compositional small" quadrant says the fitted tree predicts the directed ratio better than the local grammar candidate; it does not by itself establish that modal structure reaches beyond the grammar, only that tree composition is closer to truth than direct nearest-candidate fitting at that cell. The "direct small, compositional large" quadrant is the genuinely challenging one for tree topology.

Cycle interpretation. Predicting a held-out edge from the tree creates an implicit cycle in the equivalence-class graph. For example, the c/s prediction implies the loop t → c → s → b → t. If the fitted tree is coherent, the loop product closes to the level implied by the fitted edge errors. Cycle closure tests this for every held-out edge.

Algebraic note on independence. Held-out compositional predictions are not independent new data. The held-out ratios are algebraically dependent on the same charged-fermion masses that produced the tree's edges. Test19 is therefore an error-propagation and graph-closure test on the fitted tree, not a forward prediction in the sense of unseen data. The information it provides is whether the fitted edge set propagates consistently to held-out ratios.

Below-floor handling. Some held-out ratios (notably d/u, u/e, μ/s, d/e) fall below the grammar floor (R = 3) at most or all common scales. For such cells:

- The compositional prediction is computed normally — composition is not bounded by the grammar floor since it returns a real-valued ratio.
- The "direct grammar fit" comparison is reported as "below floor" rather than as a residual.
- The 2×2 matrix collapses to two cells (compositional residual small or large) for below-floor ratios, since direct fit isn't applicable.

A held-out ratio below the grammar floor whose composition still matches observation is informative: the fitted tree closes through a region where direct grammar fit is undefined.

PDG excluded. Same reason as test16/17/18 — convention-mixing breaks the lepton chain at PDG.

Pre-committed reading rules (per cell).

For each held-out cell, the result is classified by the size of the compositional residual:

- Strong: compositional residual ≤ 5%
- Modest: 5% < residual ≤ 15%
- Weak: 15% < residual ≤ 30%
- None: residual > 30%

The identity check discrepancy is reported separately as a bookkeeping diagnostic. It must be zero to numerical tolerance by construction; any non-zero discrepancy indicates an implementation error, not physical inconsistency.

Aggregate readings:

- Cell-level: N/63 cells in each band.
- Per-edge: N/9 cells per held-out ratio across the nine scales.

Pass criterion for the chain hypothesis. The chain hypothesis is supported if held-out compositional residuals are systematically in the strong or modest band. The hypothesis is challenged if compositional residuals are systematically weak or none — particularly if the failure pattern is structured (e.g., always at certain scales or always for certain ratio types). The identity check discrepancy is reported as an implementation diagnostic only, not as evidence for or against the hypothesis.

Specifically informative cases.

1. Held-out ratios below the grammar floor (d/u, u/e, μ/s, d/e) where the tree predicts the correct directed value with small compositional residual: shows the fitted tree closes through ratios where direct grammar fit is undefined.
2. Held-out ratios above the floor where compositional residual is significantly smaller than the direct-fit residual: shows tree composition is a better predictor than direct nearest-candidate fitting at that cell. Does not by itself establish that modal structure reaches beyond the grammar.
3. Held-out ratios above the floor where compositional residual matches direct-fit residual: tree closes consistently at the expected level of accuracy.
4. Held-out ratios where compositional residual is large AND direct-fit residual is small (≤ 3%): challenges the tree topology. The ratio fits its own grammar but does not compose with the chosen tree.

Inputs. test11 fit data at nine common scales (all within-family edges in the tree, with signed residuals for propagation). test14 fit data at nine common scales (the two cross-family bridges used in the tree, the seven held-out cross-family edges for direct-fit comparison, and signed residuals for propagation). Running mass data (for observed directed held-out ratios at each scale).

Outputs. Per-scale, per-held-out-edge table showing: directed predicted value, directed observed value, signed compositional residual, implied compositional residual (Δ_implied), identity check discrepancy, direct grammar fit (heavier/lighter orientation, fitted (a, p), residual) or "below floor" marker, and per-cell classification band. Aggregate summary across 63 cells. Per-edge summary across nine scales. Identification of any cells in the challenging quadrant (compositional large, direct small ≤ 3%).

Note on terminology. "Charged fermions" refers to the nine charged-fermion mass eigenstates (three charged leptons plus six quarks). Neutrinos are outside the scope of v1.0 and v1.1.

### fermion_mass_primes_test21.py

**CKM mixing null — random prime assignment under fixed structural template.**

What this test is. A null check on whether the specific chain-derived prime values and their slot placement produce better CKM magnitude predictions than random assignments under the same structural template. This test does not ask whether the full CKM model arises by chance. It asks the conditional question: given the discovered mode-slot topology, do the real prime values outperform random values?

##### Purpose

The impedance network model computes all nine CKM magnitudes from the mode decomposition of the nine charged fermions:

|V_ij|² ∝ exp(−Z_total), normalised per row.

Z_total = Z_bal + Z_ch + Z_R + Z_a + Z_deg + Z_tp

where:

1. Z_bal = α_s · |Σp_shed − Σp_gain| — energy balance (surplus path through substrate).

2. Z_ch = −ln(K), where K = Σ exp(−α_s · |p−q| · max(p,q)/p_max) / √N — channel conductance with frequency-dependent ring impedance. N = n_shed × n_gain. For pure shedding (no gain): Z_ch = 0.

3a. Z_R (engaging, a_tgt > a_src) = (π/2) · |Δa| · J, where J = Σln(p_shared) / (Σln(p_shed) + Σln(p_gain)) — substrate change with shared-mode rewiring load.

3b. Z_R (releasing, a_tgt < a_src) = (π/2) · |Δa| · (1 − α_s · Σp_rewired/P + α_s · Σp_gain · |Δa|/P) — cooperative relaxation with gain detuning.

4. Z_a = (α_s/P) · (Σp_shed · a_src + Σp_gain · a_tgt) — substrate throughput (pump).

5. Z_deg = (α_s/P) · |Δa| · (Σp_shed + Σp_gain) · J — channel degradation during substrate change, where J is the same shared-mode rewiring load used in the engaging term. Applies to Δa ≠ 0 transitions with both shed and gain.

6. Z_tp = (α_s / (n_src · (n_src − 1))) · Σp_shed — ring traversal for pure shedding at Δa = 0 only. n_src = number of modes on the source fermion.

The test asks: is the fit quality achieved by the real chain-derived primes {5, 17, 53, 59, 97, 127, 211, 271} in their specific slot assignments statistically unusual compared to (a) the same primes in shuffled slots, and (b) random primes in the same slots?

##### Structural template (preserved across all trials)

The following structural features are held fixed:

1. Three sectors: up-quark, down-quark, charged lepton.
2. Three generations per sector, with substrate values: a = 0 for u and d; a = 2 for all others.
3. Eight structural prime slots (A–H), assigned to sectors as:
   - Slot A: Up identity, Lepton identity, Down generation
   - Slot B: Up identity, Down identity, Lepton generation
   - Slot C: Up identity, Lepton identity, Down generation
   - Slot D: Up identity, Lepton identity, bridge (absent from Down)
   - Slot E: Up identity, Down identity, bridge (absent from Lepton)
   - Slot F: Down identity, Lepton identity, Up generation
   - Slot G: Up identity, Down identity, Lepton generation
   - Slot H: Down identity, Lepton identity, Up generation
4. Real slot assignment: A=5, B=17, C=53, D=59, E=97, F=127, G=211, H=271.
5. The scoring formula, fixed constants (α_s = 0.1180, π/2 = 1.5708), derived constants (P = sum of all 8 assigned primes, p_max = largest assigned prime), n_src per particle determined by slot count (u=6, c=7, t=8, d=5, s=6, b=7), and per-row normalisation.

##### Particle mode reconstruction from slots

Each particle's odd-prime set is determined by its sector and generation. The explicit reconstruction:

```
u = {A, B, C, D, E, G}           (up gen-1: identity only)
c = {A, B, C, D, E, F, G}        (up gen-2: identity + F)
t = {A, B, C, D, E, F, G, H}     (up gen-3: identity + F, H)

d = {B, E, F, G, H}              (down gen-1: identity only)
s = {A, B, E, F, G, H}           (down gen-2: identity + A)
b = {A, B, C, E, F, G, H}        (down gen-3: identity + A, C)

e = {A, C, D, F, H}              (lepton gen-1: identity only)
μ = {A, C, D, F, G, H}           (lepton gen-2: identity + G)
τ = {A, B, C, D, F, G, H}        (lepton gen-3: identity + B, G)
```

This table is invariant across all trials. Only the values assigned to slots A–H change.

##### Null populations

Three null populations are tested separately:

1. **Permutation null.** The same eight real primes {5, 17, 53, 59, 97, 127, 211, 271}, randomly permuted across slots A–H. This isolates whether the specific slot placement of the real primes matters. The real assignment is included as one of the 8! permutations; the minimum nonzero match fraction is 1/40320 ≈ 0.00248%, which falls in the "strong" band.

2. **Range-matched null.** 8 distinct odd primes drawn uniformly from [3, 271] (the real prime range), assigned to slots A–H in the order drawn (not sorted). This tests whether any 8 primes in the right magnitude range would work.

3. **Extended-range null.** 8 distinct odd primes drawn uniformly from [3, 997], assigned to slots A–H in the order drawn (not sorted). This tests whether the result depends on the specific magnitude range.

For populations 2 and 3, primes are assigned in draw order to prevent magnitude ordering from being baked into the slot structure.

##### Scoring

Primary score. For each trial, the nine CKM probabilities |V_ij|² are computed from the mode decomposition. The score is:

S = Σ_{i,j} [ ln( |V_ij|²_pred / |V_ij|²_obs ) ]²

summed over all 9 elements. Lower S = better fit. This is equivalent to 4 × Σ [ln(|V_ij|_pred / |V_ij|_obs)]² — a factor of 4 relative to a magnitude-log score. Rankings and null fractions are unaffected; This is the canonical score used for test21 reporting. If magnitude-log or absolute-log summaries are printed, they are labelled separately and are not used for the headline null fraction.

Observed values |V_ij|²_obs are computed from Antusch et al. (arXiv:2510.01312v2) at M_Z: sin θ₁₂ = 0.2251, sin θ₂₃ = 0.04193, sin θ₁₃ = 0.00370, δ = 1.139 rad.

Secondary score. Maximum absolute log error across all elements:

S_max = max_{i,j} | ln( |V_ij|_pred / |V_ij|_obs ) |

This treats small and large CKM elements consistently.

##### Hierarchy diagnostics

1. **Global inversion count** (Kendall tau distance): number of pairwise rank inversions between predicted and observed |V_ij| orderings across all 36 pairs. Computed from actual numeric values with no tolerance. The real model has 1 inversion (V_us versus V_cd, where the observed values differ by 0.06%). The distribution of inversion counts across random trials is reported.

2. **Row-wise ordering**: check whether the correct row hierarchy holds for all three rows:

```
|V_ud| > |V_us| > |V_ub|
|V_cs| > |V_cd| > |V_cb|
|V_tb| > |V_ts| > |V_td|
```

Report fraction of trials achieving correct ordering on all three rows simultaneously. Note: the observed |V_us| and |V_cd| are very close (0.2251 vs 0.2250), so global ordering is sensitive to small perturbations in the u→s / c→d balance.

##### Column unitarity diagnostic

Column unitarity is scored as a scalar:

U_col = max_j | Σ_i |V_ij|²_pred − 1 |

i.e. the worst column-sum deviation from unity across the three columns. The distribution of U_col across random trials is compared to the real model's value.

##### Joint diagnostic (secondary, not headline)

Report as a diagnostic (not headline reading):

fraction of trials with S ≤ S_real AND U_col ≤ U_col_real AND inversion_count ≤ inversion_count_real

This asks how often random assignments simultaneously match the model's fit quality, column balance, and hierarchy ordering. Reported for context; the headline reading rests on S alone.

##### Number of trials

- Permutation null: all 8! = 40,320 permutations (exact enumeration).
- Range-matched null: 1,000,000 random draws.
- Extended-range null: 1,000,000 random draws.

##### Pre-committed reading rules

Per null population, using the primary score S:

- **Strong:** fewer than 0.1% of assignments match or exceed the real score.
- **Modest:** fewer than 1%.
- **Weak:** fewer than 5%.
- **None:** 5% or more.

The match-or-exceed fraction is: f = (number of trials with S ≤ S_real) / N_trials

If no trial matches, reported as < 1/N.

##### Expected outcome structure

A positive result (strong or modest on all three populations) means: the CKM result depends on both the specific chain-derived prime values AND their slot placement, not merely on the structural template, row normalisation, or the exponential form of the scoring function.

A result that is strong on the range-matched/extended nulls but weak on the permutation null would mean: the real prime values are special, but their slot assignment is not tightly constrained.

A result that is strong on the permutation null but weak on the range-matched null would mean: slot assignment matters but the specific prime values do not — any primes of similar magnitude would work.

##### Inputs

Mode decomposition structural template (particle-slot table above). CKM observed values from Antusch et al. at M_Z. Fixed constants: α_s(M_Z) = 0.1180, π/2 = 1.5708. Derived per trial: P = sum of all 8 assigned primes (real value: 840), p_max = largest assigned prime (real value: 271). n_src determined by slot count per particle: u=6, c=7, t=8, d=5, s=6, b=7.

##### Outputs

Per null population: real-data score (S and S_max), random-score distribution (mean, median, 1st/5th/25th percentiles), match-or-exceed fraction with evidence band, inversion count distribution, row-ordering fraction, column unitarity distribution. Joint diagnostic fraction. Best-performing trial comparison table showing all 9 predicted |V_ij| alongside the real model and observed values.

##### Implementation cautions

1. **Observed CKM construction.** The observed |V_ij| matrix must be computed from the full standard PDG parameterisation using sin θ₁₂, sin θ₂₃, sin θ₁₃, and δ — not approximated from the three sines alone. The magnitudes of V_td and V_ts depend on the CP phase δ through complex interference terms. The Antusch values (sin θ₁₂ = 0.2251, sin θ₂₃ = 0.04193, sin θ₁₃ = 0.00370, δ = 1.139 rad) are M_Z running values.

2. **Pure-shedding convention.** For pure-shedding transitions (empty gain set): Z_ch = 0 (no channel conductance — no gain modes to couple to); Z_deg = 0 (no channel degradation without channels); Z_tp applies only at Δa = 0 (ring traversal for pure shedding without substrate change); Z_R term 3b applies when Δa ≠ 0 (releasing). Under the fixed structural template, the t row is always pure-shedding (gain set empty for t→d, t→s, t→b) and no transition is ever pure-gaining, but the code must handle both cases for robustness.

3. **Permutation ties.** If another permutation produces the same S to floating-point precision, count as match-or-exceed. Use `S_trial <= S_real + eps` with eps at machine precision. Report the real assignment's rank among all 40,320 permutations explicitly (e.g., "rank 1 of 40,320").

4. **Uniform over primes, not integers.** For range-matched and extended nulls, precompute the complete list of odd primes in each interval ([3, 271] and [3, 997]) using a sieve. Sample 8 distinct primes uniformly from these lists. Do not sample from odd integers with primality rejection.

5. **Log-space computation.** For each row of the predicted |V|² matrix: compute log-weights for each target, subtract the row maximum, exponentiate, then row-normalise. This avoids underflow when extended-range primes produce large exponent arguments.

6. **Real rank in permutation null.** Report the real assignment's rank explicitly (e.g., "real assignment is rank 1 of 40,320"). This is the single clearest result for readers.

7. **P and p_max per trial.** P = sum of all 8 assigned primes and p_max = max of all 8 assigned primes must be recomputed for each trial. They are not fixed constants.

8. **Engaging vs releasing asymmetry.** The Z_R term has two branches depending on whether a_tgt > a_src (engaging, term 3a) or a_tgt < a_src (releasing, term 3b). When Δa = 0, Z_R = 0. The code must evaluate the correct branch per transition.

9. **Shared and rewired modes.** J (the shared-mode rewiring load in terms 3a and 5) and the rewired modes (in term 3b) are computed from the intersection and symmetric difference of the source and target prime sets. These depend on the specific mode decomposition, which changes with each trial's prime assignment. The shared set = source ∩ target; shed = source − target; gain = target − source; rewired = shared (they change substrate coupling).

##### Parameters

- N_TRIALS_RANDOM = 1,000,000
- SEED = 42
- PRIME_RANGE_MATCHED = [3, 271]
- PRIME_RANGE_EXTENDED = [3, 997]
- Permutation null: exact enumeration (40,320 permutations)
