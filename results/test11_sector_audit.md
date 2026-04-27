# test11 sector audit

Re-aggregation of test11 results (within-family ratios, 2^a × p grammar, 10 scales) split by sector. Source: `logs/test11_20260425_121126.log`. No new fits computed; this is descriptive only.

PDG mixed is reported separately from the 9 common scales (M_Z through 1e16 GeV) because PDG mixes renormalization conventions.

## Up sector

Ratios: `mc/mu`, `mt/mc`

### All 10 scales

| | mc/mu | mt/mc |
|---|---|---|
| PDG | 587 (a=0) | 2³×17 (a=3) |
| M_Z | 2²×127 | 271 |
| 1 TeV | 2²×127 | 277 |
| 3 TeV | 503 | 2×139 |
| 10 TeV | 503 | 281 |
| 100 TeV | 2²×127 | 2²×71 |
| 1e7 GeV | 2²×127 | 2²×73 |
| 1e9 GeV | 2²×127 | 2³×37 |
| 1e12 GeV | 2²×127 | 2×151 |
| 1e16 GeV | 503 | 307 |

- Primes used (all scales): {17, 37, 71, 73, 127, 139, 151, 271, 277, 281, 307, 503, 587} — 13 distinct
- Substrate powers (count over 20 cells): a=0: 8, a=1: 2, a=2: 8, a=3: 2
- Bare-prime fraction (a=0): 8/20 = 40%
- Locked ratios: none
- Walking ratios: both
- Mean residual: 0.308%

### Common scale only (9 scales, PDG excluded)

- Primes used: {37, 71, 73, 127, 139, 151, 271, 277, 281, 307, 503} — 11 distinct
- Substrate powers (count over 18 cells): a=0: 7, a=1: 2, a=2: 8, a=3: 1
- Bare-prime fraction: 7/18 = 39%
- Mean residual: 0.302%

PDG-specific primes: 587 (mc/mu), 17 (mt/mc — also appears in lepton sector across all scales)

## Down sector

Ratios: `ms/md`, `mb/ms`

### All 10 scales

| | ms/md | mb/ms |
|---|---|---|
| PDG | 2²×5 | 2²×11 |
| M_Z | 2²×5 | 53 |
| 1 TeV | 2²×5 | 2²×13 |
| 3 TeV | 2²×5 | 2²×13 |
| 10 TeV | 2²×5 | 2²×13 |
| 100 TeV | 2²×5 | 2²×13 |
| 1e7 GeV | 2²×5 | 2²×13 |
| 1e9 GeV | 2²×5 | 47 |
| 1e12 GeV | 2²×5 | 47 |
| 1e16 GeV | 2²×5 | 47 |

- Primes used (all scales): {5, 11, 13, 47, 53} — 5 distinct
- Substrate powers (count over 20 cells): a=0: 4, a=2: 16
- Bare-prime fraction: 4/20 = 20%
- Locked ratios: ms/md (p=5, a=2, all 10 scales)
- Walking ratios: mb/ms ({11, 13, 47, 53})
- Mean residual: 1.171%

### Common scale only (9 scales, PDG excluded)

- Primes used: {5, 13, 47, 53} — 4 distinct
- Substrate powers (count over 18 cells): a=0: 4, a=2: 14
- Bare-prime fraction: 4/18 = 22%
- Mean residual: 1.180%

PDG-specific primes: 11 (mb/ms)

## Charged lepton sector

Ratios: `mmu/me`, `mtau/mmu`

### All 10 scales

| | mmu/me | mtau/mmu |
|---|---|---|
| PDG | 2×103 | 17 |
| M_Z | 211 | 17 |
| 1 TeV | 211 | 17 |
| 3 TeV | 211 | 17 |
| 10 TeV | 211 | 17 |
| 100 TeV | 211 | 17 |
| 1e7 GeV | 211 | 17 |
| 1e9 GeV | 211 | 17 |
| 1e12 GeV | 211 | 17 |
| 1e16 GeV | 211 | 17 |

- Primes used (all scales): {17, 103, 211} — 3 distinct
- Substrate powers (count over 20 cells): a=0: 19, a=1: 1
- Bare-prime fraction: 19/20 = 95%
- Locked ratios: mtau/mmu (p=17, all 10 scales)
- Effectively locked: mmu/me (p=211, 9/10 scales; PDG outlier at p=103)
- Mean residual: 0.179%

### Common scale only (9 scales, PDG excluded)

- Primes used: {17, 211} — 2 distinct
- Substrate powers (count over 18 cells): a=0: 18 (100%)
- Bare-prime fraction: 18/18 = 100%
- Mean residual: 0.118%

PDG-specific primes: 103 (mmu/me)

## Cross-sector observations

Across all 10 scales, prime-sector overlap occurs in exactly one case: **p=17 appears in both the lepton sector (mtau/mmu, locked across all 10 scales) and the up sector (mt/mc at PDG only, as 2³×17)**. Every other prime is unique to one sector.

At common scales (PDG excluded), there is no prime-sector overlap at all. Every prime is sector-specific.

Substrate behaviour also separates by sector at common scales:

| | a=0 | a=1 | a=2 | a=3 |
|---|---|---|---|---|
| Up (18 cells) | 39% | 11% | 44% | 6% |
| Down (18 cells) | 22% | 0% | 78% | 0% |
| Lepton (18 cells) | 100% | 0% | 0% | 0% |

Lepton uses a=0 exclusively at common scales. Down uses a=2 dominantly. Up is the most variable, spanning all four substrate values.

Mean residual by sector (common scale): lepton 0.118%, up 0.302%, down 1.180%. Lepton tightest, down loosest. The down-sector residual is dominated by mb/ms (mean 1.65%); ms/md alone is 0.71%.

## Classification

Per the audit rules: candidate physical modes are stable across scales OR repeatedly sector-specific; walking primes are nodes visited, not properties.

Applied conservatively (a prime is stable only if locked, not just sector-specific while walking):

### Candidate physical modes

| Prime | Sector | Ratio | Form | Stability |
|---|---|---|---|---|
| 5 | down | ms/md | 2²×5 | locked all 10 scales |
| 17 | lepton | mtau/mmu | 17 | locked all 10 scales (also 2³×17 in up at PDG) |
| 211 | lepton | mmu/me | 211 | 9/10 scales (PDG outlier at p=103) |

The substrate factor 2² is part of the ms/md mode signature — the lock holds on the composite (a=2, p=5) value, not on p=5 alone in any other context.

### Nodes visited (walking)

- Up sector, mc/mu: {127, 503, 587}
- Up sector, mt/mc: {17, 37, 71, 73, 139, 151, 271, 277, 281, 307}
- Down sector, mb/ms: {11, 13, 47, 53}
- Lepton sector, mmu/me: {103} (PDG outlier only; 211 lock at common scales)

### Notes on PDG behaviour

PDG mixed introduces three primes that do not appear at any other scale: 587 (up), 11 (down), 103 (lepton). It also produces the only cross-sector appearance of any prime (17 in mt/mc at PDG). All four PDG-specific primes vanish under common-scale convention.

The mtau/mmu lock at p=17 holds even at PDG, so the prime-17 cross-sector appearance is a PDG-only artefact in the up sector, not a contamination of the lepton lock.
