"""
fermion_mass_primes_test12.py

Hierarchy-matched single-prime substrate null test (multi-scale comparison).

Compares real fermion mass ratio proximity to 2^a × p forms against
hierarchy-matched random ratios. Same null generator as test07,
same metric as test11.

See scripts/README.md for full documentation.
"""

import time
import random
import math
from sympy import nextprime

RATIO_NAMES = ["mc/mu", "mt/mc", "ms/md", "mb/ms", "mmu/me", "mtau/mmu"]
SUBSTRATE_POWERS = [0, 1, 2, 3]

# ── Sampling range (same as test07) ──
RATIO_MIN = 8.4
RATIO_MAX = 1179.0
LOG_RATIO_MIN = math.log(RATIO_MIN)
LOG_RATIO_MAX = math.log(RATIO_MAX)

N_TRIALS = 1000000
SEED = 42

# ── PDG mixed-convention masses (MeV) ──
PDG_MASSES = {
    "u": 2.16, "c": 1273.0, "t": 172560.0,
    "d": 4.70, "s": 93.5, "b": 4183.0,
    "e": 0.51099895, "mu": 105.6583755, "tau": 1776.93,
}

# ── Antusch et al. Table 2 (2024 PDG, SM, MS-bar) ──
ANTUSCH_SCALES = {
    "M_Z": {
        "yu": 7.04e-6, "yc": 3.56e-3, "yt": 0.967,
        "yd": 1.54e-5, "ys": 3.06e-4, "yb": 1.630e-2,
        "ye": 2.77713e-6, "ymu": 5.85042e-4, "ytau": 0.99378e-2,
    },
    "1_TeV": {
        "yu": 6.15e-6, "yc": 3.11e-3, "yt": 0.8616,
        "yd": 1.35e-5, "ys": 2.68e-4, "yb": 1.401e-2,
        "ye": 2.8227e-6, "ymu": 5.9465e-4, "ytau": 1.0101e-2,
    },
    "3_TeV": {
        "yu": 5.84e-6, "yc": 2.95e-3, "yt": 0.8242,
        "yd": 1.28e-5, "ys": 2.54e-4, "yb": 1.321e-2,
        "ye": 2.8402e-6, "ymu": 5.9833e-4, "ytau": 1.0163e-2,
    },
    "10_TeV": {
        "yu": 5.54e-6, "yc": 2.80e-3, "yt": 0.7880,
        "yd": 1.22e-5, "ys": 2.42e-4, "yb": 1.246e-2,
        "ye": 2.8492e-6, "ymu": 6.0022e-4, "ytau": 1.0196e-2,
    },
    "100_TeV": {
        "yu": 5.07e-6, "yc": 2.57e-3, "yt": 0.7309,
        "yd": 1.12e-5, "ys": 2.22e-4, "yb": 1.130e-2,
        "ye": 2.8637e-6, "ymu": 6.0327e-4, "ytau": 1.0248e-2,
    },
    "1e7_GeV": {
        "yu": 4.39e-6, "yc": 2.22e-3, "yt": 0.6462,
        "yd": 0.97e-5, "ys": 1.93e-4, "yb": 0.961e-2,
        "ye": 2.8681e-6, "ymu": 6.0421e-4, "ytau": 1.0264e-2,
    },
    "1e9_GeV": {
        "yu": 3.90e-6, "yc": 1.98e-3, "yt": 0.5838,
        "yd": 0.87e-5, "ys": 1.72e-4, "yb": 0.844e-2,
        "ye": 2.8508e-6, "ymu": 6.0056e-4, "ytau": 1.0201e-2,
    },
    "1e12_GeV": {
        "yu": 3.37e-6, "yc": 1.71e-3, "yt": 0.5138,
        "yd": 0.75e-5, "ys": 1.50e-4, "yb": 0.719e-2,
        "ye": 2.7978e-6, "ymu": 5.8939e-4, "ytau": 1.0012e-2,
    },
    "1e16_GeV": {
        "yu": 2.87e-6, "yc": 1.45e-3, "yt": 0.4454,
        "yd": 0.65e-5, "ys": 1.29e-4, "yb": 0.606e-2,
        "ye": 2.6935e-6, "ymu": 5.6745e-4, "ytau": 0.9639e-2,
    },
}


def compute_6_ratios_from_masses(m):
    return [
        m["c"] / m["u"],
        m["t"] / m["c"],
        m["s"] / m["d"],
        m["b"] / m["s"],
        m["mu"] / m["e"],
        m["tau"] / m["mu"],
    ]


def compute_6_ratios_from_yukawa(y):
    return [
        y["yc"] / y["yu"],
        y["yt"] / y["yc"],
        y["ys"] / y["yd"],
        y["yb"] / y["ys"],
        y["ymu"] / y["ye"],
        y["ytau"] / y["ymu"],
    ]


def generate_candidates(max_val):
    """Generate all candidates 2^a × p for odd primes p,
    up to max_val. Returns sorted list of (value, a, p)."""
    candidates = []
    for a in SUBSTRATE_POWERS:
        base = 2 ** a
        p = 3
        while True:
            val = base * p
            if val > max_val:
                break
            candidates.append((val, a, p))
            p = nextprime(p)
    candidates.sort()
    return candidates


def find_best_match(target, candidates):
    """Find the candidate closest to target.
    Returns (value, a, p, residual%)."""
    best_res = 999.0
    best = None
    for val, a, p in candidates:
        res = abs((target - val) / target * 100)
        if res < best_res:
            best_res = res
            best = (val, a, p)
    return best[0], best[1], best[2], best_res


def mean_substrate_residual(ratios, candidates):
    """Mean best-fit residual across a list of ratios."""
    total = sum(find_best_match(r, candidates)[3] for r in ratios)
    return total / len(ratios)


def count_distinct_primes(ratios, candidates):
    """Count distinct primes used across a list of ratios."""
    primes = set()
    for r in ratios:
        _, _, p, _ = find_best_match(r, candidates)
        primes.add(p)
    return len(primes)


def format_match(a, p):
    if a == 0:
        return f"{p}"
    elif a == 1:
        return f"2×{p}"
    elif a == 2:
        return f"2²×{p}"
    else:
        return f"2³×{p}"


if __name__ == "__main__":
    t0 = time.time()
    rng = random.Random(SEED)

    print("=" * 90)
    print("FERMION MASS PRIMES TEST 12")
    print("Hierarchy-matched single-prime substrate null test, multi-scale comparison")
    print("=" * 90)
    print(f"Ratio range: [{RATIO_MIN}, {RATIO_MAX}] (log-uniform)")
    print(f"Grammar: 2^a × p, a in {SUBSTRATE_POWERS}, p = single odd prime")
    print(f"Random trials: {N_TRIALS}")
    print(f"Seed: {SEED}")
    print(f"Metric: mean |ratio - best(2^a × p)| / ratio × 100")

    # ── Build all real scales ──
    all_scales = []
    pdg_ratios = compute_6_ratios_from_masses(PDG_MASSES)
    all_scales.append(("PDG_mixed", pdg_ratios))
    for label, y in ANTUSCH_SCALES.items():
        all_scales.append((label, compute_6_ratios_from_yukawa(y)))

    # Pre-generate candidates
    max_ratio = max(r for _, ratios in all_scales for r in ratios)
    candidates = generate_candidates(max_ratio * 2)
    print(f"Candidates: {len(candidates)} (2^a × p up to {max_ratio*2:.0f})")

    # ── Real data at each scale ──
    print(f"\n{'━' * 90}")
    print("REAL FERMION DATA — SINGLE-PRIME SUBSTRATE RESIDUALS AT EACH SCALE")
    print(f"{'━' * 90}")

    print(f"\n  {'Scale':<12} {'Mean res':>9}", end="")
    for name in RATIO_NAMES:
        print(f"  {name:>14}", end="")
    print(f"  {'#primes':>8}")
    print(f"  {'─'*12} {'─'*9}", end="")
    for _ in RATIO_NAMES:
        print(f"  {'─'*14}", end="")
    print(f"  {'─'*8}")

    real_mean_residuals = {}
    for label, ratios in all_scales:
        mean_res = mean_substrate_residual(ratios, candidates)
        n_primes = count_distinct_primes(ratios, candidates)
        real_mean_residuals[label] = mean_res
        print(f"  {label:<12} {mean_res:>8.4f}%", end="")
        for r in ratios:
            val, a, p, res = find_best_match(r, candidates)
            form = format_match(a, p)
            print(f"  {form:>8}({res:.2f}%)", end="")
        print(f"  {n_primes:>8}")

    # ── Generate random trials ──
    print(f"\n{'━' * 90}")
    print(f"GENERATING {N_TRIALS} RANDOM TRIALS")
    print(f"{'━' * 90}")

    random_residuals = []
    random_n_primes = []

    for trial in range(N_TRIALS):
        ratios = [math.exp(rng.uniform(LOG_RATIO_MIN, LOG_RATIO_MAX))
                  for _ in range(6)]
        mean_res = mean_substrate_residual(ratios, candidates)
        n_primes = count_distinct_primes(ratios, candidates)
        random_residuals.append(mean_res)
        random_n_primes.append(n_primes)

    random_residuals_sorted = sorted(random_residuals)
    random_n_primes_sorted = sorted(random_n_primes)

    # ── Random distribution summary ──
    print(f"\n  Random mean residual distribution:")
    print(f"    Mean:   {sum(random_residuals)/len(random_residuals):.4f}%")
    print(f"    Median: {random_residuals_sorted[len(random_residuals_sorted)//2]:.4f}%")
    for pct in [5, 10, 25, 75, 90, 95]:
        idx = int(N_TRIALS * pct / 100)
        print(f"    {pct}th pctile: {random_residuals_sorted[idx]:.4f}%")
    print(f"    Min:    {min(random_residuals):.4f}%")
    print(f"    Max:    {max(random_residuals):.4f}%")

    print(f"\n  Random distinct-prime count distribution:")
    print(f"    Mean:   {sum(random_n_primes)/len(random_n_primes):.1f}")
    print(f"    Median: {random_n_primes_sorted[len(random_n_primes_sorted)//2]}")
    print(f"    Min:    {min(random_n_primes)}")
    print(f"    Max:    {max(random_n_primes)}")
    counts = {}
    for n in random_n_primes:
        counts[n] = counts.get(n, 0) + 1
    for n in sorted(counts.keys()):
        print(f"    {n} primes: {counts[n]}/{N_TRIALS}")

    # ── Comparison: real vs random at each scale ──
    print(f"\n{'━' * 90}")
    print("COMPARISON: REAL vs RANDOM AT EACH SCALE")
    print(f"{'━' * 90}")

    print(f"\n  {'Scale':<12} {'Real res':>9} {'Rand mean':>10} {'Rand med':>9} "
          f"{'Pctile':>8} {'#worse':>8} {'Real #p':>8} {'Rand #p':>8}")
    print(f"  {'─'*12} {'─'*9} {'─'*10} {'─'*9} {'─'*8} {'─'*8} {'─'*8} {'─'*8}")

    rand_mean = sum(random_residuals) / len(random_residuals)
    rand_median = random_residuals_sorted[len(random_residuals_sorted) // 2]
    rand_mean_np = sum(random_n_primes) / len(random_n_primes)

    percentiles = []
    for label, ratios in all_scales:
        real_res = real_mean_residuals[label]
        real_np = count_distinct_primes(ratios, candidates)
        count_worse = sum(1 for r in random_residuals if r > real_res)
        pctile = count_worse / N_TRIALS * 100
        percentiles.append(pctile)
        print(f"  {label:<12} {real_res:>8.4f}% {rand_mean:>9.4f}% "
              f"{rand_median:>8.4f}% {pctile:>7.1f}% {count_worse:>5}/{N_TRIALS} "
              f"{real_np:>8} {rand_mean_np:>7.1f}")

    # ── Summary across all scales ──
    print(f"\n{'━' * 90}")
    print("SUMMARY ACROSS ALL SCALES")
    print(f"{'━' * 90}")

    real_res_values = list(real_mean_residuals.values())

    print(f"\n  Real mean residual range: "
          f"{min(real_res_values):.4f}% – {max(real_res_values):.4f}%")
    print(f"  Random mean residual:     {rand_mean:.4f}%")
    print(f"  Random median residual:   {rand_median:.4f}%")
    print(f"\n  Real percentile range: "
          f"{min(percentiles):.1f}% – {max(percentiles):.1f}%")
    print(f"  Mean percentile across scales: "
          f"{sum(percentiles)/len(percentiles):.1f}%")

    scales_above_median = sum(1 for r in real_res_values if r < rand_median)
    print(f"\n  Scales where real < random median: "
          f"{scales_above_median}/{len(all_scales)}")
    print(f"  Worst scale percentile: {min(percentiles):.1f}%")
    print(f"  Best scale percentile:  {max(percentiles):.1f}%")

    total_elapsed = time.time() - t0
    print(f"\nTotal runtime: {total_elapsed:.1f}s")
    print("Done.")
