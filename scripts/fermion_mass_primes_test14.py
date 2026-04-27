"""
fermion_mass_primes_test14.py

Cross-family same-generation boundary test.

Tests whether the strict single-prime substrate grammar (2^a × p) from
test11 also describes same-generation cross-family mass ratios.

This is a boundary test, not a grammar extension.

See scripts/README.md for full documentation.
"""

import time
from sympy import nextprime

# ── Within-family ratios (test11 comparator) ──
WITHIN_FAMILY_NAMES = ["mc/mu", "mt/mc", "ms/md", "mb/ms", "mmu/me", "mtau/mmu"]

# ── Cross-family same-generation ratios ──
# (label, particle1_key, particle2_key) — label is PDG-form heavier/lighter
# R is computed as max(m1, m2) / min(m1, m2); orientation recorded separately
CROSS_FAMILY_PAIRS = [
    ("d/u",   "d",  "u"),
    ("c/s",   "c",  "s"),
    ("t/b",   "t",  "b"),
    ("u/e",   "u",  "e"),
    ("c/mu",  "c",  "mu"),
    ("t/tau", "t",  "tau"),
    ("d/e",   "d",  "e"),
    ("mu/s",  "mu", "s"),
    ("b/tau", "b",  "tau"),
]
CROSS_FAMILY_NAMES = [p[0] for p in CROSS_FAMILY_PAIRS]

SUBSTRATE_POWERS = [0, 1, 2, 3]
GRAMMAR_FLOOR = 3

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


def normalize_scale(label, raw):
    """Return mass-equivalent dict keyed by particle name (u, c, t, d, s, b, e, mu, tau)."""
    if label == "PDG_mixed":
        return raw
    return {
        "u": raw["yu"], "c": raw["yc"], "t": raw["yt"],
        "d": raw["yd"], "s": raw["ys"], "b": raw["yb"],
        "e": raw["ye"], "mu": raw["ymu"], "tau": raw["ytau"],
    }


def compute_within_family_ratios(masses):
    return [
        masses["c"]   / masses["u"],
        masses["t"]   / masses["c"],
        masses["s"]   / masses["d"],
        masses["b"]   / masses["s"],
        masses["mu"]  / masses["e"],
        masses["tau"] / masses["mu"],
    ]


def compute_cross_family_ratios(masses):
    """Return list of (R, orientation) in CROSS_FAMILY_PAIRS order.
    R = max/min; orientation = 'k1>k2' or 'k2>k1' depending on which is heavier."""
    out = []
    for label, k1, k2 in CROSS_FAMILY_PAIRS:
        m1, m2 = masses[k1], masses[k2]
        if m1 >= m2:
            out.append((m1 / m2, f"{k1}>{k2}"))
        else:
            out.append((m2 / m1, f"{k2}>{k1}"))
    return out


def generate_candidates(max_val):
    """All 2^a × p with a in SUBSTRATE_POWERS and p odd prime, up to max_val."""
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
    """Nearest-candidate fit. Returns (val, a, p, residual%).
    Ties go to the smaller candidate (sorted order)."""
    best_res = float("inf")
    best = None
    for val, a, p in candidates:
        res = abs((target - val) / target * 100)
        if res < best_res:
            best_res = res
            best = (val, a, p)
    return best[0], best[1], best[2], best_res


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
    n_pairs = len(CROSS_FAMILY_PAIRS)

    print("=" * 100)
    print("FERMION MASS PRIMES TEST 14")
    print("Cross-family same-generation boundary test — grammar: 2^a × p")
    print("=" * 100)
    print(f"Substrate powers: a in {SUBSTRATE_POWERS}")
    print(f"p = single odd prime")
    print(f"Grammar floor: R = {GRAMMAR_FLOOR} (smallest candidate, a=0 p=3)")
    print(f"Cross-family pairs: {n_pairs}")
    print(f"No composite products, no prime powers, no denominators")

    # ── Build all scales ──
    all_scales = [("PDG_mixed", normalize_scale("PDG_mixed", PDG_MASSES))]
    for label, raw in ANTUSCH_SCALES.items():
        all_scales.append((label, normalize_scale(label, raw)))
    n_scales = len(all_scales)

    # ── Pre-compute ratios at each scale ──
    cross_per_scale = {}
    within_per_scale = {}
    for label, masses in all_scales:
        cross_per_scale[label] = compute_cross_family_ratios(masses)
        within_per_scale[label] = compute_within_family_ratios(masses)

    # ── Generate candidates up to 2× global max ──
    all_R = []
    for label, _ in all_scales:
        all_R.extend([R for R, _ in cross_per_scale[label]])
        all_R.extend(within_per_scale[label])
    max_R = max(all_R)
    candidates = generate_candidates(max_R * 2)
    print(f"Candidates generated: {len(candidates)} (2^a × p up to {max_R*2:.0f})")

    # ── Cross-family ratio values table ──
    print(f"\n{'─' * 110}")
    print("CROSS-FAMILY RATIO VALUES AT EACH SCALE (heavier/lighter, * = R<3)")
    print(f"{'─' * 110}")
    header = f"\n{'Scale':<12}"
    for name in CROSS_FAMILY_NAMES:
        header += f" {name:>10}"
    print(header)
    print("─" * (12 + 11 * n_pairs))
    for label, _ in all_scales:
        row = f"{label:<12}"
        for R, _ in cross_per_scale[label]:
            mark = "*" if R < GRAMMAR_FLOOR else " "
            row += f" {R:>9.3f}{mark}"
        print(row)

    # ── Run the fits ──
    # cross_results[scale] = [(name, R, orient, below, val, a, p, res), ...]
    cross_results = {}
    for label, _ in all_scales:
        scale_results = []
        for (name, _, _), (R, orient) in zip(CROSS_FAMILY_PAIRS,
                                              cross_per_scale[label]):
            below = R < GRAMMAR_FLOOR
            val, a, p, res = find_best_match(R, candidates)
            scale_results.append((name, R, orient, below, val, a, p, res))
        cross_results[label] = scale_results

    # within_results[scale] = [(name, R, val, a, p, res), ...]
    within_results = {}
    for label, _ in all_scales:
        scale_results = []
        for name, R in zip(WITHIN_FAMILY_NAMES, within_per_scale[label]):
            val, a, p, res = find_best_match(R, candidates)
            scale_results.append((name, R, val, a, p, res))
        within_results[label] = scale_results

    # ── Per-scale detailed output ──
    for label, _ in all_scales:
        print(f"\n{'━' * 100}")
        print(f"SCALE: {label}")
        print(f"{'━' * 100}")
        print(f"  {'Ratio':<8} {'R':>10} {'orient':<10} {'<3':>4} "
              f"{'best':>6} {'form':<10} {'a':>3} {'p':>5} {'residual':>11}")
        print(f"  {'─'*8} {'─'*10} {'─'*10} {'─'*4} "
              f"{'─'*6} {'─'*10} {'─'*3} {'─'*5} {'─'*11}")
        for name, R, orient, below, val, a, p, res in cross_results[label]:
            below_str = "yes" if below else ""
            form = format_match(a, p)
            print(f"  {name:<8} {R:>10.4f} {orient:<10} {below_str:>4} "
                  f"{val:>6} {form:<10} {a:>3} {p:>5} {res:>10.4f}%")
        residuals = [r[7] for r in cross_results[label]]
        below_n = sum(1 for r in cross_results[label] if r[3])
        print(f"\n  Mean residual: {sum(residuals)/len(residuals):.4f}%   "
              f"Below floor: {below_n}/{n_pairs}")

    # ── Grand summary table ──
    print(f"\n{'━' * 160}")
    print("GRAND SUMMARY — CROSS-FAMILY BEST MATCH AT EACH SCALE")
    print(f"{'━' * 160}")
    header = f"\n{'Scale':<12}"
    for name in CROSS_FAMILY_NAMES:
        header += f" {name:>14}"
    header += f" {'Mean':>9} {'<3':>5}"
    print(header)
    sep = f"{'─'*12}"
    for _ in CROSS_FAMILY_NAMES:
        sep += f" {'─'*14}"
    sep += f" {'─'*9} {'─'*5}"
    print(sep)
    for label, _ in all_scales:
        row = f"{label:<12}"
        for r in cross_results[label]:
            _name, _R, _orient, below, val, a, p, res = r
            form = format_match(a, p)
            mark = "*" if below else " "
            cell = f"{form}({res:.1f}%){mark}"
            row += f" {cell:>14}"
        residuals = [r[7] for r in cross_results[label]]
        below_n = sum(1 for r in cross_results[label] if r[3])
        row += f" {sum(residuals)/len(residuals):>8.4f}% {below_n:>3}/{n_pairs}"
        print(row)

    # ── Per-ratio analysis ──
    print(f"\n{'━' * 100}")
    print("PER-RATIO ANALYSIS")
    print(f"{'━' * 100}")
    for name in CROSS_FAMILY_NAMES:
        cells = []  # [(label, result), ...]
        for label, _ in all_scales:
            for r in cross_results[label]:
                if r[0] == name:
                    cells.append((label, r))
                    break

        residuals = [c[1][7] for c in cells]
        primes = sorted(set(c[1][6] for c in cells))
        a_values = sorted(set(c[1][5] for c in cells))
        below_n = sum(1 for c in cells if c[1][3])
        orient_seq = [c[1][2] for c in cells]
        unity_crossing = len(set(orient_seq)) > 1

        prime_class = (f"LOCKED  at p={primes[0]}" if len(primes) == 1
                       else f"WALKING through {primes}")
        sub_class = (f"LOCKED  at a={a_values[0]}" if len(a_values) == 1
                     else f"VARIABLE a={a_values}")

        print(f"\n  {name}")
        print(f"    Mean residual:   {sum(residuals)/len(residuals):.4f}%")
        print(f"    Below floor:     {below_n}/{n_scales}")
        print(f"    Prime:           {prime_class}")
        print(f"    Substrate:       {sub_class}")
        if unity_crossing:
            print(f"    UNITY-CROSSING:  orientations seen: "
                  f"{sorted(set(orient_seq))}")
            for i in range(1, len(cells)):
                prev_label, prev_r = cells[i-1]
                cur_label, cur_r = cells[i]
                if prev_r[2] != cur_r[2]:
                    print(f"                     bracketed between "
                          f"{prev_label} ({prev_r[2]}, R={prev_r[1]:.3f}) "
                          f"and {cur_label} ({cur_r[2]}, R={cur_r[1]:.3f})")

    # ── Per-scale comparison: within-family vs cross-family ──
    print(f"\n{'━' * 100}")
    print("PER-SCALE MEAN RESIDUAL: WITHIN-FAMILY (test11) vs CROSS-FAMILY (test14)")
    print(f"{'━' * 100}")
    print(f"\n  {'Scale':<12} {'Within (6)':>12} {'Cross (9)':>12} {'<3 (cross)':>12}")
    print(f"  {'─'*12} {'─'*12} {'─'*12} {'─'*12}")
    sum_within = 0.0
    sum_cross = 0.0
    sum_below = 0
    for label, _ in all_scales:
        wres = [r[5] for r in within_results[label]]
        cres = [r[7] for r in cross_results[label]]
        bn = sum(1 for r in cross_results[label] if r[3])
        wm = sum(wres) / len(wres)
        cm = sum(cres) / len(cres)
        sum_within += wm
        sum_cross += cm
        sum_below += bn
        print(f"  {label:<12} {wm:>11.4f}% {cm:>11.4f}% {bn:>9}/{n_pairs}")
    print(f"  {'─'*12} {'─'*12} {'─'*12} {'─'*12}")
    print(f"  {'Overall':<12} {sum_within/n_scales:>11.4f}% "
          f"{sum_cross/n_scales:>11.4f}% {sum_below:>9}/{n_pairs*n_scales}")

    total_elapsed = time.time() - t0
    print(f"\nTotal runtime: {total_elapsed:.1f}s")
    print("Done.")
