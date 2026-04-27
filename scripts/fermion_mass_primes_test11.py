"""
fermion_mass_primes_test11.py

Multi-scale single-prime substrate test (minimal grammar).

Tests whether fermion mass ratios are best described as 2^a × p
where p is a single odd prime and a ∈ {0,1,2,3}.

See scripts/README.md for full documentation.
"""

import time
from sympy import nextprime

RATIO_NAMES = ["mc/mu", "mt/mc", "ms/md", "mb/ms", "mmu/me", "mtau/mmu"]
SUBSTRATE_POWERS = [0, 1, 2, 3]

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
    Returns (value, a, p, residual%, distance)."""
    best_res = 999.0
    best = None
    for val, a, p in candidates:
        res = abs((target - val) / target * 100)
        if res < best_res:
            best_res = res
            best = (val, a, p)
    if best is None:
        return None, None, None, 999.0, 999.0
    dist = target - best[0]
    return best[0], best[1], best[2], best_res, dist


def format_match(val, a, p):
    """Format a match as a string."""
    if a == 0:
        return f"{p}"
    elif a == 1:
        return f"2 × {p}"
    else:
        return f"2^{a} × {p}"


if __name__ == "__main__":
    t0 = time.time()

    print("=" * 100)
    print("FERMION MASS PRIMES TEST 11")
    print("Multi-scale single-prime substrate test — minimal grammar: 2^a × p")
    print("=" * 100)
    print(f"Substrate powers: a in {SUBSTRATE_POWERS}")
    print(f"p = any single odd prime")
    print(f"No composite products, no prime powers, no denominators")

    # ── Build all scales ──
    all_scales = []
    pdg_ratios = compute_6_ratios_from_masses(PDG_MASSES)
    all_scales.append(("PDG_mixed", pdg_ratios))
    for label, y in ANTUSCH_SCALES.items():
        all_scales.append((label, compute_6_ratios_from_yukawa(y)))

    # Pre-generate candidates up to 2× max ratio
    max_ratio = max(r for _, ratios in all_scales for r in ratios)
    candidates = generate_candidates(max_ratio * 2)
    print(f"Candidates generated: {len(candidates)} (2^a × p up to {max_ratio*2:.0f})")

    # ── Ratio values at each scale ──
    print(f"\n{'─' * 100}")
    print("RATIO VALUES AT EACH SCALE")
    print(f"{'─' * 100}")
    print(f"\n{'Scale':<12} " + " ".join(f"{r:>12}" for r in RATIO_NAMES))
    print("─" * 86)
    for label, ratios in all_scales:
        print(f"{label:<12} " + " ".join(f"{v:>12.4f}" for v in ratios))

    # ── Per-scale analysis ──
    # Store all results for cross-scale summaries
    all_results = {}  # {scale: [(ratio_name, target, val, a, p, res, dist), ...]}

    for label, ratios in all_scales:
        print(f"\n{'━' * 100}")
        print(f"SCALE: {label}")
        print(f"{'━' * 100}")
        print(f"  {'Ratio':<12} {'Target':>10} {'Match':>8} {'= form':<20} "
              f"{'a':>3} {'p':>6} {'Residual':>10} {'Distance':>10}")
        print(f"  {'─'*12} {'─'*10} {'─'*8} {'─'*20} {'─'*3} {'─'*6} "
              f"{'─'*10} {'─'*10}")

        scale_results = []
        for name, target in zip(RATIO_NAMES, ratios):
            val, a, p, res, dist = find_best_match(target, candidates)
            form = format_match(val, a, p)
            marker = "◄ bare" if a == 0 else ""
            print(f"  {name:<12} {target:>10.4f} {val:>8} {form:<20} "
                  f"{a:>3} {p:>6} {res:>9.4f}% {dist:>+10.4f} {marker}")
            scale_results.append((name, target, val, a, p, res, dist))

        all_results[label] = scale_results

        mean_res = sum(r[5] for r in scale_results) / len(scale_results)
        n_bare = sum(1 for r in scale_results if r[3] == 0)
        print(f"\n  Mean residual: {mean_res:.4f}%")
        print(f"  Bare prime (a=0): {n_bare}/6")

    # ── Grand summary table ──
    print(f"\n{'━' * 100}")
    print("GRAND SUMMARY — BEST MATCH AT EACH SCALE")
    print(f"{'━' * 100}")
    print(f"\n{'Scale':<12}", end="")
    for name in RATIO_NAMES:
        print(f" {name:>16}", end="")
    print(f" {'Mean res':>9} {'#bare':>6}")
    print(f"{'─'*12}", end="")
    for _ in RATIO_NAMES:
        print(f" {'─'*16}", end="")
    print(f" {'─'*9} {'─'*6}")

    for label, ratios in all_scales:
        results = all_results[label]
        print(f"{label:<12}", end="")
        for name, target, val, a, p, res, dist in results:
            if a == 0:
                form_short = f"{p}"
            elif a == 1:
                form_short = f"2×{p}"
            elif a == 2:
                form_short = f"2²×{p}"
            else:
                form_short = f"2³×{p}"
            print(f" {form_short:>9}({res:.2f}%)", end="")
        mean_res = sum(r[5] for r in results) / len(results)
        n_bare = sum(1 for r in results if r[3] == 0)
        print(f" {mean_res:>8.4f}% {n_bare:>5}/6")

    # ── Substrate power distribution ──
    print(f"\n{'━' * 100}")
    print("SUBSTRATE POWER USAGE (a) — PER RATIO ACROSS SCALES")
    print(f"{'━' * 100}")
    for name in RATIO_NAMES:
        a_counts = {0: 0, 1: 0, 2: 0, 3: 0}
        for label in [s[0] for s in all_scales]:
            for r in all_results[label]:
                if r[0] == name:
                    a_counts[r[3]] += 1
        total = sum(a_counts.values())
        parts = [f"a={a}:{c}/{total}" for a, c in sorted(a_counts.items()) if c > 0]
        print(f"  {name:<12} {', '.join(parts)}")

    # ── Prime stability per ratio ──
    print(f"\n{'━' * 100}")
    print("PRIME STABILITY — PER RATIO ACROSS SCALES")
    print(f"{'━' * 100}")
    for name in RATIO_NAMES:
        primes_seen = {}
        for label in [s[0] for s in all_scales]:
            for r in all_results[label]:
                if r[0] == name:
                    p = r[4]
                    primes_seen[p] = primes_seen.get(p, 0) + 1
        residuals = []
        for label in [s[0] for s in all_scales]:
            for r in all_results[label]:
                if r[0] == name:
                    residuals.append(r[5])

        if len(primes_seen) == 1:
            p = list(primes_seen.keys())[0]
            print(f"  {name:<12} LOCKED  p={p:<5} "
                  f"mean_res={sum(residuals)/len(residuals):.4f}%")
        else:
            primes_list = sorted(primes_seen.keys())
            print(f"  {name:<12} WALKING p={primes_list} "
                  f"mean_res={sum(residuals)/len(residuals):.4f}%")

    # ── Sector comparison ──
    print(f"\n{'━' * 100}")
    print("SECTOR COMPARISON — QUARK vs LEPTON")
    print(f"{'━' * 100}")
    quark_names = ["mc/mu", "mt/mc", "ms/md", "mb/ms"]
    lepton_names = ["mmu/me", "mtau/mmu"]

    for sector, names_list in [("Quark", quark_names), ("Lepton", lepton_names)]:
        a_counts = {0: 0, 1: 0, 2: 0, 3: 0}
        residuals = []
        for label in [s[0] for s in all_scales]:
            for r in all_results[label]:
                if r[0] in names_list:
                    a_counts[r[3]] += 1
                    residuals.append(r[5])
        total = sum(a_counts.values())
        print(f"\n  {sector} sector ({len(names_list)} ratios × {len(all_scales)} scales"
              f" = {total} measurements):")
        print(f"    Mean residual: {sum(residuals)/len(residuals):.4f}%")
        for a in SUBSTRATE_POWERS:
            if a_counts[a] > 0:
                print(f"    a={a}: {a_counts[a]}/{total} "
                      f"({a_counts[a]/total*100:.1f}%)")

    # ── Mean residual per scale ──
    print(f"\n{'━' * 100}")
    print("MEAN RESIDUAL — PER SCALE")
    print(f"{'━' * 100}")
    for label, ratios in all_scales:
        results = all_results[label]
        mean_res = sum(r[5] for r in results) / len(results)
        n_bare = sum(1 for r in results if r[3] == 0)
        print(f"  {label:<12} {mean_res:.4f}%  (bare: {n_bare}/6)")

    total_elapsed = time.time() - t0
    print(f"\nTotal runtime: {total_elapsed:.1f}s")
    print("Done.")
