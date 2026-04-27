"""
fermion_mass_primes_test15.py

Sector prime-set examination.

Examines the common-scale prime vocabularies produced by the test11 grammar,
split by charged-fermion sector. No new fits are computed beyond what test11
already produced; T5 alone re-derives nearest-fit per cell to inspect local
candidate density.

Six sub-tests:
  T0 — Occupancy table (per-prime frequency, ratios, substrates, lock status)
  T1 — Skip-pattern / window sparsity
  T2 — Disjointness null (window-matched and pooled-window)
  T3 — Algebraic structure scan
  T4 — Gap structure within sector
  T5 — Local candidate-density diagnostic

See scripts/README.md for full documentation.
"""

import bisect
import math
import random
import time
from sympy import nextprime, isprime, primerange

# ── Test11 within-family ratios ──
RATIO_NAMES = ["mc/mu", "mt/mc", "ms/md", "mb/ms", "mmu/me", "mtau/mmu"]

# Sector mapping
RATIO_SECTOR = {
    "mc/mu": "up",
    "mt/mc": "up",
    "ms/md": "down",
    "mb/ms": "down",
    "mmu/me": "lepton",
    "mtau/mmu": "lepton",
}
SECTORS = ["up", "down", "lepton"]

SUBSTRATE_POWERS = [0, 1, 2, 3]

# ── Antusch et al. Table 2 (common scales only — PDG excluded) ──
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

N_NULL_TRIALS_T2 = 100000
N_NULL_TRIALS_T3 = 10000
N_NULL_TRIALS_T4 = 10000
N_NULL_TRIALS_T5 = 10000
SEED = 42


# ─────────────────────────────────────────────────────────────────────
# Grammar utilities (replicating test11)
# ─────────────────────────────────────────────────────────────────────

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
    """All 2^a × p candidates with a in SUBSTRATE_POWERS and p odd prime,
    up to max_val. Returns sorted list of (val, a, p)."""
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
    """Nearest-candidate fit. Returns (val, a, p, residual%)."""
    best_res = float("inf")
    best = None
    for val, a, p in candidates:
        res = abs((target - val) / target * 100)
        if res < best_res:
            best_res = res
            best = (val, a, p)
    return best[0], best[1], best[2], best_res


def fmt_form(a, p):
    if a == 0:
        return f"{p}"
    elif a == 1:
        return f"2×{p}"
    elif a == 2:
        return f"2²×{p}"
    else:
        return f"2³×{p}"


# ─────────────────────────────────────────────────────────────────────
# Build per-cell results from common scales (re-running test11 fits for
# transparency; this is cheap and avoids hard-coding values from a log)
# ─────────────────────────────────────────────────────────────────────

def build_cells():
    """Return list of dicts: one per (ratio, scale) cell at common scales.
    Each: {ratio, sector, scale, target, val, a, p, residual}."""
    # max ratio across all common scales
    max_target = 0.0
    for label, y in ANTUSCH_SCALES.items():
        for r in compute_6_ratios_from_yukawa(y):
            max_target = max(max_target, r)
    candidates = generate_candidates(max_target * 2)

    cells = []
    for label, y in ANTUSCH_SCALES.items():
        ratios = compute_6_ratios_from_yukawa(y)
        for name, target in zip(RATIO_NAMES, ratios):
            val, a, p, res = find_best_match(target, candidates)
            cells.append({
                "ratio": name,
                "sector": RATIO_SECTOR[name],
                "scale": label,
                "target": target,
                "val": val,
                "a": a,
                "p": p,
                "residual": res,
            })
    return cells, candidates


# ─────────────────────────────────────────────────────────────────────
# T0 — Occupancy table
# ─────────────────────────────────────────────────────────────────────

def t0_occupancy(cells):
    print("=" * 90)
    print("T0 — OCCUPANCY TABLE")
    print("=" * 90)

    for sector in SECTORS:
        sector_cells = [c for c in cells if c["sector"] == sector]
        # collect per-prime info
        per_prime = {}  # p -> {count, ratios: set, substrates: set}
        for c in sector_cells:
            p = c["p"]
            if p not in per_prime:
                per_prime[p] = {"count": 0, "ratios": set(), "substrates": set()}
            per_prime[p]["count"] += 1
            per_prime[p]["ratios"].add(c["ratio"])
            per_prime[p]["substrates"].add(c["a"])

        # per-ratio lock/walk
        ratio_primes = {}  # ratio -> set of primes
        for c in sector_cells:
            ratio_primes.setdefault(c["ratio"], set()).add(c["p"])

        n_scales = len({c["scale"] for c in sector_cells})

        print(f"\n{sector.upper()} sector ({len(sector_cells)} cells, "
              f"{n_scales} common scales):")
        print(f"  {'prime':>6} {'count':>6} {'ratios':<22} {'substrates':<12} "
              f"{'lock status':<20}")
        print(f"  {'─'*6} {'─'*6} {'─'*22} {'─'*12} {'─'*20}")

        for p in sorted(per_prime.keys()):
            info = per_prime[p]
            ratios_str = ",".join(sorted(info["ratios"]))
            subs_str = ",".join(f"a={a}" for a in sorted(info["substrates"]))
            # lock status: locked if all cells of all visiting ratios use this prime
            lock_parts = []
            for r in sorted(info["ratios"]):
                if len(ratio_primes[r]) == 1:
                    lock_parts.append(f"{r}:locked")
                else:
                    lock_parts.append(f"{r}:walking")
            lock_str = ",".join(lock_parts)
            print(f"  {p:>6} {info['count']:>6} {ratios_str:<22} {subs_str:<12} "
                  f"{lock_str:<20}")

        # substrate counts for sector
        sub_counts = {a: 0 for a in SUBSTRATE_POWERS}
        for c in sector_cells:
            sub_counts[c["a"]] += 1
        total = sum(sub_counts.values())
        sub_str = ", ".join(
            f"a={a}: {sub_counts[a]}/{total} ({sub_counts[a]/total*100:.0f}%)"
            for a in SUBSTRATE_POWERS if sub_counts[a] > 0
        )
        print(f"  Substrate counts: {sub_str}")


# ─────────────────────────────────────────────────────────────────────
# T1 — Skip pattern / window sparsity
# ─────────────────────────────────────────────────────────────────────

def t1_skip_pattern(cells):
    print("\n" + "=" * 90)
    print("T1 — SKIP-PATTERN / WINDOW SPARSITY")
    print("=" * 90)

    sector_primes = {s: set() for s in SECTORS}
    for c in cells:
        sector_primes[c["sector"]].add(c["p"])

    print(f"\n{'sector':<8} {'primes':<6} {'window':<14} "
          f"{'avail':>6} {'visited':>8} {'skipped':>8} {'skip %':>8}")
    print(f"{'─'*8} {'─'*6} {'─'*14} {'─'*6} {'─'*8} {'─'*8} {'─'*8}")

    t1_results = {}
    for sector in SECTORS:
        primes = sorted(sector_primes[sector])
        if len(primes) < 1:
            continue
        lo, hi = primes[0], primes[-1]
        available = list(primerange(lo, hi + 1))  # odd primes; primerange covers all
        n_avail = len(available)
        n_visited = len(primes)
        n_skipped = n_avail - n_visited
        skip_frac = n_skipped / n_avail if n_avail > 0 else 0.0
        window_str = f"[{lo},{hi}]"
        print(f"{sector:<8} {n_visited:<6} {window_str:<14} "
              f"{n_avail:>6} {n_visited:>8} {n_skipped:>8} {skip_frac*100:>7.1f}%")
        t1_results[sector] = {
            "lo": lo, "hi": hi, "avail": n_avail,
            "visited": n_visited, "skipped": n_skipped, "skip_frac": skip_frac,
            "available_primes": available, "visited_primes": primes,
        }

    # show what was skipped
    print()
    for sector in SECTORS:
        if sector not in t1_results:
            continue
        r = t1_results[sector]
        skipped_primes = [p for p in r["available_primes"] if p not in r["visited_primes"]]
        if len(skipped_primes) <= 30:
            sk_str = ", ".join(str(p) for p in skipped_primes)
        else:
            sk_str = (", ".join(str(p) for p in skipped_primes[:15]) +
                      f", ... ({len(skipped_primes)-15} more)")
        print(f"  {sector} skipped: {sk_str}")

    return t1_results


# ─────────────────────────────────────────────────────────────────────
# T2 — Disjointness null
# ─────────────────────────────────────────────────────────────────────

def overlap_total(set_a, set_b, set_c):
    """Total number of pairwise overlaps across three sets."""
    return (len(set_a & set_b) + len(set_a & set_c) + len(set_b & set_c))


def t2_disjointness_null(cells, t1_results, rng):
    print("\n" + "=" * 90)
    print("T2 — DISJOINTNESS NULL (window-matched and pooled-window)")
    print("=" * 90)

    sector_primes = {s: set() for s in SECTORS}
    for c in cells:
        sector_primes[c["sector"]].add(c["p"])

    actual_overlap_total = overlap_total(
        sector_primes["up"], sector_primes["down"], sector_primes["lepton"]
    )
    actual_pairs = {
        "up∩down": len(sector_primes["up"] & sector_primes["down"]),
        "up∩lepton": len(sector_primes["up"] & sector_primes["lepton"]),
        "down∩lepton": len(sector_primes["down"] & sector_primes["lepton"]),
    }

    print(f"\nActual cardinalities: up={len(sector_primes['up'])}, "
          f"down={len(sector_primes['down'])}, lepton={len(sector_primes['lepton'])}")
    print(f"Actual pairwise overlaps: {actual_pairs}")
    print(f"Actual total overlap: {actual_overlap_total}")

    # ── A. Window-matched null ──
    print(f"\n--- A. Window-matched null ({N_NULL_TRIALS_T2} trials) ---")
    window_pools = {
        s: list(primerange(t1_results[s]["lo"], t1_results[s]["hi"] + 1))
        for s in SECTORS
    }
    cardinalities = {s: len(sector_primes[s]) for s in SECTORS}

    overlaps_A = []
    for _ in range(N_NULL_TRIALS_T2):
        draws = {s: set(rng.sample(window_pools[s], cardinalities[s])) for s in SECTORS}
        ov = overlap_total(draws["up"], draws["down"], draws["lepton"])
        overlaps_A.append(ov)

    n_le_actual_A = sum(1 for o in overlaps_A if o <= actual_overlap_total)
    pctile_A = n_le_actual_A / N_NULL_TRIALS_T2 * 100
    n_zero_A = sum(1 for o in overlaps_A if o == 0)
    print(f"Window pools: up={len(window_pools['up'])}, "
          f"down={len(window_pools['down'])}, lepton={len(window_pools['lepton'])}")
    print(f"Trials with overlap ≤ actual ({actual_overlap_total}): "
          f"{n_le_actual_A}/{N_NULL_TRIALS_T2} ({pctile_A:.2f}%)")
    print(f"Trials with zero overlap: {n_zero_A}/{N_NULL_TRIALS_T2} "
          f"({n_zero_A/N_NULL_TRIALS_T2*100:.2f}%)")
    print(f"Mean random overlap: {sum(overlaps_A)/N_NULL_TRIALS_T2:.3f}")

    # ── B. Pooled-window null ──
    print(f"\n--- B. Pooled-window null ({N_NULL_TRIALS_T2} trials) ---")
    max_observed = max(p for s in SECTORS for p in sector_primes[s])
    min_observed = min(p for s in SECTORS for p in sector_primes[s])
    pooled = list(primerange(min_observed, max_observed + 1))
    print(f"Pool: odd primes in [{min_observed}, {max_observed}] → {len(pooled)} primes")

    overlaps_B = []
    for _ in range(N_NULL_TRIALS_T2):
        # draw three disjoint-or-not sets independently from same pool
        draws = {s: set(rng.sample(pooled, cardinalities[s])) for s in SECTORS}
        ov = overlap_total(draws["up"], draws["down"], draws["lepton"])
        overlaps_B.append(ov)

    n_le_actual_B = sum(1 for o in overlaps_B if o <= actual_overlap_total)
    pctile_B = n_le_actual_B / N_NULL_TRIALS_T2 * 100
    n_zero_B = sum(1 for o in overlaps_B if o == 0)
    print(f"Trials with overlap ≤ actual ({actual_overlap_total}): "
          f"{n_le_actual_B}/{N_NULL_TRIALS_T2} ({pctile_B:.2f}%)")
    print(f"Trials with zero overlap: {n_zero_B}/{N_NULL_TRIALS_T2} "
          f"({n_zero_B/N_NULL_TRIALS_T2*100:.2f}%)")
    print(f"Mean random overlap: {sum(overlaps_B)/N_NULL_TRIALS_T2:.3f}")

    print()
    print("Interpretation:")
    if pctile_B < 5:
        print("  Pooled-window null: actual disjointness is unusually clean (<5%).")
    elif pctile_B > 50:
        print("  Pooled-window null: actual disjointness is typical or expected.")
    else:
        print("  Pooled-window null: actual disjointness is mildly clean but not striking.")


# ─────────────────────────────────────────────────────────────────────
# T3 — Algebraic structure scan
# ─────────────────────────────────────────────────────────────────────

def is_sophie_germain(p):
    return isprime(2 * p + 1)


def is_twin_participant(p):
    return isprime(p - 2) or isprime(p + 2)


def is_mersenne(p):
    # p = 2^n − 1 for some integer n
    n = p + 1
    return n > 0 and (n & (n - 1)) == 0


def is_fermat(p):
    # p = 2^(2^n) + 1 for some non-negative integer n
    if p < 3:
        return False
    k = p - 1
    if (k & (k - 1)) != 0:
        return False
    # k is a power of 2; check if log2(k) is itself a power of 2
    log2k = k.bit_length() - 1
    return log2k > 0 and (log2k & (log2k - 1)) == 0 or log2k == 1 or log2k == 0


def is_primorial_pm1(p):
    # primorial p_n# = product of first n primes
    primorials = []
    primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    prod = 1
    for q in primes:
        prod *= q
        primorials.append(prod)
    return any((p == prod + 1 or p == prod - 1) for prod in primorials)


def is_form_n2_n_1(p):
    # p = n^2 + n + 1 for some non-negative integer n
    # → n = (-1 + sqrt(4p - 3)) / 2
    disc = 4 * p - 3
    if disc < 0:
        return False
    s = math.isqrt(disc)
    if s * s != disc:
        return False
    if (s - 1) % 2 != 0:
        return False
    n = (s - 1) // 2
    return n * n + n + 1 == p


def is_form_n2_1(p):
    # p = n^2 + 1
    if p < 2:
        return False
    n = math.isqrt(p - 1)
    return n * n + 1 == p


def is_form_4n2_1(p):
    # p = 4n^2 + 1
    if p < 5:
        return p == 5  # 4*1+1
    if (p - 1) % 4 != 0:
        return False
    k = (p - 1) // 4
    n = math.isqrt(k)
    return n * n == k


def algebraic_properties():
    """Return list of (name, predicate) for T3."""
    props = []

    # Modular residue classes
    for m in [3, 4, 6, 8, 12]:
        for r in range(1, m):
            if math.gcd(r, m) != 1:
                continue  # only residues coprime to m can host primes (beyond p|m)
            props.append((f"≡{r} mod {m}",
                          lambda p, r=r, m=m: p % m == r))

    # Eisenstein-prime split class (rational primes that split in Z[ω])
    props.append(("Eisenstein split (≡1 mod 3 or =3)",
                  lambda p: p == 3 or p % 3 == 1))
    # Gaussian-prime structure (odd primes that remain inert in Z[i])
    props.append(("Gaussian inert (≡3 mod 4)",
                  lambda p: p % 4 == 3))

    props.append(("Sophie Germain", is_sophie_germain))
    props.append(("Twin participant", is_twin_participant))
    props.append(("Mersenne (2^n − 1)", is_mersenne))
    props.append(("Fermat (2^(2^n) + 1)", is_fermat))
    props.append(("Primorial ±1", is_primorial_pm1))
    props.append(("n² + n + 1", is_form_n2_n_1))
    props.append(("n² + 1", is_form_n2_1))
    props.append(("4n² + 1", is_form_4n2_1))

    return props


def t3_algebraic_scan(cells, t1_results, rng):
    print("\n" + "=" * 90)
    print("T3 — ALGEBRAIC STRUCTURE SCAN")
    print("=" * 90)

    props = algebraic_properties()
    sector_primes = {s: set() for s in SECTORS}
    for c in cells:
        sector_primes[c["sector"]].add(c["p"])

    window_pools = {
        s: list(primerange(t1_results[s]["lo"], t1_results[s]["hi"] + 1))
        for s in SECTORS
    }

    print(f"\nActual hit counts per sector / property "
          f"(P(rand ≥ actual) = fraction of random subsets with hits ≥ actual; "
          f"low = unusually enriched):")
    print(f"{'property':<28} ", end="")
    for s in SECTORS:
        print(f" {s+'(n=' + str(len(sector_primes[s])) + ')':<14}", end="")
    print()
    print("─" * 28 + " " + "─" * 14 * 3)

    for name, pred in props:
        row = f"{name:<28} "
        for s in SECTORS:
            actual_hits = sum(1 for p in sector_primes[s] if pred(p))
            n = len(sector_primes[s])
            # null
            pool = window_pools[s]
            if n > len(pool):
                row += f" {'n>pool':<14}"
                continue
            ge_count = 0
            for _ in range(N_NULL_TRIALS_T3):
                sample = rng.sample(pool, n)
                hits = sum(1 for p in sample if pred(p))
                if hits >= actual_hits:
                    ge_count += 1
            pctile = ge_count / N_NULL_TRIALS_T3 * 100
            row += f" {actual_hits}/{n} ({pctile:5.1f}%)  "
        print(row)


# ─────────────────────────────────────────────────────────────────────
# T4 — Gap structure
# ─────────────────────────────────────────────────────────────────────

def t4_gap_structure(cells, t1_results, rng):
    print("\n" + "=" * 90)
    print("T4 — GAP STRUCTURE WITHIN SECTOR")
    print("=" * 90)
    print("\nNote: endpoints (sector min and max prime) are fixed in both observed and")
    print("random subsets, so mean gap is identical by construction = (hi-lo)/(n-1).")
    print("Only gap variance / regularity is informative.")

    sector_primes = {s: set() for s in SECTORS}
    for c in cells:
        sector_primes[c["sector"]].add(c["p"])

    for sector in SECTORS:
        primes = sorted(sector_primes[sector])
        n = len(primes)
        if n < 2:
            print(f"\n{sector}: n={n}, insufficient primes for gaps.")
            continue
        gaps = [primes[i+1] - primes[i] for i in range(n-1)]
        mean_g = sum(gaps) / len(gaps)
        if len(gaps) > 1:
            var_g = sum((g - mean_g)**2 for g in gaps) / len(gaps)
        else:
            var_g = 0.0

        print(f"\n{sector} sector: {n} primes, gaps = {gaps}")
        print(f"  mean gap = {mean_g:.2f} (= (hi-lo)/(n-1)), variance = {var_g:.2f}")

        if n < 3:
            print(f"  (n<3, single gap reported descriptively only)")
            continue

        # null: random subsets of same cardinality with FIXED endpoints
        pool = t1_results[sector]["available_primes"]
        lo, hi = primes[0], primes[-1]
        # internal pool: primes strictly between lo and hi
        internal_pool = [p for p in pool if lo < p < hi]
        if (n - 2) > len(internal_pool):
            print(f"  Cannot run null: cardinality - 2 > internal pool size")
            continue

        null_vars = []
        for _ in range(N_NULL_TRIALS_T4):
            internal = rng.sample(internal_pool, n - 2)
            sample = sorted([lo, hi] + internal)
            sg = [sample[i+1] - sample[i] for i in range(n-1)]
            mg = sum(sg) / len(sg)  # equal to mean_g by construction
            if len(sg) > 1:
                null_vars.append(sum((g - mg)**2 for g in sg) / len(sg))
            else:
                null_vars.append(0.0)

        n_le_var = sum(1 for v in null_vars if v <= var_g)
        n_ge_var = sum(1 for v in null_vars if v >= var_g)
        print(f"  Variance: actual = {var_g:.2f}, "
              f"null mean = {sum(null_vars)/len(null_vars):.2f}")
        print(f"  P(rand var ≤ actual): {n_le_var/N_NULL_TRIALS_T4*100:.1f}% "
              f"(low = unusually regular spacing)")
        print(f"  P(rand var ≥ actual): {n_ge_var/N_NULL_TRIALS_T4*100:.1f}% "
              f"(low = unusually irregular spacing)")


# ─────────────────────────────────────────────────────────────────────
# T5 — Local candidate-density diagnostic
# ─────────────────────────────────────────────────────────────────────

def t5_density(cells, candidates, rng):
    print("\n" + "=" * 90)
    print("T5 — LOCAL CANDIDATE-DENSITY DIAGNOSTIC")
    print("=" * 90)

    cand_vals_sorted = sorted(c[0] for c in candidates)

    def fast_best_match_residual(target):
        """Return residual % of nearest candidate, using bisect."""
        idx = bisect.bisect_left(cand_vals_sorted, target)
        # check idx and idx-1
        best = float("inf")
        if idx < len(cand_vals_sorted):
            best = abs(target - cand_vals_sorted[idx]) / target * 100
        if idx > 0:
            r = abs(target - cand_vals_sorted[idx-1]) / target * 100
            if r < best:
                best = r
        return best

    print(f"\n{'cell':<26} {'target':>10} {'best':>8} {'res%':>7} "
          f"{'cand-':>7} {'cand+':>7} {'spacing':>9} {'rand pctile':>11}")
    print("─" * 26 + " " + "─" * 10 + " " + "─" * 8 + " " + "─" * 7 + " "
          + "─" * 7 + " " + "─" * 7 + " " + "─" * 9 + " " + "─" * 11)

    summary_by_sector = {s: [] for s in SECTORS}

    for c in cells:
        target = c["target"]
        val = c["val"]
        idx = cand_vals_sorted.index(val)
        below = cand_vals_sorted[idx-1] if idx > 0 else None
        above = cand_vals_sorted[idx+1] if idx + 1 < len(cand_vals_sorted) else None

        gap_below = (val - below) / val * 100 if below is not None else None
        gap_above = (above - val) / val * 100 if above is not None else None

        spacings = []
        if gap_below is not None:
            spacings.append(gap_below)
        if gap_above is not None:
            spacings.append(gap_above)
        local_spacing = sum(spacings) / len(spacings) if spacings else 0.0

        # null: random ratios in ±20% window around target
        lo = target * 0.8
        hi = target * 1.2
        ge_count = 0
        for _ in range(N_NULL_TRIALS_T5):
            r = rng.uniform(lo, hi)
            rres = fast_best_match_residual(r)
            if rres >= c["residual"]:
                ge_count += 1
        rand_pctile = ge_count / N_NULL_TRIALS_T5 * 100  # higher = real fit better than typical

        gb_str = f"{gap_below:.1f}%" if gap_below is not None else "—"
        ga_str = f"{gap_above:.1f}%" if gap_above is not None else "—"
        cell_label = f"{c['ratio']}@{c['scale']}"
        print(f"{cell_label:<26} {target:>10.3f} {val:>8} {c['residual']:>6.3f}% "
              f"{gb_str:>7} {ga_str:>7} {local_spacing:>8.2f}% {rand_pctile:>10.1f}%")

        summary_by_sector[c["sector"]].append({
            "residual": c["residual"],
            "spacing": local_spacing,
            "rand_pctile": rand_pctile,
        })

    # per-sector summary
    print()
    for sector in SECTORS:
        rows = summary_by_sector[sector]
        if not rows:
            continue
        mean_res = sum(r["residual"] for r in rows) / len(rows)
        mean_sp = sum(r["spacing"] for r in rows) / len(rows)
        mean_pct = sum(r["rand_pctile"] for r in rows) / len(rows)
        print(f"  {sector}: mean residual {mean_res:.3f}%, "
              f"mean local spacing {mean_sp:.2f}%, "
              f"mean rand-pctile {mean_pct:.1f}%")

    print()
    print("Reading: cand-, cand+ are the grammar candidates immediately below/above")
    print("the best-fit candidate. rand-pctile = fraction of random ratios (in ±20%")
    print("window of target) with worse residual than the actual fit. High values")
    print("= real fit better than typical, low values = fit no better than density.")


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    t0_start = time.time()
    rng = random.Random(SEED)

    print("=" * 90)
    print("FERMION MASS PRIMES TEST 15")
    print("Sector prime-set examination")
    print("=" * 90)
    print(f"Common scales: {list(ANTUSCH_SCALES.keys())}")
    print(f"Sectors: {SECTORS}")
    print(f"Random trials: T2={N_NULL_TRIALS_T2}, T3={N_NULL_TRIALS_T3}, "
          f"T4={N_NULL_TRIALS_T4}, T5={N_NULL_TRIALS_T5}")
    print(f"Seed: {SEED}")

    cells, candidates = build_cells()
    print(f"\nCells built: {len(cells)} (= {len(RATIO_NAMES)} ratios "
          f"× {len(ANTUSCH_SCALES)} common scales)")
    print(f"Candidates: {len(candidates)}")

    t0_occupancy(cells)
    t1_results = t1_skip_pattern(cells)
    t2_disjointness_null(cells, t1_results, rng)
    t3_algebraic_scan(cells, t1_results, rng)
    t4_gap_structure(cells, t1_results, rng)
    t5_density(cells, candidates, rng)

    print(f"\nTotal runtime: {time.time() - t0_start:.1f}s")
    print("Done.")
