"""
fermion_mass_primes_test17.py

Bridge substitution null (local enumeration).

Tests whether the cross-family bridges (t/tau, t/b) used in the chain
hypothesis carry structural information beyond their numerical magnitude.

For each common scale and each chain:
  1. Hold within-family components fixed at their v1.0 grammar fits.
  2. Compute B_required = (m_t/m_X)_observed / product(within-family fits).
  3. Enumerate every grammar candidate 2^a x p in [B_obs/2, 2 * B_obs].
  4. Substitute each into the chain; rank by absolute endpoint residual.
  5. Report rank/percentile of B_fit (the test14-selected bridge) in the pool.

Reading rules (pre-committed):
  Strong:  percentile <= 5
  Modest:  5  < percentile <= 25
  Weak:    25 < percentile <= 50
  None:    percentile >  50

This is a local bridge-density null. A strong rank supports the bridge
edge being more than magnitude-compatible. It does not establish modal
equivalence on its own; that requires the closure constraint of test19.

See scripts/README.md for full documentation.
"""

import time

import fermion_mass_primes_test11 as t11
import fermion_mass_primes_test14 as t14


# Tolerance for treating two residuals as identical (tie detection)
TIE_TOLERANCE = 1e-15

# Magnitude window factor (pre-committed)
WINDOW_FACTOR = 2.0


# ── Chain definitions (mirrors test16, but bridge-only) ──
# Each chain is a list of (edge_label, source, fitter_key) for within-family
# components, plus a bridge_key identifying the cross-family bridge.
CHAINS = {
    "lepton (m_t/m_e)": {
        "within_family": [
            ("tau/mu", "test11", "mtau/mmu"),
            ("mu/e",   "test11", "mmu/me"),
        ],
        "bridge_key": "t/tau",
        "endpoint_heavy": "t",
        "endpoint_light": "e",
    },
    "down (m_t/m_d)": {
        "within_family": [
            ("b/s", "test11", "mb/ms"),
            ("s/d", "test11", "ms/md"),
        ],
        "bridge_key": "t/b",
        "endpoint_heavy": "t",
        "endpoint_light": "d",
    },
}

COMMON_SCALES = list(t11.ANTUSCH_SCALES.keys())


def get_within_family_fits(masses, candidates):
    """Return dict {ratio_name: (R_observed, val_fit, a, p)}."""
    ratios = t11.compute_6_ratios_from_yukawa(masses)
    out = {}
    for name, R in zip(t11.RATIO_NAMES, ratios):
        val, a, p, _abs_res = t14.find_best_match(R, candidates)
        out[name] = (R, val, a, p)
    return out


def get_bridge_fit(masses_by_particle, bridge_key, candidates):
    """Return (R_observed, val_fit, a, p) for a cross-family bridge.

    bridge_key is one of the labels in t14.CROSS_FAMILY_PAIRS (e.g. "t/tau").
    """
    for label, k1, k2 in t14.CROSS_FAMILY_PAIRS:
        if label == bridge_key:
            m1, m2 = masses_by_particle[k1], masses_by_particle[k2]
            R = max(m1, m2) / min(m1, m2)
            val, a, p, _abs_res = t14.find_best_match(R, candidates)
            return R, val, a, p
    raise KeyError(f"bridge_key not in CROSS_FAMILY_PAIRS: {bridge_key}")


def enumerate_pool(B_obs, candidates, factor=WINDOW_FACTOR):
    """Enumerate all candidates 2^a x p in [B_obs/factor, factor * B_obs].

    Returns sorted list of (val, a, p) tuples.
    """
    lo = B_obs / factor
    hi = B_obs * factor
    return [(val, a, p) for (val, a, p) in candidates if lo <= val <= hi]


def assign_ranks(residuals, tol=TIE_TOLERANCE):
    """Assign ranks (1-indexed, lowest residual = rank 1) with tie averaging.

    Returns list of ranks parallel to the input residuals list.
    """
    n = len(residuals)
    indexed = sorted(enumerate(residuals), key=lambda x: x[1])
    ranks = [0.0] * n

    i = 0
    while i < n:
        j = i
        # Group ties
        while j + 1 < n and abs(indexed[j + 1][1] - indexed[i][1]) < tol:
            j += 1
        # Average rank of positions i..j (1-indexed)
        avg_rank = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            orig_idx = indexed[k][0]
            ranks[orig_idx] = avg_rank
        i = j + 1
    return ranks


def percentile(rank, n):
    """percentile = 100 * (rank - 1) / (n - 1).

    rank is 1-indexed (best = 1). For n == 1, returns 0.
    """
    if n <= 1:
        return 0.0
    return 100.0 * (rank - 1) / (n - 1)


def evidence_band(pct):
    """Pre-committed reading rules."""
    if pct <= 5:
        return "strong"
    if pct <= 25:
        return "modest"
    if pct <= 50:
        return "weak"
    return "none"


def format_match(a, p):
    if a == 0:
        return f"{p}"
    if a == 1:
        return f"2x{p}"
    if a == 2:
        return f"2^2x{p}"
    return f"2^3x{p}"


if __name__ == "__main__":
    t0 = time.time()

    print("=" * 100)
    print("FERMION MASS PRIMES TEST 17")
    print("Bridge substitution null (local enumeration)")
    print("=" * 100)
    print("Local null check: is the test14-selected bridge unusually well-aligned")
    print("with B_required among grammar candidates of similar magnitude?")
    print()
    print("Reading rules (pre-committed):")
    print("  Strong:  percentile <= 5")
    print("  Modest:  5  < percentile <= 25")
    print("  Weak:    25 < percentile <= 50")
    print("  None:    percentile >  50")
    print()
    print("Chains tested:")
    for chain_name, chain_def in CHAINS.items():
        wf = " x ".join(s[0] for s in chain_def["within_family"])
        print(f"  {chain_name}: bridge {chain_def['bridge_key']} x {wf}")
    print()
    print(f"Scales: {len(COMMON_SCALES)} common scales (PDG excluded by design)")
    print(f"Window factor: {WINDOW_FACTOR} (pool = candidates in [B_obs/{WINDOW_FACTOR}, "
          f"{WINDOW_FACTOR} * B_obs])")

    # ── Generate candidates (same upper bound as test14) ──
    all_R = []
    for scale_label in COMMON_SCALES:
        raw = t11.ANTUSCH_SCALES[scale_label]
        within = t11.compute_6_ratios_from_yukawa(raw)
        masses = {
            "u": raw["yu"], "c": raw["yc"], "t": raw["yt"],
            "d": raw["yd"], "s": raw["ys"], "b": raw["yb"],
            "e": raw["ye"], "mu": raw["ymu"], "tau": raw["ytau"],
        }
        cross = []
        for _label, k1, k2 in t14.CROSS_FAMILY_PAIRS:
            cross.append(max(masses[k1], masses[k2]) / min(masses[k1], masses[k2]))
        all_R.extend(within)
        all_R.extend(cross)
    max_R = max(all_R)
    candidates = t14.generate_candidates(max_R * 2)
    print(f"Candidates generated: {len(candidates)} (2^a x p up to {max_R*2:.0f})")

    # ── Per-scale, per-chain analysis ──
    # results[scale_label][chain_name] = dict with all per-cell data
    results = {}
    aggregate_counts = {chain: {"strong": 0, "modest": 0, "weak": 0, "none": 0}
                        for chain in CHAINS.keys()}

    for scale_label in COMMON_SCALES:
        raw = t11.ANTUSCH_SCALES[scale_label]
        masses_by_particle = {
            "u": raw["yu"], "c": raw["yc"], "t": raw["yt"],
            "d": raw["yd"], "s": raw["ys"], "b": raw["yb"],
            "e": raw["ye"], "mu": raw["ymu"], "tau": raw["ytau"],
        }
        within_fits = get_within_family_fits(raw, candidates)

        scale_results = {}
        for chain_name, chain_def in CHAINS.items():
            # Endpoint observed
            heavy = chain_def["endpoint_heavy"]
            light = chain_def["endpoint_light"]
            endpoint_obs = masses_by_particle[heavy] / masses_by_particle[light]

            # Product of within-family fitted values
            wf_prod_fit = 1.0
            wf_steps = []
            for edge_label, _src, key in chain_def["within_family"]:
                R_obs, val_fit, a, p = within_fits[key]
                wf_prod_fit *= val_fit
                wf_steps.append({
                    "label": edge_label,
                    "R_obs": R_obs,
                    "val_fit": val_fit,
                    "a": a,
                    "p": p,
                })

            # B_required and B_obs
            B_required = endpoint_obs / wf_prod_fit
            B_obs, B_fit_val, B_fit_a, B_fit_p = get_bridge_fit(
                masses_by_particle, chain_def["bridge_key"], candidates
            )

            # Enumerate the pool
            pool = enumerate_pool(B_obs, candidates, factor=WINDOW_FACTOR)
            pool_size = len(pool)

            # Pool membership assertion
            B_fit_tuple = (B_fit_val, B_fit_a, B_fit_p)
            if B_fit_tuple not in pool:
                print()
                print("=" * 100)
                print(f"HARD DIAGNOSTIC — POOL MEMBERSHIP FAILED")
                print("=" * 100)
                print(f"Scale:      {scale_label}")
                print(f"Chain:      {chain_name}")
                print(f"Bridge key: {chain_def['bridge_key']}")
                print(f"B_obs:      {B_obs:.6f}")
                print(f"Window:     [{B_obs/WINDOW_FACTOR:.6f}, "
                      f"{B_obs*WINDOW_FACTOR:.6f}]")
                print(f"B_fit:      {B_fit_val} (a={B_fit_a}, p={B_fit_p})")
                print(f"Pool size:  {pool_size}")
                print()
                print("B_fit was selected by test14 as nearest to B_obs, but lies")
                print("outside the [B_obs/2, 2*B_obs] window. The null definition")
                print("has failed for this cell.")
                print()
                print("Halting test.")
                raise SystemExit(1)

            # Score every candidate by chain-endpoint residual against observed
            scored = []
            for (val, a, p) in pool:
                pred = wf_prod_fit * val
                resid = abs(pred - endpoint_obs) / endpoint_obs
                scored.append((val, a, p, resid))

            # Assign ranks (1-indexed, with tie averaging)
            residuals = [s[3] for s in scored]
            ranks = assign_ranks(residuals, tol=TIE_TOLERANCE)

            # Find B_fit's rank and percentile
            B_fit_idx = pool.index(B_fit_tuple)
            B_fit_rank = ranks[B_fit_idx]
            B_fit_residual = scored[B_fit_idx][3]
            B_fit_pct = percentile(B_fit_rank, pool_size)
            band = evidence_band(B_fit_pct)
            aggregate_counts[chain_name][band] += 1

            # Best candidate (rank closest to 1)
            best_idx = min(range(pool_size), key=lambda i: residuals[i])
            best = scored[best_idx]

            # Detect ties involving B_fit
            tied_at_bfit = sum(1 for r in ranks if abs(r - B_fit_rank) < 1e-9)

            scale_results[chain_name] = {
                "endpoint_obs": endpoint_obs,
                "wf_steps": wf_steps,
                "wf_prod_fit": wf_prod_fit,
                "B_obs": B_obs,
                "B_required": B_required,
                "B_fit": B_fit_val,
                "B_fit_a": B_fit_a,
                "B_fit_p": B_fit_p,
                "B_fit_residual": B_fit_residual,
                "pool_size": pool_size,
                "rank": B_fit_rank,
                "percentile": B_fit_pct,
                "band": band,
                "best_val": best[0],
                "best_a": best[1],
                "best_p": best[2],
                "best_residual": best[3],
                "tied_at_bfit": tied_at_bfit,
                "scored": scored,
                "ranks": ranks,
            }

        results[scale_label] = scale_results

    # ── Per-scale, per-chain output ──
    for scale_label in COMMON_SCALES:
        print(f"\n{'━' * 100}")
        print(f"SCALE: {scale_label}")
        print(f"{'━' * 100}")

        for chain_name in CHAINS.keys():
            r = results[scale_label][chain_name]

            print(f"\n  CHAIN: {chain_name}")
            print(f"  {'─' * 96}")

            # Within-family steps (for context)
            print(f"    Within-family fits held fixed:")
            for step in r["wf_steps"]:
                form = format_match(step["a"], step["p"])
                print(f"      {step['label']:<8} R_obs={step['R_obs']:>12.4f}  "
                      f"fit={step['val_fit']:<5} {form}")
            print(f"    Within-family product: {r['wf_prod_fit']:.6e}")
            print()
            print(f"    Endpoint observed (m_t/m_X): {r['endpoint_obs']:.6e}")
            print(f"    B_obs       (observed bridge ratio):  {r['B_obs']:.6f}")
            print(f"    B_required  (chain-target bridge):    {r['B_required']:.6f}")
            print(f"    B_obs - B_required: {(r['B_obs'] - r['B_required']):+.6f} "
                  f"({(r['B_obs'] - r['B_required'])/r['B_required']*100:+.4f}%)")
            print()
            B_fit_form = format_match(r["B_fit_a"], r["B_fit_p"])
            best_form = format_match(r["best_a"], r["best_p"])
            print(f"    B_fit (test14 selection): {r['B_fit']} = {B_fit_form}")
            print(f"      residual against endpoint: {r['B_fit_residual']*100:.4f}%")
            print(f"    Pool size: {r['pool_size']} "
                  f"(window = [{r['B_obs']/WINDOW_FACTOR:.2f}, "
                  f"{r['B_obs']*WINDOW_FACTOR:.2f}])")
            print(f"    Best candidate in pool: {r['best_val']} = {best_form}")
            print(f"      residual against endpoint: {r['best_residual']*100:.4f}%")
            print()
            tie_note = ""
            if r["tied_at_bfit"] > 1:
                tie_note = f"  (tied with {r['tied_at_bfit'] - 1} other candidate(s))"
            print(f"    B_fit rank:        {r['rank']:.1f} of {r['pool_size']}{tie_note}")
            print(f"    B_fit percentile:  {r['percentile']:.2f}")
            print(f"    Evidence band:     {r['band']}")

    # ── Grand summary ──
    print(f"\n{'━' * 110}")
    print("GRAND SUMMARY — B_FIT RANK / PERCENTILE / BAND PER CHAIN PER SCALE")
    print(f"{'━' * 110}")
    header = f"\n{'Scale':<12}"
    for chain_name in CHAINS.keys():
        header += f" {chain_name:>40}"
    print(header)
    sub = f"{'':<12}"
    for _ in CHAINS.keys():
        sub += f"     {'rank/N':>10} {'pct':>8} {'band':>8}"
        sub += " " * (40 - 30)
    print(sub)
    sep = f"{'─' * 12}"
    for _ in CHAINS:
        sep += f" {'─' * 40}"
    print(sep)

    for scale_label in COMMON_SCALES:
        row = f"{scale_label:<12}"
        for chain_name in CHAINS.keys():
            r = results[scale_label][chain_name]
            cell = (f"     {r['rank']:>4.1f}/{r['pool_size']:<4} "
                    f"{r['percentile']:>7.2f} {r['band']:>8}")
            cell = cell + " " * (40 - len(cell) + 5)
            row += cell[:41]
        print(row)

    # ── Aggregate evidence summary ──
    print(f"\n{'━' * 100}")
    print("AGGREGATE EVIDENCE — N OUT OF 9 SCALES")
    print(f"{'━' * 100}")
    print()
    n_scales = len(COMMON_SCALES)
    for chain_name in CHAINS.keys():
        c = aggregate_counts[chain_name]
        print(f"  {chain_name}")
        print(f"    Strong:  {c['strong']}/{n_scales}")
        print(f"    Modest:  {c['modest']}/{n_scales}")
        print(f"    Weak:    {c['weak']}/{n_scales}")
        print(f"    None:    {c['none']}/{n_scales}")
        print()

    # ── Caution restated ──
    print(f"{'━' * 100}")
    print("INTERPRETATION CAUTION")
    print(f"{'━' * 100}")
    print()
    print("Because B_fit was selected by test14 as nearest to B_obs, and B_required")
    print("is usually close to B_obs, a strong rank partly reflects nearest-candidate")
    print("behaviour. The test is strongest at scales where B_required is noticeably")
    print("shifted from B_obs and B_fit still ranks well. Per-scale reporting above")
    print("includes both B_obs and B_required so the shift is visible.")
    print()
    print("This is a local bridge-density null. A strong rank supports the bridge")
    print("edge being more than magnitude-compatible. Modal-equivalence requires the")
    print("closure constraint of test19.")

    total_elapsed = time.time() - t0
    print(f"\nTotal runtime: {total_elapsed:.1f}s")
    print("Done.")
