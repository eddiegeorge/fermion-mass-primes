"""
fermion_mass_primes_test18.py

Hierarchy-matched 9 charged-mass null.

Tests whether the joint structure of v1.0's grammar fits — within-family
sub-percent fits, gen-3 cross-family bridges, and chain composition
residuals — occurs more often in real charged-fermion data than in random
hierarchies of comparable shape.

Per-scale only: tests whether the real per-scale co-occurrence pattern
is unusual under the hierarchy-matched charged-mass null. Does NOT test
cross-scale stability.

Random-system generator (Dirichlet log-gaps):
  1. Use real charged-fermion rank order at each scale
  2. Draw S = log(m_max/m_min) uniformly in [S_real - log 2, S_real + log 2]
  3. Draw 8 ordered log gaps from Dirichlet(1,...,1) on the 8-simplex
  4. Construct 9 ordered log masses by cumulative sum
  5. Assign labels by real rank order

Scoring features (pre-committed):
  1. Within-family sub-percent fit count (6 ratios at <= 1%)
  2. Selected cross-family fit count: t/tau, t/b, c/s at <= 3%
  3. Absolute lepton chain residual (m_t/m_e under test16 construction)
  4. Absolute down chain residual (m_t/m_d under test16 construction)

Joint score: random matches or exceeds if all four feature comparisons
satisfy the inequality.

Reading rules (pre-committed, per scale):
  Strong:  < 0.1%
  Modest:  < 1%
  Weak:    < 5%
  None:    >= 5%

Aggregate as "strong at N/9, modest at M/9, weak at K/9, none at J/9".
Not combined into a single p-value.

See scripts/README.md for full documentation.
"""

import bisect
import time

import numpy as np

import fermion_mass_primes_test11 as t11
import fermion_mass_primes_test14 as t14


# Trial count (per scale) — pre-committed
N_TRIALS = 1_000_000

# Magnitude window for log-span draw — pre-committed
LOG_SPAN_HALFWIDTH = np.log(2.0)

# Scoring thresholds — pre-committed
WITHIN_FAMILY_RESID_THRESHOLD = 1.0   # percent
CROSS_FAMILY_RESID_THRESHOLD = 3.0    # percent

# Grammar floor — same as test14
GRAMMAR_FLOOR = 3.0

# Selected cross-family ratios (t/tau, t/b, c/s)
SELECTED_CROSS = ["t/tau", "t/b", "c/s"]

COMMON_SCALES = list(t11.ANTUSCH_SCALES.keys())

# Charged-fermion labels in some canonical order
CHARGED_LABELS = ["u", "d", "s", "c", "b", "e", "mu", "tau", "t"]

# Within-family adjacent-generation pairs (heavier_label, lighter_label)
WITHIN_FAMILY_PAIRS = [
    ("c",   "u"),    # mc/mu
    ("t",   "c"),    # mt/mc
    ("s",   "d"),    # ms/md
    ("b",   "s"),    # mb/ms
    ("mu",  "e"),    # mmu/me
    ("tau", "mu"),   # mtau/mmu
]


def get_charged_masses(yukawa_dict):
    """Pull charged-fermion masses (Yukawa-equivalent) keyed by label."""
    return {
        "u": yukawa_dict["yu"],
        "c": yukawa_dict["yc"],
        "t": yukawa_dict["yt"],
        "d": yukawa_dict["yd"],
        "s": yukawa_dict["ys"],
        "b": yukawa_dict["yb"],
        "e": yukawa_dict["ye"],
        "mu": yukawa_dict["ymu"],
        "tau": yukawa_dict["ytau"],
    }


def real_rank_order(masses):
    """Return the 9 labels sorted by descending mass (heaviest first)."""
    return sorted(CHARGED_LABELS, key=lambda lab: -masses[lab])


def nearest_value_and_residual(target, candidate_vals):
    """Return (nearest_value, residual_signed_pct).

    candidate_vals must be sorted ascending. Uses bisect for O(log N) lookup.
    residual_signed_pct = (val - target) / target * 100.

    Below-floor protection: if target < grammar floor (3), returns
    (None, inf) since no grammar candidate can validly fit a ratio below 3.
    """
    if target < GRAMMAR_FLOOR:
        return None, float("inf")
    n = len(candidate_vals)
    idx = bisect.bisect_left(candidate_vals, target)
    if idx == 0:
        c = candidate_vals[0]
    elif idx == n:
        c = candidate_vals[n - 1]
    else:
        c_left = candidate_vals[idx - 1]
        c_right = candidate_vals[idx]
        c = c_left if abs(c_left - target) <= abs(c_right - target) else c_right
    return c, (c - target) / target * 100.0


def score_system(masses, candidate_vals):
    """Compute the four scoring features for a charged-mass system.

    Returns:
        within_count: number of within-family adjacent ratios at <= 1%
        cross_count:  number of selected cross-family ratios at <= 3%
        lepton_resid_abs:  |composed residual of m_t/m_e| as fraction
                           (inf if any chain component is below grammar floor)
        down_resid_abs:    |composed residual of m_t/m_d| as fraction
                           (inf if any chain component is below grammar floor)
    """
    # Within-family fits and signed residuals (for chain composition)
    within_signed = {}  # {(heavier, lighter): signed_residual_fraction or None}
    within_count = 0
    for heavier, lighter in WITHIN_FAMILY_PAIRS:
        R = masses[heavier] / masses[lighter]
        _c, signed_pct = nearest_value_and_residual(R, candidate_vals)
        # signed_pct is inf for below-floor; mark the component as invalid
        if signed_pct == float("inf"):
            within_signed[(heavier, lighter)] = None
        else:
            within_signed[(heavier, lighter)] = signed_pct / 100.0
            if abs(signed_pct) <= WITHIN_FAMILY_RESID_THRESHOLD:
                within_count += 1

    # Selected cross-family fits
    cross_signed = {}  # {pair_label: signed_residual_fraction or None}
    cross_count = 0
    for label, k1, k2 in t14.CROSS_FAMILY_PAIRS:
        if label not in SELECTED_CROSS:
            continue
        m1, m2 = masses[k1], masses[k2]
        R = max(m1, m2) / min(m1, m2)
        _c, signed_pct = nearest_value_and_residual(R, candidate_vals)
        if signed_pct == float("inf"):
            cross_signed[label] = None
        else:
            cross_signed[label] = signed_pct / 100.0
            if abs(signed_pct) <= CROSS_FAMILY_RESID_THRESHOLD:
                cross_count += 1

    # Chain composition residuals (test16 construction)
    # If any component is below the grammar floor (None), the chain is
    # not well-defined under v1.0's grammar; report inf residual.
    # Lepton: m_t/m_e via t/tau x tau/mu x mu/e
    delta_t_tau = cross_signed.get("t/tau")
    delta_tau_mu = within_signed[("tau", "mu")]
    delta_mu_e = within_signed[("mu", "e")]
    if None in (delta_t_tau, delta_tau_mu, delta_mu_e):
        lepton_resid_abs = float("inf")
    else:
        lepton_resid = (
            (1 + delta_t_tau) * (1 + delta_tau_mu) * (1 + delta_mu_e) - 1
        )
        lepton_resid_abs = abs(lepton_resid)

    # Down: m_t/m_d via t/b x b/s x s/d
    delta_t_b = cross_signed.get("t/b")
    delta_b_s = within_signed[("b", "s")]
    delta_s_d = within_signed[("s", "d")]
    if None in (delta_t_b, delta_b_s, delta_s_d):
        down_resid_abs = float("inf")
    else:
        down_resid = (
            (1 + delta_t_b) * (1 + delta_b_s) * (1 + delta_s_d) - 1
        )
        down_resid_abs = abs(down_resid)

    return (
        within_count,
        cross_count,
        lepton_resid_abs,
        down_resid_abs,
    )


def generate_random_system(rank_order, S_real, m_min, rng):
    """Generate one random charged-mass system per the test18 generator.

    Parameters:
        rank_order: list of 9 labels, heaviest first (real rank order)
        S_real:     real-data log span at this scale (log(m_max/m_min))
        m_min:      real-data minimum mass at this scale (anchors absolute
                    scale; only ratios matter for scoring)
        rng:        numpy RandomState for sampling

    Returns:
        masses dict {label: mass_value}
    """
    # Draw total log span S in [S_real - log2, S_real + log2]
    S = rng.uniform(S_real - LOG_SPAN_HALFWIDTH, S_real + LOG_SPAN_HALFWIDTH)

    # Draw 8 log gaps from Dirichlet(1,...,1), scale to sum S
    gaps = rng.dirichlet(np.ones(8)) * S

    # Build 9 ordered log masses from m_min, ascending
    log_m_min = np.log(m_min)
    log_masses_ascending = np.concatenate(
        ([log_m_min], log_m_min + np.cumsum(gaps))
    )
    masses_ascending = np.exp(log_masses_ascending)

    # rank_order is heaviest first; reverse to lightest first to match ascending masses
    rank_lightest_first = rank_order[::-1]

    return {
        label: float(mass)
        for label, mass in zip(rank_lightest_first, masses_ascending)
    }


def evidence_band(fraction):
    """Pre-committed reading rules."""
    if fraction < 0.001:
        return "strong"
    if fraction < 0.01:
        return "modest"
    if fraction < 0.05:
        return "weak"
    return "none"


def format_fraction(matches, n_trials):
    """Format match-or-exceed count as a readable fraction."""
    if matches == 0:
        return f"< {1/n_trials:.0e}  (0/{n_trials:,})"
    return f"{matches/n_trials:.4%}  ({matches:,}/{n_trials:,})"


if __name__ == "__main__":
    t0 = time.time()

    print("=" * 100)
    print("FERMION MASS PRIMES TEST 18")
    print("Hierarchy-matched 9 charged-mass null")
    print("=" * 100)
    print("Per-scale joint-structure null. Tests whether the real per-scale")
    print("co-occurrence of (within-family sub-percent fits, cross-family fits,")
    print("clean chain residuals) is unusual under random plausible hierarchies.")
    print()
    print("Reading rules (pre-committed, per scale):")
    print("  Strong:  fraction matching or exceeding < 0.1%")
    print("  Modest:  < 1%")
    print("  Weak:    < 5%")
    print("  None:    >= 5%")
    print()
    print(f"Trials per scale: {N_TRIALS:,}")
    print(f"Log-span window: real +/- log(2)")
    print(f"Within-family threshold: {WITHIN_FAMILY_RESID_THRESHOLD}%")
    print(f"Cross-family threshold:  {CROSS_FAMILY_RESID_THRESHOLD}%")
    print()
    print(f"Selected cross-family ratios: {SELECTED_CROSS}")
    print(f"Scales: {len(COMMON_SCALES)} common scales (PDG excluded by design)")

    # ── Generate candidate set ──
    # Random systems can produce scored ratios up to exp(S_real + log 2),
    # which can reach ~700,000 at the largest common-scale span. Generate
    # candidates wide enough to cover that range so random systems are
    # not artificially punished by candidate truncation.
    max_possible_R = 0.0
    for scale_label in COMMON_SCALES:
        raw = t11.ANTUSCH_SCALES[scale_label]
        masses = get_charged_masses(raw)
        m_min = min(masses.values())
        m_max = max(masses.values())
        S_real = np.log(m_max / m_min)
        # Random generator allows S up to S_real + log 2
        max_possible_R = max(max_possible_R, np.exp(S_real + LOG_SPAN_HALFWIDTH))
    candidates = t14.generate_candidates(max_possible_R * 1.01)
    candidate_vals = sorted(c[0] for c in candidates)
    print(f"Candidates generated: {len(candidate_vals)} "
          f"(2^a x p up to {max_possible_R*1.01:.0f})")
    print()

    # ── Per-scale analysis ──
    results = {}
    aggregate = {"strong": 0, "modest": 0, "weak": 0, "none": 0}

    # Single RNG, advancing across all scales — random sequences are
    # genuinely independent across scales rather than identical per scale.
    rng = np.random.RandomState(seed=42)

    for scale_label in COMMON_SCALES:
        print(f"{'━' * 100}")
        print(f"SCALE: {scale_label}")
        print(f"{'━' * 100}")

        raw = t11.ANTUSCH_SCALES[scale_label]
        masses_real = get_charged_masses(raw)
        rank = real_rank_order(masses_real)
        m_min = min(masses_real.values())
        m_max = max(masses_real.values())
        S_real = np.log(m_max / m_min)

        print(f"  Real rank order (heaviest first): "
              f"{' > '.join(rank)}")
        print(f"  Real log span S_real = {S_real:.4f}")

        # Score real data
        real_scores = score_system(masses_real, candidate_vals)
        real_within, real_cross, real_lepton, real_down = real_scores
        print(f"  Real scores:")
        print(f"    within-family sub-percent fits "
              f"(<= {WITHIN_FAMILY_RESID_THRESHOLD}%):  {real_within}/6")
        print(f"    selected cross-family fits "
              f"(<= {CROSS_FAMILY_RESID_THRESHOLD}%):     {real_cross}/3")
        print(f"    |lepton chain residual|:  {real_lepton*100:.4f}%")
        print(f"    |down   chain residual|:  {real_down*100:.4f}%")

        # Run random systems
        match_count = 0
        # Histograms
        within_hist = np.zeros(7, dtype=int)   # counts 0..6
        cross_hist = np.zeros(4, dtype=int)    # counts 0..3
        # Track tail of chain residuals for inspection
        lepton_resids = []
        down_resids = []
        n_resid_sample = min(N_TRIALS, 10000)  # subsample for histogram
        sample_stride = max(1, N_TRIALS // n_resid_sample)

        print(f"  Running {N_TRIALS:,} random systems...", flush=True)
        loop_start = time.time()

        for trial in range(N_TRIALS):
            sys_masses = generate_random_system(rank, S_real, m_min, rng)
            scores = score_system(sys_masses, candidate_vals)
            sys_within, sys_cross, sys_lepton, sys_down = scores

            # Histograms
            within_hist[sys_within] += 1
            cross_hist[sys_cross] += 1

            # Subsample chain residuals for histogram
            if trial % sample_stride == 0:
                lepton_resids.append(sys_lepton)
                down_resids.append(sys_down)

            # Joint match/exceed test
            if (sys_within >= real_within
                    and sys_cross >= real_cross
                    and sys_lepton <= real_lepton
                    and sys_down <= real_down):
                match_count += 1

        elapsed = time.time() - loop_start
        fraction = match_count / N_TRIALS
        band = evidence_band(fraction)
        aggregate[band] += 1

        print(f"  ...completed in {elapsed:.1f}s")
        print()

        # Histograms (within-family and cross-family)
        print(f"  Within-family fit count distribution (random systems):")
        for k in range(7):
            pct = within_hist[k] / N_TRIALS * 100
            bar = "#" * int(pct / 2) if pct >= 0.5 else ""
            marker = " <-- real" if k == real_within else ""
            print(f"    count={k}: {within_hist[k]:>9,} ({pct:>6.2f}%) {bar}{marker}")

        print(f"\n  Selected cross-family fit count distribution (random systems):")
        for k in range(4):
            pct = cross_hist[k] / N_TRIALS * 100
            bar = "#" * int(pct / 2) if pct >= 0.5 else ""
            marker = " <-- real" if k == real_cross else ""
            print(f"    count={k}: {cross_hist[k]:>9,} ({pct:>6.2f}%) {bar}{marker}")

        # Chain residual percentiles
        # Filter out infinities (random systems where chain components fell
        # below the grammar floor — chain undefined under v1.0 grammar)
        lepton_arr = np.array(lepton_resids)
        down_arr = np.array(down_resids)
        lepton_finite = lepton_arr[np.isfinite(lepton_arr)]
        down_finite = down_arr[np.isfinite(down_arr)]
        n_lepton_inf = len(lepton_arr) - len(lepton_finite)
        n_down_inf = len(down_arr) - len(down_finite)
        print(f"\n  Random chain residual percentiles (subsample of "
              f"{len(lepton_arr):,}; below-floor excluded):")
        if len(lepton_finite) > 0:
            print(f"    |lepton chain|:  "
                  f"5th={np.percentile(lepton_finite, 5)*100:.3f}%, "
                  f"50th={np.percentile(lepton_finite, 50)*100:.3f}%, "
                  f"95th={np.percentile(lepton_finite, 95)*100:.3f}%  "
                  f"(below-floor excluded: {n_lepton_inf}/{len(lepton_arr)})")
        else:
            print(f"    |lepton chain|:  all subsample chains had below-floor "
                  f"components ({n_lepton_inf}/{len(lepton_arr)})")
        if len(down_finite) > 0:
            print(f"    |down   chain|:  "
                  f"5th={np.percentile(down_finite, 5)*100:.3f}%, "
                  f"50th={np.percentile(down_finite, 50)*100:.3f}%, "
                  f"95th={np.percentile(down_finite, 95)*100:.3f}%  "
                  f"(below-floor excluded: {n_down_inf}/{len(down_arr)})")
        else:
            print(f"    |down   chain|:  all subsample chains had below-floor "
                  f"components ({n_down_inf}/{len(down_arr)})")

        # Joint match result
        print(f"\n  Match-or-exceed (all four features jointly):")
        print(f"    {format_fraction(match_count, N_TRIALS)}")
        print(f"    Band: {band}")
        print()

        results[scale_label] = {
            "real_within": real_within,
            "real_cross": real_cross,
            "real_lepton": real_lepton,
            "real_down": real_down,
            "match_count": match_count,
            "fraction": fraction,
            "band": band,
        }

    # ── Grand summary ──
    print(f"{'━' * 100}")
    print("GRAND SUMMARY — JOINT MATCH-OR-EXCEED FRACTION PER SCALE")
    print(f"{'━' * 100}")
    print()
    print(f"  {'scale':<12} {'real wf':>8} {'real cf':>8} "
          f"{'|lep|%':>10} {'|dn|%':>10} {'match%':>14} {'band':>10}")
    print(f"  {'─'*12} {'─'*8} {'─'*8} {'─'*10} {'─'*10} {'─'*14} {'─'*10}")
    for scale_label in COMMON_SCALES:
        r = results[scale_label]
        if r["match_count"] == 0:
            match_str = f"< {1/N_TRIALS:.0e}"
        else:
            match_str = f"{r['fraction']:.4%}"
        print(f"  {scale_label:<12} {r['real_within']:>4}/6   "
              f"{r['real_cross']:>4}/3   "
              f"{r['real_lepton']*100:>9.3f}% "
              f"{r['real_down']*100:>9.3f}% "
              f"{match_str:>14} {r['band']:>10}")

    print()
    print(f"{'━' * 100}")
    print("AGGREGATE EVIDENCE — N OUT OF 9 SCALES")
    print(f"{'━' * 100}")
    print()
    print(f"  Strong:  {aggregate['strong']}/9")
    print(f"  Modest:  {aggregate['modest']}/9")
    print(f"  Weak:    {aggregate['weak']}/9")
    print(f"  None:    {aggregate['none']}/9")
    print()
    print("Note: aggregate is descriptive only. Random systems are generated")
    print("independently per scale; this is not a combined p-value.")

    total_elapsed = time.time() - t0
    print(f"\nTotal runtime: {total_elapsed:.1f}s")
    print("Done.")
