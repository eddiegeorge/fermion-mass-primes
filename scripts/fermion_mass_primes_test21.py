#!/usr/bin/env python3
"""
fermion_mass_primes_test21.py

CKM mixing null — random prime assignment under fixed structural template.

Tests whether the specific chain-derived prime values and their slot
placement produce better CKM magnitude predictions than random assignments
under the same structural template.

Model: seven-ingredient impedance network (Z_bal + Z_ch + Z_R + Z_a +
Z_deg + Z_tp). Zero fitted parameters.

See README.md for full specification.
"""

import math
import sys
import time
import numpy as np
from itertools import permutations

# ── Constants (zero fitted parameters) ──────────────────────────────

ALPHA_S = 0.1180       # α_s(M_Z)
PI_HALF = math.pi / 2  # π/2

# ── Observed CKM from Antusch et al. at M_Z ────────────────────────
# UTfit 2023 inputs: sin θ₁₂ = 0.2251, sin θ₂₃ = 0.04193,
# sin θ₁₃ = 0.00370, δ = 1.139 rad.
# Full PDG parameterisation (magnitudes depend on δ for V_td, V_ts).

def build_observed_ckm():
    """Compute observed |V_ij| from standard PDG parameterisation."""
    s12 = 0.2251
    s23 = 0.04193
    s13 = 0.00370
    delta = 1.139

    c12 = math.sqrt(1 - s12**2)
    c23 = math.sqrt(1 - s23**2)
    c13 = math.sqrt(1 - s13**2)

    # Complex phase
    eid = complex(math.cos(delta), math.sin(delta))

    # Standard PDG parameterisation
    V = {}
    V["ud"] = c12 * c13
    V["us"] = s12 * c13
    V["ub"] = s13  # |V_ub| = s13 (phase drops out of magnitude)
    V["cd"] = abs(-s12 * c23 - c12 * s23 * s13 * eid)
    V["cs"] = abs(c12 * c23 - s12 * s23 * s13 * eid)
    V["cb"] = s23 * c13
    V["td"] = abs(s12 * s23 - c12 * c23 * s13 * eid)
    V["ts"] = abs(-c12 * s23 - s12 * c23 * s13 * eid)
    V["tb"] = c23 * c13

    return V

# ── Structural template ────────────────────────────────────────────
# Eight structural slots A-H. Each particle's prime set is determined
# by fixed slot membership. Only the values assigned to slots change.
#
# Real assignment: A=5, B=17, C=53, D=59, E=97, F=127, G=211, H=271

REAL_PRIMES = [5, 17, 53, 59, 97, 127, 211, 271]  # Slots A-H

# Particle mode reconstruction: which slots each particle uses
# Indices into the 8-element prime array [A=0, B=1, ..., H=7]
PARTICLE_SLOTS = {
    "u": [0, 1, 2, 3, 4, 6],        # A B C D E G
    "c": [0, 1, 2, 3, 4, 5, 6],     # A B C D E F G
    "t": [0, 1, 2, 3, 4, 5, 6, 7],  # A B C D E F G H
    "d": [1, 4, 5, 6, 7],           # B E F G H
    "s": [0, 1, 4, 5, 6, 7],        # A B E F G H
    "b": [0, 1, 2, 4, 5, 6, 7],     # A B C E F G H
    "e": [0, 2, 3, 5, 7],           # A C D F H
    "mu": [0, 2, 3, 5, 6, 7],       # A C D F G H
    "tau": [0, 1, 2, 3, 5, 6, 7],   # A B C D F G H
}

# Substrate values
SUBSTRATE_A = {
    "u": 0, "c": 2, "t": 2,
    "d": 0, "s": 2, "b": 2,
    "e": 2, "mu": 2, "tau": 2,
}

# Number of modes per source particle (= number of slots)
N_SRC = {
    "u": 6, "c": 7, "t": 8,
    "d": 5, "s": 6, "b": 7,
}

# CKM transitions: (source, target) for rows and elements
CKM_ROWS = [
    ("u", ["d", "s", "b"]),
    ("c", ["d", "s", "b"]),
    ("t", ["d", "s", "b"]),
]

CKM_LABELS = ["ud", "us", "ub", "cd", "cs", "cb", "td", "ts", "tb"]

# ── Model computation ──────────────────────────────────────────────

def compute_ckm_predictions(primes):
    """
    Compute predicted |V_ij| from mode decomposition using the
    seven-ingredient impedance network model.

    Z_total = Z_bal + Z_ch + Z_R + Z_a + Z_deg + Z_tp

    Args:
        primes: array of 8 values for slots A-H

    Returns:
        dict of |V_ij| magnitudes
    """
    pred = {}

    # Derived constants for this prime assignment
    P = sum(primes)            # total mode energy
    p_max = max(primes)        # largest mode, sets impedance scale

    for src, targets in CKM_ROWS:
        src_slots = set(PARTICLE_SLOTS[src])
        a_src = SUBSTRATE_A[src]
        n_src = N_SRC[src]

        # Compute log-weights (= -Z_total) for numerical stability
        log_weights = []
        for tgt in targets:
            tgt_slots = set(PARTICLE_SLOTS[tgt])
            a_tgt = SUBSTRATE_A[tgt]

            shed_idx = sorted(src_slots - tgt_slots)
            gain_idx = sorted(tgt_slots - src_slots)
            shared_idx = sorted(src_slots & tgt_slots)

            shed_vals = [primes[i] for i in shed_idx]
            gain_vals = [primes[i] for i in gain_idx]
            shared_vals = [primes[i] for i in shared_idx]

            shed_sum = sum(shed_vals)
            gain_sum = sum(gain_vals)
            da = abs(a_tgt - a_src)

            # ── Term 1: Energy balance ──
            Z_bal = ALPHA_S * abs(shed_sum - gain_sum)

            # ── Term 2: Channel conductance ──
            Z_ch = 0.0
            if shed_vals and gain_vals:
                n_shed = len(shed_vals)
                n_gain = len(gain_vals)
                total_K = sum(
                    math.exp(-ALPHA_S * abs(p - q) * max(p, q) / p_max)
                    for p in shed_vals
                    for q in gain_vals
                )
                K = total_K / math.sqrt(n_shed * n_gain)
                if K > 0:
                    Z_ch = -math.log(K)
                else:
                    Z_ch = 1e30  # effectively infinite impedance

            # ── Shared-mode rewiring load J (used in terms 3a, 5) ──
            J = 0.0
            if shared_vals and (shed_vals or gain_vals):
                ln_shared_sum = sum(math.log(p) for p in shared_vals)
                ln_shed_sum = sum(math.log(p) for p in shed_vals) if shed_vals else 0.0
                ln_gain_sum = sum(math.log(p) for p in gain_vals) if gain_vals else 0.0
                denom = ln_shed_sum + ln_gain_sum
                if denom > 0:
                    J = ln_shared_sum / denom

            # ── Term 3: Substrate change ──
            Z_R = 0.0
            if da > 0:
                if a_tgt > a_src:
                    # 3a: Engaging — ring disruption
                    Z_R = PI_HALF * da * J
                else:
                    # 3b: Releasing — cooperative relaxation + gain detuning
                    rewired_sum = sum(shared_vals)  # shared modes rewire
                    gain_detuning = ALPHA_S * gain_sum * da / P if gain_vals else 0.0
                    Z_R = PI_HALF * da * (1 - ALPHA_S * rewired_sum / P + gain_detuning)

            # ── Term 4: Substrate throughput (pump) ──
            Z_a = (ALPHA_S / P) * (shed_sum * a_src + gain_sum * a_tgt)

            # ── Term 5: Channel degradation during substrate change ──
            Z_deg = 0.0
            if da > 0 and shed_vals and gain_vals:
                Z_deg = (ALPHA_S / P) * da * (shed_sum + gain_sum) * J

            # ── Term 6: Ring traversal for pure shedding ──
            Z_tp = 0.0
            if not gain_vals and da == 0:
                # Pure shedding at Δa = 0
                Q = n_src * (n_src - 1)
                if Q > 0:
                    Z_tp = (ALPHA_S / Q) * shed_sum

            Z_total = Z_bal + Z_ch + Z_R + Z_a + Z_deg + Z_tp
            log_weights.append(-Z_total)

        # Softmax normalisation in log space
        max_lw = max(log_weights)
        weights = [math.exp(lw - max_lw) for lw in log_weights]
        total_w = sum(weights)

        for i, tgt in enumerate(targets):
            v2 = weights[i] / total_w
            label = src + tgt
            pred[label] = math.sqrt(v2)

    return pred


def score_predictions(pred, obs):
    """
    Compute primary score S = Σ [ln(|V|²_pred / |V|²_obs)]²
    and secondary score S_max = max |ln(|V|_pred / |V|_obs)|.
    """
    S = 0.0
    max_log_err = 0.0
    for label in CKM_LABELS:
        p = pred[label]
        o = obs[label]
        if p <= 0 or o <= 0:
            return 1e30, 1e30
        log_ratio_v2 = math.log(p**2 / o**2)
        S += log_ratio_v2**2
        log_err = abs(math.log(p / o))
        if log_err > max_log_err:
            max_log_err = log_err
    return S, max_log_err


def count_inversions(pred, obs):
    """
    Count pairwise rank inversions (Kendall tau distance) between
    predicted and observed |V_ij| orderings.
    """
    pred_vals = [pred[label] for label in CKM_LABELS]
    obs_vals = [obs[label] for label in CKM_LABELS]

    # Get rank orderings
    pred_order = sorted(range(9), key=lambda i: -pred_vals[i])
    obs_order = sorted(range(9), key=lambda i: -obs_vals[i])

    # Build rank arrays
    pred_rank = [0] * 9
    obs_rank = [0] * 9
    for rank, idx in enumerate(pred_order):
        pred_rank[idx] = rank
    for rank, idx in enumerate(obs_order):
        obs_rank[idx] = rank

    # Count inversions
    inversions = 0
    for i in range(9):
        for j in range(i + 1, 9):
            if (pred_rank[i] - pred_rank[j]) * (obs_rank[i] - obs_rank[j]) < 0:
                inversions += 1
    return inversions


def check_row_ordering(pred):
    """
    Check row-wise hierarchy:
      |V_ud| > |V_us| > |V_ub|
      |V_cs| > |V_cd| > |V_cb|
      |V_tb| > |V_ts| > |V_td|
    Returns True if all three rows are correctly ordered.
    """
    ok = True
    ok &= pred["ud"] > pred["us"] > pred["ub"]
    ok &= pred["cs"] > pred["cd"] > pred["cb"]
    ok &= pred["tb"] > pred["ts"] > pred["td"]
    return ok


def column_unitarity(pred):
    """
    Compute U_col = max_j |Σ_i |V_ij|² - 1|.
    """
    col_d = pred["ud"]**2 + pred["cd"]**2 + pred["td"]**2
    col_s = pred["us"]**2 + pred["cs"]**2 + pred["ts"]**2
    col_b = pred["ub"]**2 + pred["cb"]**2 + pred["tb"]**2
    return max(abs(col_d - 1), abs(col_s - 1), abs(col_b - 1))


# ── Prime generation ───────────────────────────────────────────────

def sieve_odd_primes(max_val):
    """Return list of odd primes up to max_val."""
    if max_val < 3:
        return []
    sieve = [True] * (max_val + 1)
    sieve[0] = sieve[1] = False
    for i in range(2, int(max_val**0.5) + 1):
        if sieve[i]:
            for j in range(i*i, max_val + 1, i):
                sieve[j] = False
    return [p for p in range(3, max_val + 1) if sieve[p] and p % 2 == 1]


# ── Null populations ───────────────────────────────────────────────

def run_permutation_null(obs, real_pred, real_score, real_S_max, real_ucol, real_inv):
    """Exact enumeration of all 8! permutations of real primes."""
    print("\n" + "=" * 70)
    print("NULL POPULATION 1: PERMUTATION (exact enumeration)")
    print("Same 8 primes, all 8! = 40,320 slot permutations")
    print("=" * 70)

    n_total = 0
    n_match_S = 0
    n_match_joint = 0
    scores = []
    smax_vals = []
    inv_counts = []
    row_order_count = 0
    ucol_vals = []
    best_score = 1e30
    best_perm = None
    best_pred = None
    real_rank = None

    EPS = 1e-12

    for perm in permutations(range(8)):
        primes = [REAL_PRIMES[i] for i in perm]
        pred = compute_ckm_predictions(primes)
        S, S_max = score_predictions(pred, obs)
        inv = count_inversions(pred, obs)
        row_ok = check_row_ordering(pred)
        ucol = column_unitarity(pred)

        n_total += 1
        scores.append(S)
        smax_vals.append(S_max)
        inv_counts.append(inv)
        ucol_vals.append(ucol)
        if row_ok:
            row_order_count += 1

        if S <= real_score + EPS:
            n_match_S += 1

        if S <= real_score + EPS and ucol <= real_ucol + EPS and inv <= real_inv:
            n_match_joint += 1

        if S < best_score:
            best_score = S
            best_perm = primes[:]
            best_pred = pred.copy()

    # Compute real rank
    scores_arr = np.array(scores)
    real_rank = int(np.sum(scores_arr < real_score - EPS)) + 1

    frac_S = n_match_S / n_total
    frac_joint = n_match_joint / n_total

    band = classify_band(frac_S)

    print(f"\n  Real score S = {real_score:.6f}")
    print(f"  Real assignment rank: {real_rank} of {n_total}")
    print(f"  Match-or-exceed (S): {n_match_S}/{n_total} = {frac_S:.6f} [{band}]")
    print(f"  Joint diagnostic: {n_match_joint}/{n_total} = {frac_joint:.6f}")
    print(f"\n  Score distribution:")
    print(f"    Mean:   {np.mean(scores_arr):.3f}")
    print(f"    Median: {np.median(scores_arr):.3f}")
    print(f"    Min:    {np.min(scores_arr):.6f}")
    print(f"    P1:     {np.percentile(scores_arr, 1):.3f}")
    print(f"    P5:     {np.percentile(scores_arr, 5):.3f}")
    print(f"    P25:    {np.percentile(scores_arr, 25):.3f}")

    smax_arr = np.array(smax_vals)
    print(f"\n  Secondary score S_max distribution:")
    print(f"    Real S_max: {real_S_max:.6f}")
    print(f"    Mean:   {np.mean(smax_arr):.3f}")
    print(f"    Median: {np.median(smax_arr):.3f}")
    print(f"    S_max <= real:"
          f" {np.sum(smax_arr <= real_S_max + EPS)}/{n_total}")

    inv_arr = np.array(inv_counts)
    print(f"\n  Inversion count distribution:")
    print(f"    Zero inversions: {np.sum(inv_arr == 0)}/{n_total}"
          f" ({100*np.sum(inv_arr == 0)/n_total:.2f}%)")
    print(f"    Mean inversions: {np.mean(inv_arr):.1f}")
    print(f"    Row-order correct: {row_order_count}/{n_total}"
          f" ({100*row_order_count/n_total:.2f}%)")

    ucol_arr = np.array(ucol_vals)
    print(f"\n  Column unitarity U_col distribution:")
    print(f"    Real U_col: {real_ucol:.4f}")
    print(f"    Mean: {np.mean(ucol_arr):.4f}")
    print(f"    Median: {np.median(ucol_arr):.4f}")
    print(f"    U_col <= real: {np.sum(ucol_arr <= real_ucol + EPS)}/{n_total}")

    if best_pred:
        print(f"\n  Best permutation: {best_perm}")
        print(f"  Best score: {best_score:.6f}")
        print_comparison(best_pred, obs, "Best perm", real_pred=real_pred)

    return frac_S, band


def run_random_null(obs, real_pred, real_score, real_S_max, real_ucol, real_inv,
                    prime_pool, n_trials, label, rng):
    """Monte Carlo null with random prime draws."""
    print(f"\n{'=' * 70}")
    print(f"NULL POPULATION: {label}")
    print(f"8 distinct primes from pool of {len(prime_pool)} odd primes, "
          f"{n_trials:,} trials")
    print(f"{'=' * 70}")

    n_match_S = 0
    n_match_joint = 0
    scores = np.empty(n_trials)
    smax_vals = np.empty(n_trials)
    inv_counts = np.empty(n_trials, dtype=int)
    row_order_count = 0
    ucol_vals = np.empty(n_trials)
    best_score = 1e30
    best_primes = None
    best_pred = None

    EPS = 1e-12
    pool_arr = np.array(prime_pool)

    t0 = time.time()
    report_interval = n_trials // 10

    for trial in range(n_trials):
        # Draw 8 distinct primes uniformly from pool (not sorted)
        idx = rng.choice(len(pool_arr), size=8, replace=False)
        primes = [int(pool_arr[i]) for i in idx]

        pred = compute_ckm_predictions(primes)
        S, S_max = score_predictions(pred, obs)
        inv = count_inversions(pred, obs)
        row_ok = check_row_ordering(pred)
        ucol = column_unitarity(pred)

        scores[trial] = S
        smax_vals[trial] = S_max
        inv_counts[trial] = inv
        ucol_vals[trial] = ucol
        if row_ok:
            row_order_count += 1

        if S <= real_score + EPS:
            n_match_S += 1

        if S <= real_score + EPS and ucol <= real_ucol + EPS and inv <= real_inv:
            n_match_joint += 1

        if S < best_score:
            best_score = S
            best_primes = primes[:]
            best_pred = pred.copy()

        if report_interval > 0 and (trial + 1) % report_interval == 0:
            elapsed = time.time() - t0
            rate = (trial + 1) / elapsed
            remaining = (n_trials - trial - 1) / rate
            print(f"    ... {trial+1:,}/{n_trials:,} "
                  f"({elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)",
                  flush=True)

    elapsed = time.time() - t0

    frac_S = n_match_S / n_trials
    frac_joint = n_match_joint / n_trials
    band = classify_band(frac_S)

    if frac_S == 0:
        frac_str = f"< {1/n_trials:.1e}"
    else:
        frac_str = f"{frac_S:.6f}"

    if frac_joint == 0:
        joint_str = f"< {1/n_trials:.1e}"
    else:
        joint_str = f"{frac_joint:.6f}"

    print(f"\n  Completed in {elapsed:.1f}s")
    print(f"\n  Real score S = {real_score:.6f}")
    print(f"  Match-or-exceed (S): {n_match_S}/{n_trials:,} = {frac_str}"
          f" [{band}]")
    print(f"  Joint diagnostic: {n_match_joint}/{n_trials:,} = {joint_str}")
    print(f"\n  Score distribution:")
    print(f"    Mean:   {np.mean(scores):.3f}")
    print(f"    Median: {np.median(scores):.3f}")
    print(f"    Min:    {np.min(scores):.6f}")
    print(f"    P1:     {np.percentile(scores, 1):.3f}")
    print(f"    P5:     {np.percentile(scores, 5):.3f}")
    print(f"    P25:    {np.percentile(scores, 25):.3f}")

    print(f"\n  Secondary score S_max distribution:")
    print(f"    Real S_max: {real_S_max:.6f}")
    print(f"    Mean:   {np.mean(smax_vals):.3f}")
    print(f"    Median: {np.median(smax_vals):.3f}")
    print(f"    S_max <= real:"
          f" {np.sum(smax_vals <= real_S_max + EPS)}/{n_trials:,}")

    print(f"\n  Inversion count distribution:")
    print(f"    Zero inversions: {np.sum(inv_counts == 0)}/{n_trials:,}"
          f" ({100*np.sum(inv_counts == 0)/n_trials:.4f}%)")
    print(f"    Mean inversions: {np.mean(inv_counts):.1f}")
    print(f"    Row-order correct: {row_order_count}/{n_trials:,}"
          f" ({100*row_order_count/n_trials:.4f}%)")

    print(f"\n  Column unitarity U_col distribution:")
    print(f"    Real U_col: {real_ucol:.4f}")
    print(f"    Mean: {np.mean(ucol_vals):.4f}")
    print(f"    Median: {np.median(ucol_vals):.4f}")
    print(f"    U_col <= real:"
          f" {np.sum(ucol_vals <= real_ucol + EPS)}/{n_trials:,}")

    if best_pred:
        print(f"\n  Best random primes: {best_primes}")
        print(f"  Best score: {best_score:.6f}")
        print_comparison(best_pred, obs, "Best rand", real_pred=real_pred)

    return frac_S, band


# ── Reporting ──────────────────────────────────────────────────────

def classify_band(frac):
    if frac < 0.001:
        return "STRONG"
    elif frac < 0.01:
        return "MODEST"
    elif frac < 0.05:
        return "WEAK"
    else:
        return "NONE"


def print_comparison(pred, obs, label, real_pred=None):
    """Print comparison table. Optionally include real model for three-way."""
    if real_pred:
        print(f"\n  {label} vs real model vs observed:")
        print(f"  {'Element':<10} {label:>10} {'Real':>10} {'Observed':>10}")
        print(f"  {'─'*10} {'─'*10} {'─'*10} {'─'*10}")
        for lbl in CKM_LABELS:
            p = pred[lbl]
            r = real_pred[lbl]
            o = obs[lbl]
            print(f"  |V_{lbl}|{'':<4} {p:>10.5f} {r:>10.5f} {o:>10.5f}")
    else:
        print(f"\n  {label} vs observed:")
        print(f"  {'Element':<10} {'Predicted':>10} {'Observed':>10} {'%err':>8}")
        print(f"  {'─'*10} {'─'*10} {'─'*10} {'─'*8}")
        for lbl in CKM_LABELS:
            p = pred[lbl]
            o = obs[lbl]
            pct = (p - o) / o * 100
            print(f"  |V_{lbl}|{'':<4} {p:>10.5f} {o:>10.5f} {pct:>+7.1f}%")


# ── Main ───────────────────────────────────────────────────────────

def main():
    SEED = 42
    N_TRIALS = 1_000_000
    #N_TRIALS = 100_000
    PRIME_RANGE_MATCHED = 271
    PRIME_RANGE_EXTENDED = 997

    rng = np.random.default_rng(SEED)

    print("=" * 70)
    print("TEST 21: CKM MIXING NULL")
    print("Random prime assignment under fixed structural template")
    print("=" * 70)
    print(f"\n  Fixed constants: α_s(M_Z) = {ALPHA_S}, π/2 = {PI_HALF:.4f}")
    print(f"  Real primes: {REAL_PRIMES}")
    print(f"  Derived: P = {sum(REAL_PRIMES)}, p_max = {max(REAL_PRIMES)}")
    print(f"  n_src: {dict((k, v) for k, v in N_SRC.items())}")
    print(f"  Seed: {SEED}")
    print(f"  Trials per random null: {N_TRIALS:,}")

    # Build observed CKM
    obs = build_observed_ckm()
    print(f"\n  Observed CKM (Antusch M_Z, full PDG parameterisation):")
    for lbl in CKM_LABELS:
        print(f"    |V_{lbl}| = {obs[lbl]:.5f}")

    # Compute real model predictions
    real_pred = compute_ckm_predictions(REAL_PRIMES)
    real_S, real_S_max = score_predictions(real_pred, obs)
    real_inv = count_inversions(real_pred, obs)
    real_row_ok = check_row_ordering(real_pred)
    real_ucol = column_unitarity(real_pred)

    print(f"\n  Real model results:")
    print(f"    Primary score S = {real_S:.6f}")
    print(f"    Secondary score S_max = {real_S_max:.6f}")
    print(f"    Inversions: {real_inv}")
    print(f"    Row ordering correct: {real_row_ok}")
    print(f"    Column unitarity U_col: {real_ucol:.4f}")
    print_comparison(real_pred, obs, "Real model")

    # ── Population 1: Permutation null ──
    frac1, band1 = run_permutation_null(
        obs, real_pred, real_S, real_S_max, real_ucol, real_inv
    )

    # ── Population 2: Range-matched null ──
    pool_matched = sieve_odd_primes(PRIME_RANGE_MATCHED)
    print(f"\n  Range-matched prime pool: {len(pool_matched)} odd primes"
          f" in [3, {PRIME_RANGE_MATCHED}]")
    if len(pool_matched) < 8:
        print("  ERROR: insufficient primes in range-matched pool")
        sys.exit(1)

    frac2, band2 = run_random_null(
        obs, real_pred, real_S, real_S_max, real_ucol, real_inv,
        pool_matched, N_TRIALS,
        f"RANGE-MATCHED [3, {PRIME_RANGE_MATCHED}]", rng
    )

    # ── Population 3: Extended-range null ──
    pool_extended = sieve_odd_primes(PRIME_RANGE_EXTENDED)
    print(f"\n  Extended prime pool: {len(pool_extended)} odd primes"
          f" in [3, {PRIME_RANGE_EXTENDED}]")

    frac3, band3 = run_random_null(
        obs, real_pred, real_S, real_S_max, real_ucol, real_inv,
        pool_extended, N_TRIALS,
        f"EXTENDED [3, {PRIME_RANGE_EXTENDED}]", rng
    )

    # ── Summary ──
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print(f"{'=' * 70}")
    print(f"\n  Real model score: S = {real_S:.6f}")
    print(f"\n  {'Population':<35} {'Match frac':>12} {'Band':>8}")
    print(f"  {'─'*35} {'─'*12} {'─'*8}")
    print(f"  {'Permutation (same primes)':<35}"
          f" {frac1:>12.6f} {band1:>8}")

    if frac2 == 0:
        f2_str = f"< {1/N_TRIALS:.1e}"
    else:
        f2_str = f"{frac2:.6f}"
    print(f"  {f'Range-matched [3,{PRIME_RANGE_MATCHED}]':<35}"
          f" {f2_str:>12} {band2:>8}")

    if frac3 == 0:
        f3_str = f"< {1/N_TRIALS:.1e}"
    else:
        f3_str = f"{frac3:.6f}"
    print(f"  {f'Extended [3,{PRIME_RANGE_EXTENDED}]':<35}"
          f" {f3_str:>12} {band3:>8}")

    print(f"\n  Interpretation:")
    if all(b in ("STRONG", "MODEST") for b in (band1, band2, band3)):
        print("  The CKM result depends on both the specific prime values")
        print("  AND their slot placement.")
    elif band2 in ("STRONG", "MODEST") and band1 in ("WEAK", "NONE"):
        print("  The real prime values are special, but their slot")
        print("  assignment is not tightly constrained.")
    elif band1 in ("STRONG", "MODEST") and band2 in ("WEAK", "NONE"):
        print("  Slot assignment matters but the specific prime values")
        print("  do not — any primes of similar magnitude would work.")
    else:
        print("  No statistically significant structural result.")

    # Verdict
    all_strong = all(b in ("STRONG", "MODEST")
                     for b in (band1, band2, band3))
    if all_strong:
        print(f"\n  VERDICT: The chain-derived prime values and their slot")
        print(f"  placement carry structural information about quark mixing")
        print(f"  that random assignments do not replicate.")
    else:
        print(f"\n  VERDICT: Result is mixed. See per-population bands.")


if __name__ == "__main__":
    main()
