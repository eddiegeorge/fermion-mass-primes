"""
fermion_mass_primes_test16.py

Multi-scale chain composition (accounting check).

Diagnostic test on whether existing grammar-fitted edges from test11 and
test14 compose coherently along a chosen tree path. Confirms that chosen
edge fits propagate correctly under multiplicative composition.

This is not independent evidence for modal equivalence. The diagnostic
criterion is an arithmetic identity: composed residuals must equal the
residual implied by signed component errors, to numerical tolerance.

Three chains tested at each common scale (PDG excluded):
  Lepton: m_t / m_e = (t/tau) x (tau/mu) x (mu/e)
  Down:   m_t / m_d = (t/b) x (b/s) x (s/d)
  Up:     m_t / m_u = (t/c) x (c/u)

The lepton and down chains pass through cross-family bridges (t/tau, t/b)
that v1.0 found to fit individually at sub-1%. The up chain uses two
within-family edges only.

PDG excluded because the lepton chain breaks at PDG (mu/e fits 2x103
rather than the locked p=211; v1.0 sec 3.1).

See scripts/README.md for full documentation.
"""

import time

import fermion_mass_primes_test11 as t11
import fermion_mass_primes_test14 as t14


# Tolerance for the arithmetic identity check (numerical, not physical)
IDENTITY_TOLERANCE = 1e-9


# ── Chain definitions ──
# Each chain is a list of (edge_label, source, fitter) where:
#   edge_label: human-readable name for the edge (e.g. "t/tau")
#   source: "test11" or "test14" (which fit table to draw from)
#   fitter: key into that table identifying the fitted ratio
#
# Within-family edges live in test11's six-ratio table.
# Cross-family edges live in test14's nine-ratio table.
CHAINS = {
    "lepton (m_t/m_e)": [
        ("t/tau",  "test14", "t/tau"),
        ("tau/mu", "test11", "mtau/mmu"),
        ("mu/e",   "test11", "mmu/me"),
    ],
    "down (m_t/m_d)": [
        ("t/b",    "test14", "t/b"),
        ("b/s",    "test11", "mb/ms"),
        ("s/d",    "test11", "ms/md"),
    ],
    "up (m_t/m_u)": [
        ("t/c",    "test11", "mt/mc"),
        ("c/u",    "test11", "mc/mu"),
    ],
}

# Endpoint observed ratio for each chain — composed from running masses
# at each scale to give the "true" m_t/m_X to compare predictions against.
CHAIN_ENDPOINTS = {
    "lepton (m_t/m_e)": ("t",  "e"),
    "down (m_t/m_d)":   ("t",  "d"),
    "up (m_t/m_u)":     ("t",  "u"),
}

# Common scales only — PDG excluded by design (see module docstring)
COMMON_SCALES = list(t11.ANTUSCH_SCALES.keys())


def get_within_family_fits(scale_label, masses, candidates):
    """Return dict {ratio_name: (R_observed, val_fit, a, p, signed_residual)}
    for all six within-family ratios at the given scale.

    signed_residual = (val_fit - R_observed) / R_observed
    """
    ratios = t11.compute_6_ratios_from_yukawa(masses)
    out = {}
    for name, R in zip(t11.RATIO_NAMES, ratios):
        val, a, p, _abs_res = t14.find_best_match(R, candidates)
        signed_res = (val - R) / R
        out[name] = (R, val, a, p, signed_res)
    return out


def get_cross_family_fits(scale_label, masses_by_particle, candidates):
    """Return dict {pair_label: (R_observed, val_fit, a, p, signed_residual)}
    for all nine cross-family pairs at the given scale.

    R is computed as max(m1,m2)/min(m1,m2) per test14 convention.
    signed_residual = (val_fit - R_observed) / R_observed
    """
    out = {}
    for label, k1, k2 in t14.CROSS_FAMILY_PAIRS:
        m1, m2 = masses_by_particle[k1], masses_by_particle[k2]
        R = max(m1, m2) / min(m1, m2)
        val, a, p, _abs_res = t14.find_best_match(R, candidates)
        signed_res = (val - R) / R
        out[label] = (R, val, a, p, signed_res)
    return out


def compose_chain(chain_steps, within_fits, cross_fits):
    """Walk a chain and return per-edge data plus chain-level totals.

    Returns:
        steps: list of dicts, one per edge, with keys:
            label, source, R_obs, val_fit, a, p, signed_res
        prod_pred:    product of fitted values (predicted m_heavy / m_light)
        prod_obs:     product of observed component ratios
        prod_factor:  product of (1 + delta_i) over all components
        composed_signed_res_direct:  (prod_pred - prod_obs) / prod_obs
        composed_signed_res_implied: prod_factor - 1
        identity_residual: (direct - implied) — should be ~0
    """
    steps = []
    prod_pred = 1.0
    prod_obs = 1.0
    prod_factor = 1.0

    for edge_label, source, fitter_key in chain_steps:
        if source == "test11":
            R_obs, val_fit, a, p, signed_res = within_fits[fitter_key]
        elif source == "test14":
            R_obs, val_fit, a, p, signed_res = cross_fits[fitter_key]
        else:
            raise ValueError(f"unknown source: {source}")

        steps.append({
            "label": edge_label,
            "source": source,
            "R_obs": R_obs,
            "val_fit": val_fit,
            "a": a,
            "p": p,
            "signed_res": signed_res,
        })
        prod_pred *= val_fit
        prod_obs *= R_obs
        prod_factor *= (1.0 + signed_res)

    composed_direct = (prod_pred - prod_obs) / prod_obs
    composed_implied = prod_factor - 1.0
    identity_residual = composed_direct - composed_implied

    return {
        "steps": steps,
        "prod_pred": prod_pred,
        "prod_obs": prod_obs,
        "prod_factor": prod_factor,
        "composed_signed_res_direct": composed_direct,
        "composed_signed_res_implied": composed_implied,
        "identity_residual": identity_residual,
    }


def classify_pattern(steps):
    """Classify the cancellation pattern across a chain's signed residuals.

    Returns one of:
        "uniformly small" — all |delta| < 0.01 (1%)
        "cancelling"      — signs mixed AND at least one |delta| >= 0.01
        "reinforcing"     — all signs same AND at least one |delta| >= 0.01
        "trivial"         — fewer than 2 steps
    """
    if len(steps) < 2:
        return "trivial"
    deltas = [s["signed_res"] for s in steps]
    max_abs = max(abs(d) for d in deltas)
    if max_abs < 0.01:
        return "uniformly small"
    signs = set(1 if d > 0 else -1 for d in deltas if d != 0)
    if len(signs) > 1:
        return "cancelling"
    return "reinforcing"


def format_signed_res(delta):
    """Format a signed residual as a percentage with sign."""
    return f"{delta * 100:+.4f}%"


if __name__ == "__main__":
    t0 = time.time()

    print("=" * 100)
    print("FERMION MASS PRIMES TEST 16")
    print("Multi-scale chain composition (accounting check)")
    print("=" * 100)
    print("Diagnostic check: composed residuals must equal residuals implied")
    print("by signed component errors, to numerical tolerance.")
    print()
    print("A successful diagnostic confirms only that chosen edge fits propagate")
    print("correctly along the tree. It is not evidence for modal equivalence.")
    print()
    print("Chains tested:")
    for chain_name, chain_steps in CHAINS.items():
        edge_labels = " x ".join(s[0] for s in chain_steps)
        print(f"  {chain_name}: {edge_labels}")
    print()
    print(f"Scales: {len(COMMON_SCALES)} common scales (PDG excluded by design)")
    print(f"Identity tolerance: {IDENTITY_TOLERANCE:.0e}")

    # ── Generate candidates (large enough for all scales) ──
    # We need to cover both within-family ratios (largest ~600) and
    # cross-family ratios (largest ~120). Use the same upper bound as
    # test14 to ensure consistency.
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
            m1, m2 = masses[k1], masses[k2]
            cross.append(max(m1, m2) / min(m1, m2))
        all_R.extend(within)
        all_R.extend(cross)
    max_R = max(all_R)
    candidates = t14.generate_candidates(max_R * 2)
    print(f"Candidates generated: {len(candidates)} (2^a x p up to {max_R*2:.0f})")

    # ── Per-scale chain composition ──
    # results[scale][chain_name] = compose_chain output
    results = {}

    for scale_label in COMMON_SCALES:
        raw = t11.ANTUSCH_SCALES[scale_label]
        masses_by_particle = {
            "u": raw["yu"], "c": raw["yc"], "t": raw["yt"],
            "d": raw["yd"], "s": raw["ys"], "b": raw["yb"],
            "e": raw["ye"], "mu": raw["ymu"], "tau": raw["ytau"],
        }
        within_fits = get_within_family_fits(scale_label, raw, candidates)
        cross_fits = get_cross_family_fits(
            scale_label, masses_by_particle, candidates
        )

        scale_results = {}
        for chain_name, chain_steps in CHAINS.items():
            scale_results[chain_name] = compose_chain(
                chain_steps, within_fits, cross_fits
            )

            # Verify chain endpoint matches the direct mass ratio.
            # The product of observed component ratios should equal
            # (m_heavy / m_light) computed directly from the running masses.
            heavy_key, light_key = CHAIN_ENDPOINTS[chain_name]
            direct_endpoint = (
                masses_by_particle[heavy_key] / masses_by_particle[light_key]
            )
            chain_endpoint = scale_results[chain_name]["prod_obs"]
            endpoint_check = (chain_endpoint - direct_endpoint) / direct_endpoint
            scale_results[chain_name]["endpoint_check"] = endpoint_check

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

            # Per-edge breakdown
            print(f"    {'edge':<10} {'source':<8} {'R_obs':>12} "
                  f"{'val_fit':>9} {'form':<10} {'signed_res':>12}")
            print(f"    {'─'*10} {'─'*8} {'─'*12} {'─'*9} {'─'*10} {'─'*12}")
            for step in r["steps"]:
                form = t14.format_match(step["a"], step["p"])
                src_short = step["source"].replace("test", "t")
                print(f"    {step['label']:<10} {src_short:<8} "
                      f"{step['R_obs']:>12.4f} {step['val_fit']:>9} "
                      f"{form:<10} {format_signed_res(step['signed_res']):>12}")

            # Chain-level totals
            print(f"\n    Product of fitted values (predicted): "
                  f"{r['prod_pred']:.6e}")
            print(f"    Product of observed values:           "
                  f"{r['prod_obs']:.6e}")
            print(f"    Product of (1 + delta_i):             "
                  f"{r['prod_factor']:.6f}")
            print()
            print(f"    Composed signed residual (direct):    "
                  f"{format_signed_res(r['composed_signed_res_direct'])}")
            print(f"    Composed signed residual (implied):   "
                  f"{format_signed_res(r['composed_signed_res_implied'])}")
            print(f"    Identity check (direct - implied):    "
                  f"{r['identity_residual']:+.2e}")

            # Diagnostic verdict
            ok = abs(r["identity_residual"]) < IDENTITY_TOLERANCE
            verdict = "OK" if ok else "FAIL"
            print(f"    Diagnostic:                           {verdict}")

            # Endpoint sanity
            ep_ok = abs(r["endpoint_check"]) < IDENTITY_TOLERANCE
            ep_verdict = "OK" if ep_ok else "FAIL"
            print(f"    Endpoint check (chain vs direct):     "
                  f"{r['endpoint_check']:+.2e}  {ep_verdict}")

            # Cancellation pattern
            pattern = classify_pattern(r["steps"])
            print(f"    Cancellation pattern:                 {pattern}")

    # ── Grand summary ──
    print(f"\n{'━' * 100}")
    print("GRAND SUMMARY — COMPOSED RESIDUAL PER CHAIN PER SCALE")
    print(f"{'━' * 100}")
    header = f"\n{'Scale':<12}"
    for chain_name in CHAINS.keys():
        header += f" {chain_name:>22}"
    header += f" {'all OK?':>10}"
    print(header)
    sep = f"{'─' * 12}"
    for _ in CHAINS:
        sep += f" {'─' * 22}"
    sep += f" {'─' * 10}"
    print(sep)

    for scale_label in COMMON_SCALES:
        row = f"{scale_label:<12}"
        all_ok = True
        for chain_name in CHAINS.keys():
            r = results[scale_label][chain_name]
            cell = format_signed_res(r['composed_signed_res_direct'])
            row += f" {cell:>22}"
            if abs(r["identity_residual"]) >= IDENTITY_TOLERANCE:
                all_ok = False
            if abs(r["endpoint_check"]) >= IDENTITY_TOLERANCE:
                all_ok = False
        verdict = "yes" if all_ok else "NO"
        row += f" {verdict:>10}"
        print(row)

    # ── Cancellation pattern summary ──
    print(f"\n{'━' * 100}")
    print("CANCELLATION PATTERN — PER CHAIN PER SCALE")
    print(f"{'━' * 100}")
    header = f"\n{'Scale':<12}"
    for chain_name in CHAINS.keys():
        header += f" {chain_name:>22}"
    print(header)
    sep = f"{'─' * 12}"
    for _ in CHAINS:
        sep += f" {'─' * 22}"
    print(sep)

    for scale_label in COMMON_SCALES:
        row = f"{scale_label:<12}"
        for chain_name in CHAINS.keys():
            pattern = classify_pattern(
                results[scale_label][chain_name]["steps"]
            )
            row += f" {pattern:>22}"
        print(row)

    # ── Identity-check summary ──
    print(f"\n{'━' * 100}")
    print("IDENTITY CHECK SUMMARY")
    print(f"{'━' * 100}")
    print(f"\nTolerance: {IDENTITY_TOLERANCE:.0e}")
    print()

    n_total = len(COMMON_SCALES) * len(CHAINS)
    n_id_ok = 0
    n_ep_ok = 0
    max_id_resid = 0.0
    max_ep_resid = 0.0
    for scale_label in COMMON_SCALES:
        for chain_name in CHAINS.keys():
            r = results[scale_label][chain_name]
            if abs(r["identity_residual"]) < IDENTITY_TOLERANCE:
                n_id_ok += 1
            if abs(r["endpoint_check"]) < IDENTITY_TOLERANCE:
                n_ep_ok += 1
            max_id_resid = max(max_id_resid, abs(r["identity_residual"]))
            max_ep_resid = max(max_ep_resid, abs(r["endpoint_check"]))

    print(f"  Identity check (direct = implied):  {n_id_ok}/{n_total}  "
          f"max |residual| = {max_id_resid:.2e}")
    print(f"  Endpoint check (chain = direct):    {n_ep_ok}/{n_total}  "
          f"max |residual| = {max_ep_resid:.2e}")

    if n_id_ok == n_total and n_ep_ok == n_total:
        print()
        print("  All diagnostics OK: chosen edge fits propagate correctly")
        print("  along the tree. This is a precondition for test19; it is not")
        print("  evidence for modal equivalence.")
    else:
        print()
        print("  ONE OR MORE DIAGNOSTICS FAILED.")
        print("  Investigate bookkeeping, orientation, sign convention, or")
        print("  implementation. Modal-equivalence interpretation requires")
        print("  these diagnostics to be clean before test19 is meaningful.")

    total_elapsed = time.time() - t0
    print(f"\nTotal runtime: {total_elapsed:.1f}s")
    print("Done.")
