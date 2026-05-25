"""
fermion_mass_primes_test19.py

Held-out cycle closure.

Tests whether held-out cross-family edges — those NOT used to construct
the chain tree — close consistently when predicted by the tree's
composition. If the modal-equivalence reading holds, the tree predicts
a value for every cross-family ratio; comparing tree-predicted values
against observed at every common scale tests whether the modal structure
is internally consistent beyond the edges it was built from.

Tree bridges (used in test16):
  t/tau (lepton chain)
  t/b   (down chain)

Held-out cross-family edges (predicted via composition):
  c/s   = (t/b x b/s) / (t/c)
  d/u   = (t/c x c/u) / (t/b x b/s x s/d)
  u/e   = (t/tau x tau/mu x mu/e) / (t/c x c/u)
  c/mu  = (t/tau x tau/mu) / (t/c)
  d/e   = (t/tau x tau/mu x mu/e) / (t/b x b/s x s/d)
  mu/s  = (t/b x b/s) / (t/tau x tau/mu)
  b/tau = (t/tau) / (t/b)

Seven held-out edges x nine common scales = 63 cells.

Per cell:
  - Compute predicted directed ratio from tree composition
  - Compute observed directed ratio from running masses
  - Signed compositional residual = (pred - obs) / obs
  - Implied compositional residual: prod(1 + delta_i)^s_i - 1
    (s_i = +1 numerator, -1 denominator; delta_i signed component fits)
  - Identity check discrepancy: signed_residual - implied_residual
    (must be ~0 by construction; non-zero indicates implementation bug)
  - Direct grammar fit residual (test14-style nearest 2^a x p),
    or "below floor" if R < 3.

Reading rules (per cell, on |compositional residual|):
  Strong:  <= 5%
  Modest:  5  < x <= 15%
  Weak:    15 < x <= 30%
  None:    > 30%

PDG excluded — convention-mixing breaks the lepton chain at PDG.

See scripts/README.md for full documentation.
"""

import time

import fermion_mass_primes_test11 as t11
import fermion_mass_primes_test14 as t14


# Tolerance for the identity check (numerical, not physical)
IDENTITY_TOLERANCE = 1e-9

# Threshold for "direct grammar fit small" in 2x2 matrix
DIRECT_FIT_SMALL_THRESHOLD = 3.0  # percent

# Grammar floor (matches test14)
GRAMMAR_FLOOR = 3.0


# ── Held-out edges, each as (signed factors) where each factor is
#    (edge_label, source, fitter_key, sign)
#    sign = +1 if numerator, -1 if denominator
#    source = "test11" or "test14"
#
# Each edge is identified by its directed label "h/l" — m_h / m_l where
# h is the first particle in the label. (e.g., "b/tau" = m_b / m_tau.)
#
# The fitter_keys reference the within-family fits (test11 RATIO_NAMES)
# and cross-family fits (test14 CROSS_FAMILY_PAIRS labels). Test11 keys
# always use the heavier/lighter convention "mc/mu"; test14 keys are
# the labels in CROSS_FAMILY_PAIRS, which are also heavier/lighter.
HELD_OUT_EDGES = {
    "c/s":   {
        "directed": ("c", "s"),
        "factors": [
            ("t/b",   "test14", "t/b",     +1),
            ("b/s",   "test11", "mb/ms",   +1),
            ("t/c",   "test11", "mt/mc",   -1),
        ],
    },
    "d/u":   {
        "directed": ("d", "u"),
        "factors": [
            ("t/c",   "test11", "mt/mc",   +1),
            ("c/u",   "test11", "mc/mu",   +1),
            ("t/b",   "test14", "t/b",     -1),
            ("b/s",   "test11", "mb/ms",   -1),
            ("s/d",   "test11", "ms/md",   -1),
        ],
    },
    "u/e":   {
        "directed": ("u", "e"),
        "factors": [
            ("t/tau",   "test14", "t/tau",      +1),
            ("tau/mu",  "test11", "mtau/mmu",   +1),
            ("mu/e",    "test11", "mmu/me",     +1),
            ("t/c",     "test11", "mt/mc",      -1),
            ("c/u",     "test11", "mc/mu",      -1),
        ],
    },
    "c/mu":  {
        "directed": ("c", "mu"),
        "factors": [
            ("t/tau",   "test14", "t/tau",      +1),
            ("tau/mu",  "test11", "mtau/mmu",   +1),
            ("t/c",     "test11", "mt/mc",      -1),
        ],
    },
    "d/e":   {
        "directed": ("d", "e"),
        "factors": [
            ("t/tau",   "test14", "t/tau",      +1),
            ("tau/mu",  "test11", "mtau/mmu",   +1),
            ("mu/e",    "test11", "mmu/me",     +1),
            ("t/b",     "test14", "t/b",        -1),
            ("b/s",     "test11", "mb/ms",      -1),
            ("s/d",     "test11", "ms/md",      -1),
        ],
    },
    "mu/s":  {
        "directed": ("mu", "s"),
        "factors": [
            ("t/b",     "test14", "t/b",        +1),
            ("b/s",     "test11", "mb/ms",      +1),
            ("t/tau",   "test14", "t/tau",      -1),
            ("tau/mu",  "test11", "mtau/mmu",   -1),
        ],
    },
    "b/tau": {
        "directed": ("b", "tau"),
        "factors": [
            ("t/tau",   "test14", "t/tau",      +1),
            ("t/b",     "test14", "t/b",        -1),
        ],
    },
}

HELD_OUT_ORDER = ["c/s", "d/u", "u/e", "c/mu", "d/e", "mu/s", "b/tau"]

# Map test14 cross-family labels back to (heavier, lighter) particle keys.
# We need this to compute the heavier/lighter ratio at each scale that
# test14 actually fitted (since CROSS_FAMILY_PAIRS uses max/min orientation).
T14_LABEL_TO_KEYS = {label: (k1, k2) for label, k1, k2 in t14.CROSS_FAMILY_PAIRS}

# Within-family pair convention (test11 names): "{heavier_short}/{lighter_short}"
# Map each test11 ratio name to (heavier_key, lighter_key).
T11_NAME_TO_KEYS = {
    "mc/mu":     ("c",   "u"),
    "mt/mc":     ("t",   "c"),
    "ms/md":     ("s",   "d"),
    "mb/ms":     ("b",   "s"),
    "mmu/me":    ("mu",  "e"),
    "mtau/mmu":  ("tau", "mu"),
}

COMMON_SCALES = list(t11.ANTUSCH_SCALES.keys())


def get_masses(yukawa_dict):
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


def fit_within_family(masses, candidates):
    """Return dict {test11_name: (R_obs, val_fit, a, p, signed_res)}.

    `masses` is keyed by particle short name (u, c, t, d, s, b, e, mu, tau).
    R_obs is the heavier/lighter ratio per test11 convention.
    signed_res = (val_fit - R_obs) / R_obs.
    """
    out = {}
    for name, (h, l) in T11_NAME_TO_KEYS.items():
        R = masses[h] / masses[l]
        val, a, p, _abs_res = t14.find_best_match(R, candidates)
        signed_res = (val - R) / R
        out[name] = (R, val, a, p, signed_res)
    return out


def fit_cross_family(masses, candidates):
    """Return dict {label: (R_obs, val_fit, a, p, signed_res, orient)}.

    R_obs is max(m1,m2)/min(m1,m2) per test14 convention.
    orient is "k1>k2" or "k2>k1".
    signed_res = (val_fit - R_obs) / R_obs.
    """
    out = {}
    for label, k1, k2 in t14.CROSS_FAMILY_PAIRS:
        m1, m2 = masses[k1], masses[k2]
        if m1 >= m2:
            R = m1 / m2
            orient = f"{k1}>{k2}"
        else:
            R = m2 / m1
            orient = f"{k2}>{k1}"
        val, a, p, _abs_res = t14.find_best_match(R, candidates)
        signed_res = (val - R) / R
        out[label] = (R, val, a, p, signed_res, orient)
    return out


def directed_signed_factor(label, source, fitter_key, sign,
                           within_fits, cross_fits, masses):
    """Get the signed-residual contribution of one factor in a tree composition.

    The held-out edges work in directed (heavier/lighter or mixed) form, but
    test11 and test14 fit ratios in heavier/lighter form. For each factor
    we therefore need:
      - the directed value the composition uses (e.g., "t/c" = m_t/m_c)
      - the corresponding test11/test14 fit value
      - the signed residual relative to whichever ratio we're using

    All within-family fits store mh/ml ratios where mh > ml, by construction
    (test11 RATIO_NAMES are written heavier-over-lighter and the masses
    enforce this). Cross-family fits store max/min, also heavier/lighter.

    The held-out edge labels happen to use heavier/lighter convention too
    (e.g., "t/c", "b/s"), so the directed factor in the composition matches
    the test11/test14 fit orientation directly. We just need to know:
      - factor_value = m[h] / m[l]  (matching the fit)
      - factor_fit   = val_fit
      - factor_delta = signed_res

    sign tells us whether this factor enters as numerator (+1) or
    denominator (-1) in the composition.
    """
    if source == "test11":
        R, val, _a, _p, signed_res = within_fits[fitter_key]
        return R, val, signed_res
    elif source == "test14":
        R, val, _a, _p, signed_res, _orient = cross_fits[fitter_key]
        return R, val, signed_res
    else:
        raise ValueError(f"unknown source: {source}")


def compose_held_out(edge_label, edge_def, within_fits, cross_fits, masses):
    """Compute everything for one held-out edge at one scale.

    Returns dict with:
        directed_obs:       observed directed ratio m_h / m_l
        directed_pred:      tree-predicted directed ratio
        signed_resid:       (pred - obs) / obs
        implied_resid:      prod((1+delta_i)^s_i) - 1
        identity_disc:      signed_resid - implied_resid  (~0 by construction)
        factors:            list of dicts (per factor) for output
        direct_fit:         dict with R_obs, val_fit, a, p, signed_res, orient
                            from test14 (or None if below grammar floor)
        below_floor:        True if observed directed ratio in test14 form < 3
    """
    h, l = edge_def["directed"]

    # Observed directed ratio
    directed_obs = masses[h] / masses[l]

    # Compose predicted ratio and propagate signed errors
    directed_pred = 1.0
    prod_factor = 1.0
    factor_records = []
    for (label, source, fitter_key, sign) in edge_def["factors"]:
        R_obs_factor, val_fit, signed_res = directed_signed_factor(
            label, source, fitter_key, sign,
            within_fits, cross_fits, masses
        )
        if sign == +1:
            directed_pred *= val_fit
            prod_factor *= (1.0 + signed_res)
        else:
            directed_pred /= val_fit
            prod_factor /= (1.0 + signed_res)
        factor_records.append({
            "label": label,
            "source": source,
            "fitter_key": fitter_key,
            "sign": sign,
            "R_obs": R_obs_factor,
            "val_fit": val_fit,
            "signed_res": signed_res,
        })

    signed_resid = (directed_pred - directed_obs) / directed_obs
    implied_resid = prod_factor - 1.0
    identity_disc = signed_resid - implied_resid

    # Direct grammar fit comparison (from test14's existing fit)
    # Look up the cross-family label in test14: it's the held-out edge label
    # if oriented heavier/lighter, but might be inverted otherwise.
    # All seven held-out labels are also test14 labels by direct match.
    if edge_label in T14_LABEL_TO_KEYS:
        df_R, df_val, df_a, df_p, df_signed, df_orient = cross_fits[edge_label]
        below_floor = df_R < GRAMMAR_FLOOR
        direct_fit = {
            "R_obs": df_R,
            "val_fit": df_val,
            "a": df_a,
            "p": df_p,
            "signed_res": df_signed,
            "abs_res_pct": abs(df_signed) * 100.0,
            "orient": df_orient,
            "below_floor": below_floor,
        }
    else:
        direct_fit = None
        below_floor = False

    return {
        "directed_obs": directed_obs,
        "directed_pred": directed_pred,
        "signed_resid": signed_resid,
        "implied_resid": implied_resid,
        "identity_disc": identity_disc,
        "factors": factor_records,
        "direct_fit": direct_fit,
        "below_floor": below_floor,
    }


def evidence_band(abs_resid_pct):
    """Pre-committed reading rules on |compositional residual|."""
    if abs_resid_pct <= 5.0:
        return "strong"
    if abs_resid_pct <= 15.0:
        return "modest"
    if abs_resid_pct <= 30.0:
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


def format_signed_pct(delta):
    return f"{delta * 100:+.4f}%"


if __name__ == "__main__":
    t0 = time.time()

    print("=" * 100)
    print("FERMION MASS PRIMES TEST 19")
    print("Held-out cycle closure")
    print("=" * 100)
    print("Graph-closure test on the modal-equivalence reading. Tree built from")
    print("two cross-family bridges (t/tau, t/b) plus within-family edges.")
    print("Predicts seven held-out cross-family edges; tests whether predictions")
    print("match observation across nine common scales.")
    print()
    print("Reading rules (per cell, |compositional residual|):")
    print("  Strong:  <= 5%")
    print("  Modest:  5  < x <= 15%")
    print("  Weak:    15 < x <= 30%")
    print("  None:    > 30%")
    print()
    print(f"Identity check tolerance: {IDENTITY_TOLERANCE:.0e}")
    print(f"Direct fit small threshold: {DIRECT_FIT_SMALL_THRESHOLD}%")
    print(f"Grammar floor: R = {GRAMMAR_FLOOR}")
    print()
    print("Held-out edges:")
    for label in HELD_OUT_ORDER:
        parts = []
        for i, (lab, _src, _key, s) in enumerate(HELD_OUT_EDGES[label]["factors"]):
            if i == 0:
                parts.append(lab)
            else:
                op = "x" if s == +1 else "/"
                parts.append(f"{op} {lab}")
        print(f"  {label:<8} = {' '.join(parts)}")
    print()
    print(f"Scales: {len(COMMON_SCALES)} common scales (PDG excluded by design)")

    # ── Generate candidates (same upper bound as test14) ──
    all_R = []
    for scale_label in COMMON_SCALES:
        raw = t11.ANTUSCH_SCALES[scale_label]
        within = t11.compute_6_ratios_from_yukawa(raw)
        masses = get_masses(raw)
        cross = []
        for _label, k1, k2 in t14.CROSS_FAMILY_PAIRS:
            cross.append(max(masses[k1], masses[k2]) / min(masses[k1], masses[k2]))
        all_R.extend(within)
        all_R.extend(cross)
    max_R = max(all_R)
    candidates = t14.generate_candidates(max_R * 2)
    print(f"Candidates generated: {len(candidates)} (2^a x p up to {max_R*2:.0f})")

    # ── Per-scale, per-edge analysis ──
    # results[scale][edge_label] = compose_held_out output
    results = {}
    for scale_label in COMMON_SCALES:
        raw = t11.ANTUSCH_SCALES[scale_label]
        masses = get_masses(raw)
        within_fits = fit_within_family(masses, candidates)
        cross_fits = fit_cross_family(masses, candidates)
        scale_results = {}
        for edge_label in HELD_OUT_ORDER:
            scale_results[edge_label] = compose_held_out(
                edge_label, HELD_OUT_EDGES[edge_label],
                within_fits, cross_fits, masses
            )
        results[scale_label] = scale_results

    # ── Per-scale, per-edge detailed output ──
    for scale_label in COMMON_SCALES:
        print(f"\n{'━' * 100}")
        print(f"SCALE: {scale_label}")
        print(f"{'━' * 100}")

        for edge_label in HELD_OUT_ORDER:
            r = results[scale_label][edge_label]
            band = evidence_band(abs(r["signed_resid"]) * 100)
            df = r["direct_fit"]

            print(f"\n  EDGE: {edge_label}")
            print(f"  {'─' * 96}")

            # Per-factor breakdown
            print(f"    Factors used in composition:")
            print(f"    {'label':<10} {'source':<8} {'sign':>5} {'R_obs':>12} "
                  f"{'val_fit':>9} {'signed_res':>13}")
            print(f"    {'─'*10} {'─'*8} {'─'*5} {'─'*12} {'─'*9} {'─'*13}")
            for f in r["factors"]:
                src_short = "t" + f["source"][-2:]  # "t11" or "t14"
                sign_str = "+1" if f["sign"] == +1 else "-1"
                print(f"    {f['label']:<10} {src_short:<8} {sign_str:>5} "
                      f"{f['R_obs']:>12.4f} {f['val_fit']:>9} "
                      f"{format_signed_pct(f['signed_res']):>13}")

            # Predicted vs observed
            print(f"\n    Directed observed:  {r['directed_obs']:.6e}")
            print(f"    Directed predicted: {r['directed_pred']:.6e}")
            print(f"    Signed compositional residual:  "
                  f"{format_signed_pct(r['signed_resid'])}")
            print(f"    Implied compositional residual: "
                  f"{format_signed_pct(r['implied_resid'])}")
            print(f"    Identity check discrepancy:     "
                  f"{r['identity_disc']:+.2e}")

            id_ok = abs(r["identity_disc"]) < IDENTITY_TOLERANCE
            print(f"    Identity check: {'OK' if id_ok else 'FAIL — implementation error'}")

            # Direct fit comparison
            if df is None:
                print(f"    Direct grammar fit: not in test14 set")
            elif df["below_floor"]:
                print(f"    Direct grammar fit: below floor "
                      f"(R_obs={df['R_obs']:.4f} < {GRAMMAR_FLOOR})")
            else:
                form = format_match(df["a"], df["p"])
                print(f"    Direct grammar fit: orient={df['orient']}, "
                      f"R_obs={df['R_obs']:.4f}, fit={df['val_fit']} ({form}), "
                      f"residual={df['abs_res_pct']:.4f}%")

            # 2x2 quadrant
            comp_small = abs(r["signed_resid"]) * 100 <= 5.0
            if df is None:
                quadrant_note = "(direct fit not applicable)"
            elif df["below_floor"]:
                quadrant_note = (
                    "(direct fit below floor; only compositional column applies)"
                )
            elif df["abs_res_pct"] <= DIRECT_FIT_SMALL_THRESHOLD:
                if comp_small:
                    quadrant_note = "[direct small AND compositional small]"
                else:
                    quadrant_note = (
                        "[direct small AND compositional large — challenges tree]"
                    )
            else:
                if comp_small:
                    quadrant_note = (
                        "[direct large/none AND compositional small — "
                        "tree predicts better than local grammar candidate]"
                    )
                else:
                    quadrant_note = "[direct large AND compositional large]"

            print(f"    Compositional band: {band}  "
                  f"(|residual| = {abs(r['signed_resid'])*100:.4f}%)")
            print(f"    Quadrant: {quadrant_note}")

    # ── Grand summary: per-cell band table ──
    print(f"\n{'━' * 110}")
    print("GRAND SUMMARY — COMPOSITIONAL RESIDUAL BAND PER CELL")
    print(f"{'━' * 110}")
    header = f"\n{'Scale':<12}"
    for edge in HELD_OUT_ORDER:
        header += f" {edge:>13}"
    print(header)
    sep = f"{'─'*12}"
    for _ in HELD_OUT_ORDER:
        sep += f" {'─'*13}"
    print(sep)

    for scale_label in COMMON_SCALES:
        row = f"{scale_label:<12}"
        for edge in HELD_OUT_ORDER:
            r = results[scale_label][edge]
            abs_resid_pct = abs(r["signed_resid"]) * 100
            band = evidence_band(abs_resid_pct)
            cell = f"{abs_resid_pct:>5.1f}% {band:<6}"
            row += f" {cell:>13}"
        print(row)

    # ── Identity check summary ──
    print(f"\n{'━' * 100}")
    print("IDENTITY CHECK SUMMARY (bookkeeping diagnostic)")
    print(f"{'━' * 100}")
    n_total = len(COMMON_SCALES) * len(HELD_OUT_ORDER)
    n_id_ok = 0
    max_id_disc = 0.0
    for scale_label in COMMON_SCALES:
        for edge in HELD_OUT_ORDER:
            r = results[scale_label][edge]
            if abs(r["identity_disc"]) < IDENTITY_TOLERANCE:
                n_id_ok += 1
            max_id_disc = max(max_id_disc, abs(r["identity_disc"]))
    print(f"  Identity check: {n_id_ok}/{n_total}  "
          f"max |discrepancy| = {max_id_disc:.2e}")
    print(f"  Tolerance: {IDENTITY_TOLERANCE:.0e}")
    if n_id_ok != n_total:
        print()
        print("  ONE OR MORE IDENTITY CHECKS FAILED.")
        print("  Investigate bookkeeping, orientation, or sign convention.")
    else:
        print("  All identity checks pass — bookkeeping is correct.")

    # ── Cell-level aggregate (N/63) ──
    print(f"\n{'━' * 100}")
    print("AGGREGATE — CELL-LEVEL BANDS (N/63)")
    print(f"{'━' * 100}")
    band_counts = {"strong": 0, "modest": 0, "weak": 0, "none": 0}
    for scale_label in COMMON_SCALES:
        for edge in HELD_OUT_ORDER:
            r = results[scale_label][edge]
            band_counts[evidence_band(abs(r["signed_resid"]) * 100)] += 1
    print(f"  Strong:  {band_counts['strong']}/{n_total}")
    print(f"  Modest:  {band_counts['modest']}/{n_total}")
    print(f"  Weak:    {band_counts['weak']}/{n_total}")
    print(f"  None:    {band_counts['none']}/{n_total}")

    # ── Per-edge aggregate (N/9) ──
    print(f"\n{'━' * 100}")
    print("AGGREGATE — PER-EDGE BANDS (N/9)")
    print(f"{'━' * 100}")
    print()
    print(f"  {'edge':<8} {'strong':>8} {'modest':>8} {'weak':>8} {'none':>8}  "
          f"{'mean |resid|':>14}")
    print(f"  {'─'*8} {'─'*8} {'─'*8} {'─'*8} {'─'*8}  {'─'*14}")
    for edge in HELD_OUT_ORDER:
        edge_bands = {"strong": 0, "modest": 0, "weak": 0, "none": 0}
        edge_resids = []
        for scale_label in COMMON_SCALES:
            r = results[scale_label][edge]
            abs_pct = abs(r["signed_resid"]) * 100
            edge_bands[evidence_band(abs_pct)] += 1
            edge_resids.append(abs_pct)
        mean_resid = sum(edge_resids) / len(edge_resids)
        print(f"  {edge:<8} {edge_bands['strong']:>7}/9 "
              f"{edge_bands['modest']:>7}/9 "
              f"{edge_bands['weak']:>7}/9 "
              f"{edge_bands['none']:>7}/9  "
              f"{mean_resid:>13.2f}%")

    # ── Challenging quadrant cells ──
    print(f"\n{'━' * 100}")
    print("CHALLENGING QUADRANT — direct fit small (<= 3%) AND compositional large (> 5%)")
    print(f"{'━' * 100}")
    challenging = []
    for scale_label in COMMON_SCALES:
        for edge in HELD_OUT_ORDER:
            r = results[scale_label][edge]
            df = r["direct_fit"]
            if df is None or df["below_floor"]:
                continue
            comp_pct = abs(r["signed_resid"]) * 100
            if df["abs_res_pct"] <= DIRECT_FIT_SMALL_THRESHOLD and comp_pct > 5.0:
                challenging.append((scale_label, edge, comp_pct, df["abs_res_pct"]))
    if not challenging:
        print("  None.")
    else:
        print(f"\n  {'scale':<12} {'edge':<8} {'comp resid':>12} {'direct resid':>14}")
        print(f"  {'─'*12} {'─'*8} {'─'*12} {'─'*14}")
        for scale_label, edge, comp_pct, direct_pct in challenging:
            print(f"  {scale_label:<12} {edge:<8} "
                  f"{comp_pct:>11.2f}% {direct_pct:>13.2f}%")

    # ── Pass criterion verdict ──
    print(f"\n{'━' * 100}")
    print("CHAIN HYPOTHESIS — INTERPRETATION")
    print(f"{'━' * 100}")
    n_strong = band_counts["strong"]
    n_modest = band_counts["modest"]
    n_weak = band_counts["weak"]
    n_none = band_counts["none"]
    print()
    print(f"  Strong+modest: {n_strong + n_modest}/{n_total}")
    print(f"  Weak+none:     {n_weak + n_none}/{n_total}")
    print()
    if n_weak + n_none == 0:
        print("  All cells in strong/modest band. Chain hypothesis supported")
        print("  on closure across the held-out cross-family graph.")
    elif n_strong + n_modest >= 2 * (n_weak + n_none):
        print("  Cells are predominantly strong/modest. Chain hypothesis")
        print("  supported with caveats; failure pattern needs inspection.")
    else:
        print("  Cells are systematically weak/none. Chain hypothesis is")
        print("  challenged on closure across the held-out graph.")
        print("  Inspect failure pattern (per-edge or per-scale clustering).")

    total_elapsed = time.time() - t0
    print(f"\nTotal runtime: {total_elapsed:.1f}s")
    print("Done.")
