"""
Exploratory OPC refinements (RECD-compatible primitives).

These are design/reference implementations for proposals in
`docs/OPC_RECD_LEVEL3_REFINEMENT_PROPOSALS.md`.

Constraints:
- Discrete / ordinal only — no basal μ/σ z-scores on continuous series.
- Does NOT replace frozen abs-z.
- Does NOT claim clinical utility.
- OPS (ordinal synergy surplus) is a Level-3–consistent *proxy*, not excess3 itself.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from .opc_detector import _validate_symbols


def opc_detect_gap_tolerant(
    sigma: np.ndarray,
    *,
    L: int = 50,
    theta_D: float = 0.35,
    theta_R: int = 5,
    G: int = 1,
    K: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    OPC with gap-tolerant persistence (Proposal A / OPC-G).

    Persistence credit increments on low-diversity windows. Up to G consecutive
    high-diversity windows are tolerated without zeroing credit; more than G
    consecutive high-div windows reset credit to 0.

    Alarm when credit >= theta_R and the current window is low-diversity
    (alarm only on collapse, not on pure gap fill).
    """
    if L < 2:
        raise ValueError("L must be >= 2")
    if not (0.0 < theta_D <= 1.0):
        raise ValueError("theta_D must be in (0, 1]")
    if theta_R < 1:
        raise ValueError("theta_R must be >= 1")
    if G < 0:
        raise ValueError("G must be >= 0")

    if K is None:
        arr = np.asarray(sigma)
        K = int(np.max(arr)) + 1 if arr.size else 2
        K = max(K, 2)

    sigma = _validate_symbols(sigma, K)
    T = sigma.size
    alarm = np.zeros(T, dtype=np.int8)
    diversity = np.full(T, np.nan, dtype=float)
    credit = np.zeros(T, dtype=np.int64)
    low_div = np.zeros(T, dtype=np.int8)

    run_credit = 0
    gap_count = 0
    for t in range(T):
        if t < L - 1:
            continue
        window = sigma[t - L + 1 : t + 1]
        D = len(np.unique(window)) / float(K)
        diversity[t] = D
        if D <= theta_D:
            low_div[t] = 1
            run_credit += 1
            gap_count = 0
        else:
            if run_credit > 0 and gap_count < G:
                gap_count += 1
                # credit held (gap-tolerant)
            else:
                run_credit = 0
                gap_count = 0
        credit[t] = run_credit
        if low_div[t] == 1 and run_credit >= theta_R:
            alarm[t] = 1

    return {
        "alarm": alarm,
        "diversity": diversity,
        "persistence": credit,
        "low_div": low_div,
    }


def opc_detect_basal_relative(
    sigma: np.ndarray,
    *,
    L: int = 50,
    theta_D: float = 0.35,
    theta_R: int = 5,
    basal_end: int,
    rho: float = 0.85,
    K: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    OPC with basal-relative collapse gate (Proposal C / OPC-BR).

    Low-diversity requires both absolute D_t <= theta_D and
    D_t <= rho * D_basal, where D_basal is support diversity on sigma[0:basal_end].
    """
    if basal_end < L:
        raise ValueError("basal_end must be >= L")
    if not (0.0 < rho <= 1.0):
        raise ValueError("rho must be in (0, 1]")

    if K is None:
        arr = np.asarray(sigma)
        K = int(np.max(arr)) + 1 if arr.size else 2
        K = max(K, 2)

    sigma = _validate_symbols(sigma, K)
    basal = sigma[:basal_end]
    # diversity of entire basal segment as reference repertoire richness
    D_basal = len(np.unique(basal)) / float(K)
    rel_cap = min(theta_D, rho * D_basal)

    T = sigma.size
    alarm = np.zeros(T, dtype=np.int8)
    diversity = np.full(T, np.nan, dtype=float)
    persistence = np.zeros(T, dtype=np.int64)
    low_div = np.zeros(T, dtype=np.int8)

    run = 0
    for t in range(T):
        if t < max(L - 1, basal_end):
            # search masked during basal (same spirit as SDD)
            continue
        window = sigma[t - L + 1 : t + 1]
        D = len(np.unique(window)) / float(K)
        diversity[t] = D
        if D <= rel_cap:
            low_div[t] = 1
            run += 1
        else:
            run = 0
        persistence[t] = run
        if low_div[t] == 1 and run >= theta_R:
            alarm[t] = 1

    return {
        "alarm": alarm,
        "diversity": diversity,
        "persistence": persistence,
        "low_div": low_div,
        "D_basal": np.array([D_basal]),
        "rel_cap": np.array([rel_cap]),
    }


def ordinal_synergy_surplus(
    pi1: np.ndarray,
    pi2: np.ndarray,
    *,
    L: int,
    k1: int,
    k2: int,
) -> np.ndarray:
    """
    Windowed TV distance between joint law and product of margins (Proposal B / OPS).

    S_t = TV(P_joint, P1 ⊗ P2) on the last L samples ending at t.
    Returns array of length T with NaN where undefined.
    Purely discrete — no continuous moments.
    """
    pi1 = np.asarray(pi1, dtype=np.int64).ravel()
    pi2 = np.asarray(pi2, dtype=np.int64).ravel()
    if pi1.shape != pi2.shape:
        raise ValueError("pi1 and pi2 must have the same shape")
    if L < 2:
        raise ValueError("L must be >= 2")
    if k1 < 2 or k2 < 2:
        raise ValueError("k1 and k2 must be >= 2")
    if pi1.size and (pi1.min() < 0 or pi1.max() >= k1):
        raise ValueError("pi1 symbols out of range")
    if pi2.size and (pi2.min() < 0 or pi2.max() >= k2):
        raise ValueError("pi2 symbols out of range")

    T = pi1.size
    S = np.full(T, np.nan, dtype=float)
    if T < L:
        return S

    # Rolling joint histogram: O(T · k1·k2) instead of O(T · L · k1·k2).
    # Mathematically identical to rebuilding counts each window.
    joint_codes = pi1 * int(k2) + pi2
    counts = np.zeros(int(k1) * int(k2), dtype=np.float64)
    inv_L = 1.0 / float(L)
    for i in range(L):
        counts[joint_codes[i]] += 1.0
    for t in range(L - 1, T):
        if t >= L:
            counts[joint_codes[t - L]] -= 1.0
            counts[joint_codes[t]] += 1.0
        joint = counts.reshape(int(k1), int(k2)) * inv_L
        p1 = joint.sum(axis=1)
        p2 = joint.sum(axis=0)
        prod = np.outer(p1, p2)
        S[t] = 0.5 * np.abs(joint - prod).sum()
    return S


def _mad_scale(x: np.ndarray) -> float:
    """Robust scale: 1.4826 * median(|x - median(x)|). Floor tiny values."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return 0.0
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale < 1e-12:
        # Fall back to std if MAD collapses (near-constant basal)
        std = float(np.std(x))
        return std if np.isfinite(std) and std > 1e-12 else 1e-6
    return scale


def compute_surplus_threshold(
    basal_vals: np.ndarray,
    *,
    basal_mode: str = "mean_delta",
    theta_delta_S: Optional[float] = None,
    basal_q: float = 90.0,
    mad_kappa: float = 2.5,
) -> Tuple[float, float, str]:
    """
    Map basal surplus samples → (reference_stat, absolute_threshold, mode).

    Modes
    -----
    mean_delta : thr = mean(basal) + theta_delta_S  (legacy I0)
    percentile : thr = percentile(basal, basal_q)
    mad        : thr = median(basal) + mad_kappa * MAD

    Returns (S_ref, thr, mode_resolved). For mean_delta, S_ref is the mean;
    for percentile/mad, S_ref is median (diagnostic anchor).
    """
    mode = str(basal_mode).lower().strip()
    vals = np.asarray(basal_vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        raise ValueError("no finite basal S values")
    if mode == "mean_delta":
        if theta_delta_S is None:
            raise ValueError("theta_delta_S required when basal_mode=mean_delta")
        s_ref = float(np.mean(vals))
        thr = s_ref + float(theta_delta_S)
        return s_ref, thr, mode
    if mode in ("percentile", "quantile", "q"):
        q = float(basal_q)
        if not (0.0 < q < 100.0):
            raise ValueError("basal_q must be in (0, 100)")
        s_ref = float(np.median(vals))
        thr = float(np.percentile(vals, q))
        return s_ref, thr, "percentile"
    if mode == "mad":
        if float(mad_kappa) <= 0:
            raise ValueError("mad_kappa must be > 0")
        s_ref = float(np.median(vals))
        thr = s_ref + float(mad_kappa) * _mad_scale(vals)
        return s_ref, thr, mode
    raise ValueError(
        "basal_mode must be one of: mean_delta, percentile, mad"
    )


def ops_detect(
    pi1: np.ndarray,
    pi2: np.ndarray,
    *,
    L: int = 50,
    theta_S: float = 0.15,
    theta_R: int = 5,
    k1: int = 6,
    k2: int = 6,
    basal_end: Optional[int] = None,
    theta_delta_S: Optional[float] = None,
    basal_mode: str = "mean_delta",
    basal_q: float = 90.0,
    mad_kappa: float = 2.5,
    hop: int = 1,
) -> Dict[str, np.ndarray]:
    """
    Ordinal Persistence of Synergy surplus (Proposal B).

    Basal modes (when basal_end is set)
    -----------------------------------
    mean_delta (default / legacy I0):
        ell_t = 1{ S_t - mean(S_basal) >= theta_delta_S }
    percentile:
        ell_t = 1{ S_t >= percentile(S_basal, basal_q) }
    mad:
        ell_t = 1{ S_t >= median(S_basal) + mad_kappa * MAD }

    Without basal_end: absolute ell_t = 1{ S_t >= theta_S }.

    Persistence
    -----------
    hop=1 (default): consecutive active samples for theta_R steps (legacy).
    hop>1: run-length credit only on semi-independent indices
        t where (t - basal_end) % hop == 0 (or every hop from t=0 if no basal).
        theta_R is counted in *hop credits*, not raw sample steps.
    """
    if theta_R < 1:
        raise ValueError("theta_R must be >= 1")
    if hop < 1:
        raise ValueError("hop must be >= 1")
    S = ordinal_synergy_surplus(pi1, pi2, L=L, k1=k1, k2=k2)
    T = S.size
    alarm = np.zeros(T, dtype=np.int8)
    high = np.zeros(T, dtype=np.int8)
    persistence = np.zeros(T, dtype=np.int64)

    S_basal = np.nan
    thr = float(theta_S)
    use_absolute_thr = True
    mode_resolved = "absolute"
    if basal_end is not None:
        if basal_end < L:
            raise ValueError("basal_end must be >= L")
        basal_vals = S[L - 1 : basal_end]
        S_basal, thr, mode_resolved = compute_surplus_threshold(
            basal_vals,
            basal_mode=basal_mode,
            theta_delta_S=theta_delta_S,
            basal_q=basal_q,
            mad_kappa=mad_kappa,
        )
        # mean_delta keeps legacy (S - mean) >= delta form for bit-identical tests;
        # percentile/mad use absolute thr on S.
        use_absolute_thr = mode_resolved != "mean_delta"

    origin = int(basal_end) if basal_end is not None else 0
    run = 0
    for t in range(T):
        if not np.isfinite(S[t]):
            run = 0
            continue
        if basal_end is not None and t < basal_end:
            run = 0
            continue

        if use_absolute_thr:
            active = S[t] >= thr
        else:
            active = (S[t] - S_basal) >= float(theta_delta_S)

        if active:
            high[t] = 1

        on_hop = hop == 1 or ((t - origin) % hop == 0)
        if on_hop:
            if active:
                run += 1
            else:
                run = 0
        persistence[t] = run
        if high[t] == 1 and run >= theta_R and on_hop:
            alarm[t] = 1

    out: Dict[str, np.ndarray] = {
        "alarm": alarm,
        "synergy": S,
        "high_syn": high,
        "persistence": persistence,
        "S_thr": np.array([float(thr)]),
        "hop": np.array([int(hop)]),
    }
    if basal_end is not None:
        out["S_basal"] = np.array([S_basal])
        out["basal_mode"] = np.array([mode_resolved], dtype=object)
    return out


def joint_symbols_from_factors(pi1: np.ndarray, pi2: np.ndarray, k2: int = 6) -> np.ndarray:
    """Encode bivariate ordinal factors as joint codes in {0,...,k1*k2-1}."""
    pi1 = np.asarray(pi1, dtype=np.int64).ravel()
    pi2 = np.asarray(pi2, dtype=np.int64).ravel()
    return pi1 * int(k2) + pi2


def _windowed_diversity(sigma: np.ndarray, *, L: int, K: int) -> np.ndarray:
    """Support diversity D_t = |supp(W_t)| / K on joint symbols (ordinal only)."""
    sigma = np.asarray(sigma, dtype=np.int64).ravel()
    T = sigma.size
    D = np.full(T, np.nan, dtype=float)
    for t in range(L - 1, T):
        window = sigma[t - L + 1 : t + 1]
        D[t] = len(np.unique(window)) / float(K)
    return D


def opsp_integrated_detect(
    pi1: np.ndarray,
    pi2: np.ndarray,
    *,
    L: int = 50,
    theta_R: int = 5,
    k1: int = 6,
    k2: int = 6,
    basal_end: Optional[int] = None,
    theta_delta_S: Optional[float] = None,
    theta_S: float = 0.15,
    collapse_role: str = "none",
    theta_D: float = 0.35,
    confirm_window: int = 5,
    modulate_delta_R: int = 2,
    modulate_delta_S: float = 0.02,
    basal_mode: str = "mean_delta",
    basal_q: float = 90.0,
    mad_kappa: float = 2.5,
    hop: int = 1,
) -> Dict[str, np.ndarray]:
    """
    Integrated ordinal detector: Nivel-2 *persistence* on Nivel-3 *synergy surplus*,
    with collapse/locking only complementary (never the sole alarm source).

    Primary predicate (all roles)
    ----------------------------
    High synergistic surplus indicator ell^S_t (absolute S_t, mean+ΔS, percentile,
    or MAD basal), then persistence for theta_R hop-credits (Φ₂ run-length logic on
    an OPS / Level-3–consistent proxy — not on diversity collapse).

    Structural options (Track S1)
    -----------------------------
    basal_mode: mean_delta | percentile | mad  (see compute_surplus_threshold)
    hop:        hop=1 legacy consecutive samples; hop>1 semi-independent credits

    Collapse roles (subordinate)
    ----------------------------
    - ``none``: pure surplus-persistence core (I0). Collapse ignored for alarm.
    - ``tag``: same alarm as ``none``; also flags collapse co-occurrence near alarms.
    - ``modulate``: when D_t ≤ theta_D, slightly relax surplus threshold and/or
      required run length — collapse *helps* but cannot create an alarm without
      surplus mass. (modulate still uses mean_delta-style thr shifts when applicable)
    - ``confirm``: core surplus-persistence alarm is *kept* only if a low-div
      window exists within ±confirm_window samples. Collapse filters; never
      invents alarms. May silence pure high-support high-synergy paths — use
      only as a FAR-leaning secondary variant.

    Constraints
    -----------
    - Pure ordinal (joint TV surplus + support diversity). No μ/σ or abs-z.
    - OPS surplus is a Level-3–consistent *proxy*, not continuous excess3.
    - Collapse-only streams must not alarm under any role.
    """
    role = str(collapse_role).lower().strip()
    if role not in ("none", "tag", "modulate", "confirm"):
        raise ValueError(
            "collapse_role must be one of: none, tag, modulate, confirm"
        )
    if theta_R < 1:
        raise ValueError("theta_R must be >= 1")
    if hop < 1:
        raise ValueError("hop must be >= 1")
    if not (0.0 < theta_D <= 1.0):
        raise ValueError("theta_D must be in (0, 1]")
    if confirm_window < 0:
        raise ValueError("confirm_window must be >= 0")
    if modulate_delta_R < 0:
        raise ValueError("modulate_delta_R must be >= 0")
    if modulate_delta_S < 0:
        raise ValueError("modulate_delta_S must be >= 0")

    S = ordinal_synergy_surplus(pi1, pi2, L=L, k1=k1, k2=k2)
    joint = joint_symbols_from_factors(pi1, pi2, k2=k2)
    K = int(k1) * int(k2)
    diversity = _windowed_diversity(joint, L=L, K=K)
    T = S.size

    low_div = np.zeros(T, dtype=np.int8)
    for t in range(T):
        if np.isfinite(diversity[t]) and diversity[t] <= theta_D:
            low_div[t] = 1

    S_basal = np.nan
    thr_abs = float(theta_S)
    use_relative = basal_end is not None
    use_mean_delta = False
    mode_resolved = "absolute"
    if use_relative:
        if basal_end < L:
            raise ValueError("basal_end must be >= L")
        basal_vals = S[L - 1 : basal_end]
        S_basal, thr_abs, mode_resolved = compute_surplus_threshold(
            basal_vals,
            basal_mode=basal_mode,
            theta_delta_S=theta_delta_S,
            basal_q=basal_q,
            mad_kappa=mad_kappa,
        )
        use_mean_delta = mode_resolved == "mean_delta"

    origin = int(basal_end) if basal_end is not None else 0
    high_syn = np.zeros(T, dtype=np.int8)
    persistence = np.zeros(T, dtype=np.int64)
    core_alarm = np.zeros(T, dtype=np.int8)
    effective_theta_R = np.full(T, theta_R, dtype=np.int64)

    run = 0
    for t in range(T):
        if not np.isfinite(S[t]):
            run = 0
            continue
        if use_relative and t < basal_end:
            run = 0
            continue

        need_R = int(theta_R)
        on_hop = hop == 1 or ((t - origin) % hop == 0)

        if use_mean_delta:
            # Legacy mean+ΔS; modulate may soften the *delta* gate
            base_thr = float(theta_delta_S)  # type: ignore[arg-type]
            thr_delta = base_thr
            if role == "modulate" and low_div[t] == 1:
                if base_thr > 0.0:
                    thr_delta = max(0.5 * base_thr, base_thr - float(modulate_delta_S))
                need_R = max(1, int(theta_R) - int(modulate_delta_R))
            active = (S[t] - S_basal) >= thr_delta
        elif use_relative:
            # percentile / mad: absolute thr on S; modulate can lower thr slightly
            thr_use = float(thr_abs)
            if role == "modulate" and low_div[t] == 1:
                thr_use = thr_use - float(modulate_delta_S)
                need_R = max(1, int(theta_R) - int(modulate_delta_R))
            active = S[t] >= thr_use
        else:
            thr_use = float(theta_S)
            if role == "modulate" and low_div[t] == 1:
                thr_use = max(0.0, thr_use - float(modulate_delta_S))
                need_R = max(1, int(theta_R) - int(modulate_delta_R))
            active = S[t] >= thr_use

        effective_theta_R[t] = need_R

        if active:
            high_syn[t] = 1

        if on_hop:
            if active:
                run += 1
            else:
                run = 0
        persistence[t] = run
        if high_syn[t] == 1 and run >= need_R and on_hop:
            core_alarm[t] = 1

    # Collapse never creates alarms alone: start from core surplus-persistence
    alarm = core_alarm.copy()
    collapse_coincident = np.zeros(T, dtype=np.int8)

    if role in ("tag", "confirm"):
        for t in range(T):
            if core_alarm[t] != 1:
                continue
            lo = max(0, t - confirm_window)
            hi = min(T, t + confirm_window + 1)
            if int(low_div[lo:hi].sum()) > 0:
                collapse_coincident[t] = 1
        if role == "confirm":
            # Keep only surplus-core alarms that co-occur with collapse nearby
            alarm = (core_alarm == 1) & (collapse_coincident == 1)
            alarm = alarm.astype(np.int8)
    elif role == "modulate":
        # Co-occurrence flag for diagnostics (alarm already may use relaxed gates)
        for t in range(T):
            if alarm[t] == 1 and low_div[t] == 1:
                collapse_coincident[t] = 1
    # role == "none": collapse_coincident stays zero

    out: Dict[str, np.ndarray] = {
        "alarm": alarm,
        "core_alarm": core_alarm,
        "synergy": S,
        "high_syn": high_syn,
        "persistence": persistence,
        "diversity": diversity,
        "low_div": low_div,
        "collapse_coincident": collapse_coincident,
        "effective_theta_R": effective_theta_R,
        "S_thr": np.array([float(thr_abs if use_relative else theta_S)]),
        "hop": np.array([int(hop)]),
    }
    if use_relative:
        out["S_basal"] = np.array([S_basal])
        out["basal_mode"] = np.array([mode_resolved], dtype=object)
    return out
