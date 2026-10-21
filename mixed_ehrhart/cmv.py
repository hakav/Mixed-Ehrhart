"""Combinatorial mixed volumes and exact inequality checks."""

from itertools import product
from operator import index

from sage.all import ZZ

from .ehrhart import forward_diff_coeffs
from .geometry import (
    _MinkowskiCache,
    _group_polytopes,
    _lattice_count_validated,
    validate_polytopes,
)

__all__ = [
    "combinatorial_mixed_volume",
    "precompute_counts_bivariate",
    "cmv_bivariate",
    "check_af_inequality",
    "degree3_slice_values",
    "check_degree2_inequality",
    "check_degree3_inequalities",
]


def _nonnegative_integer(value, name):
    """Accept integer indices without rounding."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a nonnegative integer.")
    try:
        result = index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be a nonnegative integer.") from exc
    if result < 0:
        raise ValueError(f"{name} must be nonnegative.")
    return result


def _group_scales(scales, groups, size):
    result = [0] * size
    for scale, group in zip(scales, groups):
        result[group] += scale
    return tuple(result)


def _counts_for_scales(polys, scales, count):
    """Count each distinct dilation sum once; inputs are already validated."""
    counter = _lattice_count_validated if count is None else count
    if not callable(counter):
        raise TypeError("count must be a callable or None.")
    unique, groups = _group_polytopes(polys)
    sums = _MinkowskiCache(unique)
    cached = {(0,) * len(unique): ZZ(1)}
    counts = {}
    for mult in scales:
        if mult in counts:
            continue
        key = _group_scales(mult, groups, len(unique))
        if key not in cached:
            value = counter(sums.sum(key))
            cached[key] = ZZ(_nonnegative_integer(value, "lattice count"))
        counts[mult] = cached[key]
    return counts


def combinatorial_mixed_volume(polys, *, count=None):
    """Return the CMV, grouping repeated polytopes to reuse lattice counts.

    An empty family has CMV one; an optional exact counter receives each sum.
    """
    polys = tuple(polys)
    if not polys:
        if count is not None and not callable(count):
            raise TypeError("count must be a callable or None.")
        return ZZ(1)
    polys = validate_polytopes(polys)
    if count is not None and not callable(count):
        raise TypeError("count must be a callable or None.")
    unique, groups = _group_polytopes(polys)
    if len(polys) > _MinkowskiCache(unique).dimension():
        return ZZ(0)
    order = _group_scales((1,) * len(polys), groups, len(unique))
    scales = product(*(range(m + 1) for m in order))
    counts = _counts_for_scales(unique, scales, count)
    return forward_diff_coeffs(counts, order)


def precompute_counts_bivariate(P, Q, max_a=3, max_b=3, *, count=None):
    """Return counts of ``a*P + b*Q`` on the requested rectangular grid.

    Bounds are inclusive nonnegative integers; the origin count is one.
    """
    polys = validate_polytopes((P, Q))
    max_a = _nonnegative_integer(max_a, "max_a")
    max_b = _nonnegative_integer(max_b, "max_b")
    scales = product(range(max_a + 1), range(max_b + 1))
    return _counts_for_scales(polys, scales, count)


def cmv_bivariate(counts, mP, mQ):
    """Return the CMV of ``mP`` copies of P and ``mQ`` copies of Q.

    The count table must include the whole rectangle up to ``(mP,mQ)``.
    """
    mP = _nonnegative_integer(mP, "mP")
    mQ = _nonnegative_integer(mQ, "mQ")
    if counts[(0, 0)] != 1:
        raise ValueError("counts[(0, 0)] must be one for the singleton origin.")
    return forward_diff_coeffs(counts, (mP, mQ))


def check_af_inequality(P, Q, others=(), *, count=None):
    """Check ``CMV(P,Q,*others)^2 >= CMV(P,P,*others)*CMV(Q,Q,*others)``.

    Return the three CMVs, both sides, and ``holds`` for these inputs.
    """
    polys = validate_polytopes((P, Q) + tuple(others))
    if count is not None and not callable(count):
        raise TypeError("count must be a callable or None.")
    unique, groups = _group_polytopes(polys)
    tail = (1,) * (len(polys) - 2)
    orders = tuple(
        _group_scales(order, groups, len(unique))
        for order in ((1, 1) + tail, (2, 0) + tail, (0, 2) + tail)
    )
    if len(polys) > _MinkowskiCache(unique).dimension():
        mixed = repeated_P = repeated_Q = ZZ(0)
    else:
        needed = set()
        for order in set(orders):
            needed.update(product(*(range(m + 1) for m in order)))
        counts = _counts_for_scales(unique, sorted(needed), count)
        values = {order: forward_diff_coeffs(counts, order) for order in set(orders)}
        mixed, repeated_P, repeated_Q = (values[order] for order in orders)
    lhs = mixed**2
    rhs = repeated_P * repeated_Q
    return {
        "mixed": mixed,
        "repeated_P": repeated_P,
        "repeated_Q": repeated_Q,
        "lhs": lhs,
        "rhs": rhs,
        "holds": bool(lhs >= rhs),
    }


def check_degree2_inequality(P, Q, *, count=None):
    """Return the CMVs, both sides, and result of the degree-two comparison."""
    result = check_af_inequality(P, Q, count=count)
    return {
        "CMV_PP": result["repeated_P"],
        "CMV_PQ": result["mixed"],
        "CMV_QQ": result["repeated_Q"],
        "lhs": result["lhs"],
        "rhs": result["rhs"],
        "holds": result["holds"],
    }


def degree3_slice_values(P, Q, *, count=None):
    """Return CMVs for multiplicities ``(3,0), (2,1), (1,2), (0,3)``."""
    polys = validate_polytopes((P, Q))
    scales = ((a, b) for a in range(4) for b in range(4 - a))
    counts = _counts_for_scales(polys, scales, count)
    return tuple(cmv_bivariate(counts, 3 - b, b) for b in range(4))


def check_degree3_inequalities(P, Q, *, count=None):
    """Check ``A*C <= B^2`` and ``B*D <= C^2`` for the degree-three slice.

    Return the four CMVs, the two comparisons, and their exact sides.
    """
    A, B, C, D = degree3_slice_values(P, Q, count=count)
    ineq1_lhs, ineq1_rhs = A * C, B**2
    ineq2_lhs, ineq2_rhs = B * D, C**2
    return {
        "CMV_P3": A,
        "CMV_P2Q": B,
        "CMV_PQ2": C,
        "CMV_Q3": D,
        "ineq1_holds": bool(ineq1_lhs <= ineq1_rhs),
        "ineq2_holds": bool(ineq2_lhs <= ineq2_rhs),
        "ineq1_lhs": ineq1_lhs,
        "ineq1_rhs": ineq1_rhs,
        "ineq2_lhs": ineq2_lhs,
        "ineq2_rhs": ineq2_rhs,
    }
