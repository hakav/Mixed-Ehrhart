from itertools import product
from math import comb
from operator import index

from sage.all import PolynomialRing, QQ, ZZ

from .geometry import (
    _MinkowskiCache,
    _group_polytopes,
    _lattice_count_validated,
    poly_from_vertices,
    validate_polytopes,
)


def _nonnegative_integer(value, label):
    """Read an integer index without rounding."""
    try:
        result = index(value)
    except TypeError as exc:
        raise TypeError(f"{label} must be a nonnegative integer") from exc
    if result < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return result


def multiindices_leq_degree(k, degree):
    """List nonnegative k-tuples by total degree, then lexicographically"""
    k = _nonnegative_integer(k, "k")
    degree = _nonnegative_integer(degree, "degree")

    def fixed_degree(length, total):
        if length == 0:
            if total == 0:
                yield ()
            return
        if length == 1:
            yield (total,)
            return
        for first in range(total + 1):
            for tail in fixed_degree(length - 1, total - first):
                yield (first,) + tail

    return [
        alpha
        for total in range(degree + 1)
        for alpha in fixed_degree(k, total)
    ]


def forward_diff_coeffs(counts, alpha):
    """Return Delta^alpha E(0)"""
    alpha = tuple(_nonnegative_integer(a, "index entry") for a in alpha)
    weights = [
        [(-1) ** (a - b) * comb(a, b) for b in range(a + 1)]
        for a in alpha
    ]
    total = ZZ(0)
    for beta in product(*(range(a + 1) for a in alpha)):
        multiplier = 1
        for row, b in zip(weights, beta):
            multiplier *= row[b]
        total += multiplier * ZZ(counts[beta])
    return total


def _forward_differences(counts, indices, degree):
    """Apply ordinary difference tables along each axis of the simplex"""
    coefficients = dict(counts)
    for axis in range(len(indices[0])):
        for alpha in indices:
            if alpha[axis]:
                continue
            length = degree - sum(alpha) + 1
            keys = [alpha[:axis] + (j,) + alpha[axis + 1:] for j in range(length)]
            values = [coefficients[key] for key in keys]
            for order in range(1, length):
                for j in range(length - 1, order - 1, -1):
                    values[j] -= values[j - 1]
            coefficients.update(zip(keys, values))
    return coefficients


def _checked_count(counter, poly):
    """Check that a counter returned an exact nonnegative integer"""
    value = counter(poly)
    if isinstance(value, (float, complex)):
        raise ValueError("The counting function must return an exact integer.")
    parent = getattr(value, "parent", None)
    if callable(parent):
        is_exact = getattr(parent(), "is_exact", None)
        if callable(is_exact) and not is_exact():
            raise ValueError("The counting function must return an integer")
    try:
        result = ZZ(value)
        valid = value == result and result >= 0
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(
            "The counting function must return a nonnegative integer"
        ) from exc
    if not valid:
        raise ValueError(
            "The counting function must return a nonnegative integer"
        )
    return result


def _variable_names(names, k):
    if names is None:
        return tuple(f"n{i + 1}" for i in range(k))
    if isinstance(names, str):
        names = tuple(part.strip() for part in names.split(","))
    else:
        names = tuple(names)
    if len(names) != k:
        raise ValueError(f"Expected {k} variable names, received {len(names)}")
    if any(not isinstance(name, str) or not name.isidentifier() for name in names):
        raise ValueError("Variable names must be nonempty Python identifiers")
    if len(set(names)) != k:
        raise ValueError("Variable names must be distinct")
    return names


def multivariate_ehrhart(polys, *, names=None, count=None):
    """Return the polynomial over QQ and its binomial-basis coefficient dict.

    An optional count(poly) supplies exact lattice counts
    """
    polys = validate_polytopes(polys)
    k = len(polys)
    variable_names = _variable_names(names, k)
    if count is not None and not callable(count):
        raise TypeError("count must be callable")

    unique, groups = _group_polytopes(polys)
    sums = _MinkowskiCache(unique)
    degree = sums.dimension()
    indices = multiindices_leq_degree(k, degree)
    grouped_indices = (
        indices if len(unique) == k else multiindices_leq_degree(len(unique), degree)
    )
    counter = (
        _lattice_count_validated
        if count is None else lambda poly: _checked_count(count, poly)
    )
    counts = {alpha: counter(sums.sum(alpha)) for alpha in grouped_indices}
    grouped_coefficients = _forward_differences(counts, grouped_indices, degree)
    if len(unique) == k:
        coefficients = grouped_coefficients
    else:
        coefficients = {}
        for alpha in indices:
            totals = [0] * len(unique)
            for group, scale in zip(groups, alpha):
                totals[group] += scale
            coefficients[alpha] = grouped_coefficients[tuple(totals)]

    # Keep exponent tuples consistent when there is only one variable.
    ring = PolynomialRing(QQ, k, names=variable_names)
    variables = ring.gens()
    grouped_variables = [ring.zero() for _ in unique]
    for variable, group in zip(variables, groups):
        grouped_variables[group] += variable
    binomial_polynomials = []
    for variable in grouped_variables:
        basis = [ring.one()]
        for order in range(1, degree + 1):
            basis.append(basis[-1] * (variable - order + 1) / order)
        binomial_polynomials.append(basis)

    polynomial = ring.zero()
    for alpha, coefficient in grouped_coefficients.items():
        if coefficient == 0:
            continue
        term = ring(coefficient)
        for i, order in enumerate(alpha):
            term *= binomial_polynomials[i][order]
        polynomial += term
    return polynomial, coefficients


def multivariate_ehrhart_from_vertices(vertex_lists, *, names=None, count=None):
    """Compute the polynomial from one vertex list per polytope."""
    polys = [poly_from_vertices(vertices) for vertices in vertex_lists]
    return multivariate_ehrhart(polys, names=names, count=count)


def homogeneous_part(polynomial, degree):
    """Return terms of the requested degree, keeping the original ring."""
    degree = _nonnegative_integer(degree, "degree")
    ring = polynomial.parent()
    result = ring.zero()
    for exponent, coefficient in polynomial.dict().items():
        # Sage uses integer keys for ordinary univariate polynomials and
        # exponent tuples for multivariate polynomials (including k == 1).
        try:
            powers = tuple(exponent)
        except TypeError:
            powers = (exponent,)
        if sum(powers) == degree:
            term = ring(coefficient)
            for variable, power in zip(ring.gens(), powers):
                term *= variable ** power
            result += term
    return result
