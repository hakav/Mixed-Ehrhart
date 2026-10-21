from fractions import Fraction
from math import factorial, prod
from operator import index

from sage.all import AA, PolynomialRing, QQ, ZZ, parent, zero_matrix
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.polynomial.polynomial_element import Polynomial


def _nonnegative_integer(value, label):
    try:
        if isinstance(value, bool):
            raise TypeError
        value = index(value)
    except TypeError as error:
        raise ValueError(f"{label} must be a nonnegative integer.") from error
    if value < 0:
        raise ValueError(f"{label} must be a nonnegative integer.")
    return value


def _support_vectors(support):
    vectors = {tuple(_nonnegative_integer(x, "Exponents") for x in row) for row in support}
    if vectors:
        k = len(next(iter(vectors)))
        if any(len(row) != k for row in vectors):
            raise ValueError("Exponent vectors must have the same length.")
    return vectors


def is_m_convex_support(support):
    # Check simultaneous exchange for nonnegative integer exponent vectors.
    support = _support_vectors(support)
    if not support:
        return True
    if len({sum(row) for row in support}) != 1:
        return False
    for u in support:
        for v in support:
            for i in range(len(u)):
                if u[i] <= v[i]:
                    continue
                found = False
                for j in range(len(u)):
                    if u[j] >= v[j]:
                        continue
                    left, right = list(u), list(v)
                    left[i] -= 1
                    left[j] += 1
                    right[i] += 1
                    right[j] -= 1
                    if tuple(left) in support and tuple(right) in support:
                        found = True
                        break
                if not found:
                    return False
    return True


def polynomial_slice(coefficients, degree, *, names=None, normalize=True):
    """Build a homogeneous slice from the Ehrhart forward differences.

    Divide by exponent factorials unless ``normalize=False``.
    """
    degree = _nonnegative_integer(degree, "degree")
    if not coefficients:
        raise ValueError("Empty coefficient dictionary")
    support = _support_vectors(coefficients)
    k = len(next(iter(support)))
    if names is None:
        names = [f"y{i + 1}" for i in range(k)]
    elif isinstance(names, str):
        names = [name.strip() for name in names.split(",")]
    else:
        names = list(names)
    if len(names) != k or len(set(names)) != k:
        raise ValueError("Provide one distinct variable name per coordinate")
    ring = PolynomialRing(QQ, k, names=names)
    variables = ring.gens()
    result = ring.zero()
    for alpha, value in coefficients.items():
        if not isinstance(value, (int, Fraction)) and not getattr(parent(value), "is_exact", lambda: False)():
            raise ValueError("Coefficients must be rational numbers")
        value = QQ(value)
        if sum(alpha) != degree or not value:
            continue
        denominator = prod(factorial(a) for a in alpha) if normalize else 1
        monomial = prod(y ** a for y, a in zip(variables, alpha))
        result += value / denominator * monomial
    return result


def _compositions(total, parts):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, parts - 1):
            yield (first,) + rest


def lorentzian_slice_test(polynomial):
    """Return ``(holds, failures)`` for a homogeneous QQ/ZZ polynomial.

    Check coefficients, M-convex support, and exact Hessian signatures.
    """
    if not isinstance(polynomial, (Polynomial, MPolynomial)):
        raise TypeError("Expected a Sage polynomial over QQ or ZZ.")
    original_ring = polynomial.parent()
    if original_ring.base_ring() not in (QQ, ZZ):
        raise ValueError("Use exact rational or integer polynomial coefficients.")
    ring = PolynomialRing(QQ, original_ring.ngens(), names=original_ring.variable_names())
    polynomial = ring(polynomial)
    if polynomial.is_zero():
        return True, []
    degrees = {sum(alpha) for alpha in polynomial.dict()}
    if len(degrees) != 1:
        raise ValueError("The polynomial must be homogeneous")
    degree = degrees.pop()
    failures = []
    if any(value < 0 for value in polynomial.coefficients()):
        failures.append({"condition": "nonnegative_coefficients"})
    if not is_m_convex_support(polynomial.dict()):
        failures.append({"condition": "m_convex_support"})
    if failures or degree <= 1:
        return not failures, failures
    variables = ring.gens()
    k = len(variables)
    for beta in _compositions(degree - 2, k):
        quadratic = polynomial
        for variable, order in zip(variables, beta):
            for _ in range(order):
                quadratic = quadratic.derivative(variable)
        if quadratic.is_zero():
            continue
        hessian = zero_matrix(QQ, k, k)
        for i, first in enumerate(variables):
            for j, second in enumerate(variables):
                hessian[i, j] = QQ(quadratic.derivative(first).derivative(second))
        eigenvalues = hessian.charpoly().roots(ring=AA)
        positive = sum(multiplicity for value, multiplicity in eigenvalues if value > 0)
        if positive > 1:
            failures.append({"condition": "hessian_signature", "derivative": beta,
                             "positive_eigenvalues": positive, "eigenvalues": eigenvalues})
    return not failures, failures
