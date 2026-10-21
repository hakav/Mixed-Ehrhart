"""Lattice polytopes, Minkowski sums, and exact point counts."""

from collections import OrderedDict
from fractions import Fraction
from functools import reduce
from math import gcd, prod
from operator import index

from sage.all import Polyhedron, QQ, ZZ, matrix, parent
from sage.features import FeatureNotPresentError
from sage.geometry.polyhedron.base import Polyhedron_base


def _integer_coordinate(value):
    """Read an integer coordinate without rounding."""
    ring = parent(value)
    if not isinstance(value, (int, Fraction)) and not getattr(ring, "is_exact", lambda: False)():
        raise ValueError("Coordinates must be integers")
    try:
        rational = QQ(value)
        if rational.denominator() != 1:
            raise ValueError("Coordinates must be integers")
        return ZZ(rational)
    except (TypeError, ValueError) as error:
        raise ValueError("Coordinates must be integers") from error


def _validate_polytope(poly, *, allow_empty=False):
    if not isinstance(poly, Polyhedron_base):
        raise TypeError("Expected a SageMath Polyhedron.")
    if poly.is_empty():
        if allow_empty:
            return poly
        raise ValueError("Polytopes must be nonempty.")
    if not poly.is_compact():
        raise ValueError("Polytopes must be bounded.")
    if not poly.base_ring().is_exact():
        raise ValueError("Polytopes must use exact coordinates")
    for vertex in poly.vertices():
        for coordinate in vertex:
            _integer_coordinate(coordinate)
    return poly


def validate_polytopes(polys):
    """Check that the inputs are nonempty lattice polytopes in one space."""
    polys = tuple(polys)
    if not polys:
        raise ValueError("Provide at least one polytope")
    for poly in polys:
        _validate_polytope(poly)
    ambient_dimension = polys[0].ambient_dim()
    if any(poly.ambient_dim() != ambient_dimension for poly in polys):
        raise ValueError("Polytopes must have the same ambient dimension")
    return polys


def poly_from_vertices(vertices):
    """Take the convex hull of a nonempty list of integer vertices."""
    rows = [tuple(_integer_coordinate(c) for c in vertex) for vertex in vertices]
    if not rows:
        raise ValueError("Provide at least one vertex")
    ambient_dimension = len(rows[0])
    if any(len(row) != ambient_dimension for row in rows):
        raise ValueError("All vertices must have the same number of coordinates")
    return Polyhedron(vertices=rows, base_ring=QQ)


def minkowski_sum_scaled(polys, scales):
    """Form a Minkowski sum with nonnegative integer scales"""
    polys = validate_polytopes(polys)
    scales = tuple(scales)
    if len(scales) != len(polys):
        raise ValueError("Provide exactly one scale per polytope")
    normalized = []
    for scale in scales:
        try:
            if isinstance(scale, bool):
                raise TypeError
            scale = index(scale)
        except TypeError as error:
            raise ValueError("Scales must be nonnegative integers") from error
        if scale < 0:
            raise ValueError("Scales must be nonnegative integers")
        normalized.append(scale)
    return _MinkowskiCache(polys).sum(normalized)


def _group_polytopes(polys):
    """Group equal polytopes, keeping the position of every original input"""
    unique, groups = [], []
    identities, vertices = {}, {}
    for poly in polys:
        group = identities.get(id(poly))
        if group is None:
            key = frozenset(tuple(vertex) for vertex in poly.vertices())
            group = vertices.get(key)
            if group is None:
                group = len(unique)
                unique.append(poly)
                vertices[key] = group
            identities[id(poly)] = group
        groups.append(group)
    return tuple(unique), tuple(groups)


class _MinkowskiCache:
    """Reuse dilations and a bounded number of partial sums within one call."""

    def __init__(self, polys):
        self.polys = polys
        self._dilations = {}
        self._prefixes = OrderedDict()
        self._origin = None

    def dimension(self):
        ambient = self.polys[0].ambient_dim()
        if any(poly.dimension() == ambient for poly in self.polys):
            return int(ambient)
        directions = []
        for poly in self.polys:
            vertices = poly.vertices_list()
            anchor = vertices[0]
            directions.extend([[b - a for a, b in zip(anchor, vertex)]
                               for vertex in vertices[1:]])
        return int(matrix(QQ, directions).rank()) if directions else 0

    def _scaled(self, i, scale):
        if scale == 1:
            return self.polys[i]
        key = (i, scale)
        if key not in self._dilations:
            self._dilations[key] = self.polys[i].dilation(scale)
        return self._dilations[key]

    def sum(self, scales):
        scales = tuple(scales)
        active = [i for i, scale in enumerate(scales) if scale]
        if not active:
            if self._origin is None:
                self._origin = Polyhedron(
                    vertices=[(0,) * self.polys[0].ambient_dim()], base_ring=QQ
                )
            return self._origin
        result = None
        for i in active:
            key = scales[:i + 1]
            if key in self._prefixes:
                result = self._prefixes[key]
                self._prefixes.move_to_end(key)
                continue
            part = self._scaled(i, scales[i])
            result = part if result is None else result.minkowski_sum(part)
            # Full sums are counted once; only partial sums are worth keeping.
            if i < len(self.polys) - 1:
                self._prefixes[key] = result
                if len(self._prefixes) > 128:
                    self._prefixes.popitem(last=False)
        return result


def lattice_count(poly, *, max_enumeration_box=100_000):
    """Count lattice points, with limit on enumeration."""
    try:
        if isinstance(max_enumeration_box, bool):
            raise TypeError
        limit = index(max_enumeration_box)
    except TypeError as error:
        raise ValueError("max_enumeration_box must be a nonnegative integer.") from error
    if limit < 0:
        raise ValueError("max_enumeration_box must be a nonnegative integer.")
    _validate_polytope(poly, allow_empty=True)
    return _lattice_count_validated(poly, max_enumeration_box=limit)


def _lattice_count_validated(poly, *, max_enumeration_box=100_000):
    """Count a polytope whose integer vertices have already been checked."""
    if poly.is_empty():
        return ZZ(0)
    dimension = poly.dimension()
    if dimension == 0:
        return ZZ(1)
    if dimension == 1:
        first, last = poly.vertices_list()
        return ZZ(reduce(gcd, (abs(ZZ(a) - ZZ(b)) for a, b in zip(first, last)), 0) + 1)
    if poly.ambient_dim() == 2:
        boundary = ZZ(0)
        for edge in poly.faces(1):
            first, last = edge.vertices()
            boundary += gcd(abs(ZZ(first[0]) - ZZ(last[0])),
                            abs(ZZ(first[1]) - ZZ(last[1])))
        return ZZ(poly.volume() + QQ(boundary) / 2 + 1)
    lower, upper = poly.bounding_box(integral_hull=True)
    box_size = prod(ZZ(b) - ZZ(a) + 1 for a, b in zip(lower, upper))
    if box_size <= max_enumeration_box:
        return ZZ(len(poly.integral_points()))
    # The QQ counter respects the request to avoid enumeration.
    rational_poly = Polyhedron(vertices=poly.vertices_list(), base_ring=QQ)
    try:
        return ZZ(rational_poly.integral_points_count(explicit_enumeration_threshold=0))
    except (FeatureNotPresentError, NotImplementedError) as error:
        raise RuntimeError(
            f"The integer bounding box has {box_size} points, exceeding the "
            f"enumeration limit {max_enumeration_box}. Install Sage's optional latte_int "
            "package, supply count= with your own exact counter, or deliberately "
            "raise max_enumeration_box using functools.partial(lattice_count, ...)."
        ) from error
