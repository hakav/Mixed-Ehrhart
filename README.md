# Mixed Ehrhart

SageMath code for exact multivariate Ehrhart polynomials, combinatorial mixed volumes (CMVs), and Lorentzian slice checks. 

## Install and run

Requires SageMath 10.6 and Python 3.11 or later. Run with Conda by

```sh
conda env create -f environment.yaml
conda activate mixed-ehrhart
python -m pip install --no-build-isolation -e .
sage -n jupyterlab
```

Open a notebook in `notebooks/` and select its SageMath kernel. Each runs from the first cell. Installing this package with pip does not install SageMath.

## Background

For lattice polytopes $P_1,\ldots,P_k\subseteq\mathbb R^d$, the library computes

$$
E(n)=\#\bigl((n_1P_1+\cdots+n_kP_k)\cap\mathbb Z^d\bigr)
=\sum_{|\alpha|\le D}c_\alpha\prod_i\binom{n_i}{\alpha_i},
\qquad c_\alpha=\Delta^\alpha E(0),
$$

where $D=\dim(P_1+\cdots+P_k)$. `multivariate_ehrhart` returns both the ordinary polynomial over `QQ` and the dictionary of integer forward differences $c_\alpha$. These are coefficients of the multivariate Ehrhart polynomial in the binomial basis. 

Combinatorial mixed volume is the inclusion–exclusion sum

$$
\operatorname{CMV}(P_1,\ldots,P_k)
=\sum_{I\subseteq\{1,\ldots,k\}}(-1)^{k-|I|}
\#\bigl((\sum_{i\in I}P_i)\cap\mathbb Z^d\bigr).
$$

The empty sum is the origin, with count 1. Repeated inputs give higher forward differences. The bivariate helpers reuse counts for repeated copies of two polytopes.

`polynomial_slice` forms $\sum_{|\alpha|=m}c_\alpha y^\alpha/\alpha!$. This normalized slice agrees with the ordinary homogeneous part in top degree, but can differ in lower degrees. `lorentzian_slice_test` checks nonnegative coefficients, M-convex support, and, for degree at least two, at most one positive eigenvalue for the Hessian of every derivative of order $m-2$. Signatures are computed exactly. The characterization is [Brändén–Huh, Theorem 2.25](https://arxiv.org/abs/1902.03719).

## Examples

- [Ehrhart polynomials](notebooks/01_ehrhart_polynomials.ipynb): perpendicular segments, repeated triangles, and a Reeve tetrahedron.
- [Mixed volumes and slices](notebooks/02_mixed_volumes_and_slices.ipynb): inclusion–exclusion, repeated inputs, inequalities, and Lorentzian slice checks.

