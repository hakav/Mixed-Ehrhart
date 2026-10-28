# Mixed Ehrhart

SageMath code for computing multivariate Ehrhart polynomials, combinatorial mixed volumes (CMVs), and checking Lorentzianity properties in this context. 

## Install and run

Requires SageMath 10.6 and Python 3.11 or later. Run with Conda by

```sh
conda env create -f environment.yaml
conda activate mixed-ehrhart
python -m pip install --no-build-isolation -e .
sage -n jupyterlab
```

Open a notebook in `notebooks/` and select its SageMath kernel. Each runs from the first cell.

## Background

For lattice polytopes $P_1,\ldots,P_k\subseteq\mathbb R^d$, the library computes

$$
E(n)=\\bigl((n_1P_1+\cdots+n_kP_k)\cap\mathbb Z^d\bigr)
=\sum_{|\alpha|\le D}c_\alpha\prod_i\binom{n_i}{\alpha_i}, \quad c_\alpha=\Delta^\alpha E(0),
$$

where $D=\dim(P_1+\cdots+P_k)$. `multivariate_ehrhart` returns both the ordinary polynomial over `QQ` and the dictionary of integer forward differences $c_\alpha$. These are coefficients of the multivariate Ehrhart polynomial in the binomial basis. 

The coefficients $c_\alpha=\Delta^\alpha E(0)$ are the combinatorial mixed volumes, given by the inclusion–exclusion sum

$$
CMV(P_1,\ldots,P_k)
=\sum_{I\subseteq\{1,\ldots,k\}}(-1)^{k-|I|}
\bigl((\sum_{i\in I}P_i)\cap\mathbb Z^d\bigr).
$$

See https://arxiv.org/pdf/1605.07431 or https://hal.science/hal-05774751v1 for more background.

The CMVs are computed using inclusion-exclusion. Equal polytopes are grouped by their vertex sets to avoid repeating unnecessary computations of lattice point counts, and these counts are stored in a dictionary. Dilations and partial Minkowski sums are cached to avoid repeated computations. 
For the full mixed Ehrhart polynomial expansion, we compute the coefficients $c_\alpha$ for all $| \alpha | \leq D$. All intermediate computations are reused for computing these CMVs. To get the polynomials in the monomial basis, the binomial basis polynomials are expanded using exact rational arithmetic.

Notebooks contain small examples of usage.
