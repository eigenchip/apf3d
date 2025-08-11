# Fast Multi-Source Nanophotonic Simulations

A computational method for solving Maxwell's equations in complex nanophotonic systems using augmented partial factorization. Achieves 1,000-30,000,000x speedup over conventional methods for multi-channel electromagnetic simulations.

## Theory

### Maxwell's Equations in Frequency Domain
After discretization, Maxwell's equations become a system of linear equations:

```math
Ax_m = b_m
```

Where A is the sparse Maxwell differential operator, `b_m` specifies the mth input source, and `x_m` contains the full-basis solution.

### Generalized Scattering Matrix
The linear response of any system is described by a generalized scattering matrix S that relates input vector v to output vector u:

```math
u = Sv
```

The M columns of S correspond to M distinct inputs (different angles, beam profiles, or waveguide modes).

### Augmented Partial Factorization 
Instead of computing the full solutions `X = A^(-1)B`, we directly calculate the scattering matrix:

```math
S = CA^{-1}B - D
```

Where:
- C projects solutions onto outputs of interest
- B contains all input source profiles
- D subtracts baseline contribution

This is achieved by building an augmented sparse matrix K and performing partial factorization:

```math
K = \begin{bmatrix}
A & B \\
C & D
\end{bmatrix} = LU = \begin{bmatrix}
L_{11} & 0 \\
L_{21} & L_{22}
\end{bmatrix} \begin{bmatrix}
U_{11} & U_{12} \\
0 & H
\end{bmatrix}
```

The Schur complement H yields the scattering matrix: S = -H.

## Methodology

1. Discretize Maxwell's equations using finite-difference on Yee grid
2. Build augmented sparse matrix K containing Maxwell operator A, source profiles B, and projection profiles C
3. Perform single partial factorization to extract Schur complement
4. Apply compression to matrices B and C when needed (APF-c)
5. Extended the 2D code presented in the paper to 3D metasurfaces and nanophotonic systems

The method bypasses full-basis solutions and eliminates repetition over inputs, computing the entire scattering matrix in one operation.

## Applications
- Disordered media (500λ wide × 100λ)
- Metalenses with high NA (11+ million pixels with 3,761 input channels)


## Tools
- MUMPS package for sparse matrix factorization
- MESTI (Maxwell's Equations Solver with Thousands of Inputs)
- MATLAB


## References

1. Lin, H.-C., Wang, Z. & Hsu, C.W. Fast multi-source nanophotonic simulations using augmented partial factorization. Nature Computational Science 2, 815–822 (2022).
2. Popoff, S. M. et al. Measuring the transmission matrix in optics. Physical Review Letters 104, 100601 (2010).
3. Rotter, S. & Gigan, S. Light fields in complex media: mesoscopic scattering meets wave control. Reviews of Modern Physics 89, 015005 (2017).
4. Fisher, D. S. & Lee, P. A. Relation between conductivity and transmission matrix. Physical Review B 23, 6851–6854 (1981).
5. Zhang, F. The Schur Complement and Its Applications (Springer, 2015).

---

**Code Repository**: https://github.com/complexphoton/MESTI.m  
**DOI**: https://doi.org/10.1038/s43588-022-00370-6
