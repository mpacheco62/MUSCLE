title: Math Modules

# The math modules (`src/math/`)

Small dense linear algebra and spectral routines for the point kernels, in `sp` and `dp`. Every routine
marked *kernel* is `pure`, has `@:DEVICE()` and can run inside an OpenMP target region (no LAPACK, no
allocation, no IEEE modules).

## Two layers

The tensor types use part of the math modules (their inverses), and the rest of the math modules use the
tensor types. A routine therefore belongs to one of two layers, and a module of the lower layer never uses a
tensor:

```
 models (yield criteria, kinematics, solvers)
        │
 math on tensors     muscle_math_eigen, muscle_math_spectral_derivs, muscle_math_matrix_functions
        │
 tensor types        muscle_tensor_* (muscle_tensors)
        │
 math below tensors  muscle_math_inverses, muscle_math_lu_solve, muscle_math_fp

 host only, beside:  muscle_math_linsolvers
```

## What is where

| I need | Routine | Module |
|---|---|---|
| Inverse of a 3x3 … 9x9 matrix, with its condition number | `M33INV`, `M44INV`, `M55INV`, `M66INV`, `M99INV`, `COND_INF` | `muscle_math_inverses` |
| The infinity norm (maximum absolute row sum) of a matrix, in device code | `NORM_INF` | `muscle_math_inverses` |
| A small linear system `A x = b` in a kernel (n = 2 … 9) | `LU<n>SOLVE` (`lu_solve`: host front end) | `muscle_math_lu_solve` |
| A quiet NaN, `exponent(x)` or `scale(x, n)` in device code | `QUIET_NAN_<k>`, `fp_exponent`, `fp_scale` | `muscle_math_fp` |
| Eigenvalues of a symmetric tensor (values only, cheap) | `eigenvals` | `muscle_math_eigen` |
| Eigenvalues and eigenvectors (for spectral gradients/Hessians) | `spectral_decomposition` | `muscle_math_eigen` |
| Eigenvalues and eigenprojections accurate for close eigenvalues | `eigenprojections`, `eigenvals_refined` | `muscle_math_eigen` |
| Derivatives of the eigenvalues (Cayley-Hamilton, no eigenvectors) | `dEigenvalues_dTensor`, `d2Eigenvalues_dTensor2` | `muscle_math_spectral_derivs` |
| Gradient and Hessian of a symmetric function of the eigenvalues | `spectral_gradient`, `spectral_hessian` | `muscle_math_spectral_derivs` |
| Derivative of the square of a symmetric tensor | `dsquare` | `muscle_math_spectral_derivs` |
| Rotation of the polar decomposition F = R U | `polar_rotation` | `muscle_math_matrix_functions` |
| (ln a − ln b)/(a − b), accurate for close a, b | `log_divided_difference` | `muscle_math_matrix_functions` |
| A dense system of any size on the host (Gauss, Gauss-Seidel, LAPACK `?gesv`) | `gauss_solve`, `gauss_seidel`, `mklsolve` | `muscle_math_linsolvers` |

Tensor conversions (`to_matrix`, `from_upper`), norms, determinants and inverses of tensors are not here: they
are free functions of the tensor types (`muscle_tensor_functions`, [Tensor Algebra Guide](tensors.html)).
Numerical derivatives are the macros of `tests/common/muscle_fd.fypp`, for tests only.

## One spectral decomposition

`spectral_decomposition` is the spectral decomposition of the library (D-032, "K"): the deviator scaled by a
power of 2, its closed-form eigenvalues, the eigenvector of the most separated eigenvalue by a cross product of
rows, the remaining pair by one plane rotation without overflow and the eigenvalues as Rayleigh quotients. It is
accurate to a few eps for every eigenvalue gap, at any scale, and is what the anisotropic yield criteria use for
their equivalent stress and its derivatives, and what the finite strain kinematics use through
`eigenprojections` (\( P_a = v_a \otimes v_a \)). `eigenvals` (closed form, values only) and the
Cayley-Hamilton derivatives of `muscle_math_spectral_derivs` remain public routines.

## Shared code

`muscle_math_lu_bodies.fypp` is an include, not a module: the FYPP macros with the LU factorization and the
inverse and solve bodies shared by `muscle_math_inverses` and `muscle_math_lu_solve`, instantiated with a
compile-time size so the loops unroll.
