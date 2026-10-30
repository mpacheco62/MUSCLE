title: Tensor Algebra Guide

# Using the Tensor Engine

MUSCLE provides a high-level syntax for tensor operations to make the code look similar to
mathematical notation. Everything is reached with a single `use muscle_tensors`.

The tensor modules are generated with [FYPP](https://fypp.readthedocs.io) from one source per
family: the same template produces the 3D and 2D types in single (`sp`) and double (`dp`)
precision, so both dimensions and both kinds share exactly the same code and operators.

## Available Types

Type names follow `ten_<dim><family>_<kind>`: `ten_3D2Osym_sp`, `ten_3D2Osym_dp`,
`ten_2D4O3sym_sp`, ... The historical names without suffix (`ten_3D2Osym`, `ten_2D4O3sym`, ...)
are the double precision types, so existing code compiles unchanged.

| Family   | 3D (components)              | 2D (components)              | Description                          |
|----------|------------------------------|------------------------------|--------------------------------------|
| `2Osym`  | `ten_3D2Osym` (6)            | `ten_2D2Osym` (4)            | Symmetric 2nd order                  |
| `2O`     | `ten_3D2O` (9)               | `ten_2D2O` (5)               | General 2nd order                    |
| `4O2sym` | `ten_3D4O2sym` (6x6)         | `ten_2D4O2sym` (4x4)         | 4th order with minor symmetries      |
| `4O3sym` | `ten_3D4O3sym` (21)          | `ten_2D4O3sym` (10)          | Fully symmetric 4th order            |
| `4O`     | `ten_3D4O` (9x9)             | `ten_2D4O` (5x5)             | General 4th order                    |

The components are stored in `vals` in tensorial Voigt notation (the shear is not doubled; the
factor 2 goes into the contractions): 3D symmetric `(11, 22, 33, 12, 23, 13)`, 2D symmetric
`(11, 22, 33, 12)`, 2D general `(11, 21, 12, 22, 33)`.

2D tensors (plane strain / plane stress / axisymmetric) are 3D tensors whose out-of-plane
coupling components (13, 23, 31, 32) vanish; the out-of-plane normal component (33) is kept.
They expose the same operators as their 3D counterparts.

### Identities

The identities do not depend on the dimension: they combine with 3D and 2D tensors.

- Unscaled: `iden_2O` (\(\delta_{ij}\)), `iden_4O3T` (\(\delta_{ij}\delta_{kl}\)) and `iden_4O4T`
  (symmetric identity \(\mathbb{I}^S\)). They store nothing, have no kind and serve both kinds.
- Scaled (\(c\,\mathbf{I}\), ...): `iden_2OS_sp/_dp`, `iden_4O3TS_sp/_dp`, `iden_4O4TS_sp/_dp`, one
  type per kind; `iden_2OS`, `iden_4O3TS`, `iden_4O4TS` are the double precision ones.
- A scalar times an unscaled identity gives the scaled identity of the scalar's kind
  (`2.0_sp * iden_2O()` is an `iden_2OS_sp`). Operations between unscaled identities only
  (`iden_2O() + iden_2O()`, `iden_4O3T() + iden_4O4T()`, ...) have no scalar to fix the kind and
  return double precision.

### Kinds

There are no mixed-precision operators: both operands of an operator must have the same kind
(a scalar too: `2.0_sp * a` for a single precision `a`), and `a_dp + b_sp` does not compile.
Convert explicitly through the components:

```fortran
use muscle_kinds, only: sp, dp
type(ten_3D2Osym)    :: eps        ! = ten_3D2Osym_dp
type(ten_3D2Osym_sp) :: eps_sp
eps%vals = real(eps_sp%vals, dp)
```

## Common Operations

The operators are generic over dimension and kind, also between different types (a 4th order
tensor and a 2nd order one, a tensor and an identity, ...): `+`, `-`, `*` (by a scalar, or single
contraction of 2nd order tensors), `/` (by a scalar), `=` (between compatible types),
`.ddot.`, `.tdot.`, `.tdotsym.`, `.transform.`, `.dev.`, `.inv.` and `.approx.`. The result
type is the natural one (e.g. `ten_3D4O3sym .ddot. ten_3D2Osym` is a `ten_3D2Osym`).

### Double Contraction (`.ddot.`)
Used for \(\pmb{\sigma} = \mathbb{C} : \pmb{\varepsilon}\) or \(\text{tr}(\pmb{\sigma}) = \pmb{\sigma} : \pmb{I}\).

```fortran
use muscle_kinds, only: dp
type(ten_3D2Osym)  :: stress, strain
type(ten_3D4O3sym) :: C
real(dp) :: lambda = 100.0_dp, mu = 80.0_dp

C = lambda * iden_4O3T() + 2.0_dp * mu * iden_4O4T()     ! isotropic stiffness
call strain%init(1.0e-3_dp, 0.0_dp, 0.0_dp, 2.0e-4_dp, 0.0_dp, 0.0_dp)
stress = C .ddot. strain
stress = lambda * (strain .ddot. iden_2O()) * iden_2O() + 2.0_dp * mu * strain   ! same result
```

The sums of two 4th-order identities exist only in 3D (their result is a `ten_3D4O3sym`); in 2D
add the identities to a 2D tensor:

```fortran
type(ten_2D4O3sym) :: C2
C2 = 2.0_dp * mu * iden_4O4T()
C2 = C2 + lambda * iden_4O3T()
```

### Symmetrized Dyadic Product (`.tdotsym.`)
Computes the symmetric part of the outer product \(\mathbf{a} \otimes \mathbf{b}\)
(`.tdot.` is the plain product \(a_{ij} b_{kl}\)).

```fortran
type(ten_3D2Osym)  :: n, m
type(ten_3D4O3sym) :: Rank4
Rank4 = n .tdotsym. m
Rank4 = .tdotsym. n     ! n (x) n
```

### Transformation (`.transform.`)
Congruence transformation \(\mathbf{F}\cdot\mathbf{S}\cdot\mathbf{F}^T\) (and its 4th-order
counterpart), e.g. a push-forward in 2D:

```fortran
type(ten_2D2O)    :: F
type(ten_2D2Osym) :: S, sigma
call F%init(xx=1.1_dp, xy=0.2_dp, yx=0.0_dp, yy=0.95_dp, zz=1.0_dp)
call S%init([100.0_dp, 50.0_dp, 0.0_dp, 10.0_dp])
sigma = (F .transform. S) / F%det()   ! Cauchy stress from the PK2 stress
```

### Inverse
Every tensor type can be inverted in three equivalent ways. For 4th order tensors the inverse is
taken with respect to the double contraction (\( \mathbb{A}^{-1} : \mathbb{A} = \mathbb{I}^S \)
for the symmetric types).

```fortran
type(ten_3D4O3sym) :: Stiff, Compl
real(dp) :: kappa
logical :: ok

Compl = .inv. Stiff                     ! operator
Compl = Stiff%inv()                     ! method
call inverse(Stiff, Compl, ok, kappa)   ! success flag and condition number (optional)
```

If the tensor is singular, `.inv.` and `%inv()` return quiet NaN in every component, and
`inverse` returns `ok = .false.` (and `kappa = huge(1.0_dp)`). `kappa` is the condition number
\( \kappa_\infty \) of the matrix that is inverted. When compiling with `-ffast-math`, check `ok`
rather than testing for NaN.

## Free functions and type-bound procedures

Two equivalent APIs exist for the named operations:

- **Free functions** (`pure`, `type` arguments, usable in GPU kernels): `norm`, `det`, `square`,
  `transpose`, `is_approx`, `inverse`, `convert_3sym`, `convert_2sym`. They are not exported by
  `muscle_tensors` (to avoid clashing with user variables named `norm`, `det`, ...): import them
  from `muscle_tensor_functions`, renaming if needed.
- **Type-bound procedures** (historical API, host code only): `%init()`, the accessors `%xx()`,
  `%xy()`, ..., `%norm()`, `%det()`, `%inv()`, `%is_approx()`, ... They are thin wrappers over the
  free functions; being type-bound they take a polymorphic `class` argument, so kernels must not
  call them.

```fortran
use muscle_tensors
use muscle_tensor_functions, only: norm, det, inverse
type(ten_3D2Osym) :: s
print *, norm(s), det(s)      ! free functions (kernels)
print *, s%norm(), s%det()    ! type-bound (host only)
```

The operators are free `pure` procedures too and can be used in kernels.

## Formatted output

`write(unit, *) t` and `write(unit, '(DT"LIST")') t` print a tensor only when the library is
built with `ENABLE_UDTIO` (CMake option `MUSCLE_UDTIO`, `AUTO` by default: enabled except with
gfortran < 16.3, whose UDTIO is not thread-safe under OpenMP). Otherwise print `t%vals`.
Formatted output is host-only.
