title: Tensor Algebra Guide

# Using the Tensor Engine

UMatLib provides a high-level syntax for tensor operations to make the code look similar to mathematical notation.

## Available Types
- `ten_3D2Osym`: Symmetric 2nd order tensor (6 components).
- `ten_3D2O`: General 2nd order tensor (9 components).
- `ten_3D4O3sym`: Fully symmetric 4th order tensor (21 components).

## Common Operations

### Double Contraction (`.ddot.`)
Used for \(\pmb{\sigma} = \mathbb{C} : \pmb{\varepsilon}\) or \(\text{tr}(\pmb{\sigma}) = \pmb{\sigma} : \pmb{I}\).

```fortran
type(ten_3D2Osym) :: stress, strain
type(ten_3D4O3sym) :: C
stress = C .ddot. strain
```

### Symmetrized Dyadic Product (`.tdotsym.`)
Computes the symmetric part of the outer product \(\mathbf{a} \otimes \mathbf{b}\).

```fortran
type(ten_3D2Osym) :: n, m
type(ten_3D4O3sym) :: Rank4
Rank4 = n .tdotsym. m
```

### Inverse of 4th Order Tensors
You can invert the 6x6 matrix representation of fully symmetric tensors:

```fortran
type(ten_3D4O3sym) :: Stiff, Compl
S = .inv. Stiff
```
