title: Tensor Algebra Guide

# Using the Tensor Engine

MUSCLE provides a high-level syntax for tensor operations to make the code look similar to mathematical notation.

## Available Types

### 3D
- `ten_3D2Osym`: Symmetric 2nd order tensor (6 components).
- `ten_3D2O`: General 2nd order tensor (9 components).
- `ten_3D4O2sym`: 4th order tensor with minor symmetries (6x6 Voigt matrix).
- `ten_3D4O3sym`: Fully symmetric 4th order tensor (21 components).
- `ten_3D4O`: General 4th order tensor (9x9 matrix).

### 2D (plane strain / plane stress / axisymmetric)
2D tensors are 3D tensors whose out-of-plane coupling components (13, 23, 31, 32)
vanish; the out-of-plane normal component (33) is kept. They expose the same
operators as their 3D counterparts.

- `ten_2D2Osym`: Symmetric 2nd order tensor (4 components: 11, 22, 33, 12).
- `ten_2D2O`: General 2nd order tensor (5 components: 11, 21, 12, 22, 33).
- `ten_2D4O2sym`: 4th order tensor with minor symmetries (4x4 Voigt matrix).
- `ten_2D4O3sym`: Fully symmetric 4th order tensor (10 components).
- `ten_2D4O`: General 4th order tensor (5x5 matrix).

```fortran
type(ten_2D2O)     :: F
type(ten_2D2Osym)  :: S, sigma
type(ten_2D4O3sym) :: C, Cinv
call F%init(xx=1.1D0, xy=0.2D0, yx=0D0, yy=0.95D0, zz=1D0)
sigma = (F .transform. S) / F%det()   ! push-forward of the PK2 stress
Cinv  = .inv. C                       ! compliance tensor
```

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
