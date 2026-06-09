title: Theoretical Background
summary: Mathematical foundation of the library.

# Theoretical Foundation

This page describes the mathematical models implemented in MUSCLE.

## 1. Plasticity and Return Mapping
The library uses an implicit integration scheme. Given a strain increment \(\Delta\varepsilon\), we find the state at \(n+1\) by solving:

\[ \pmb{\sigma}_{n+1} = \mathbb{C} : ( \pmb{\varepsilon}_{n+1} - \pmb{\varepsilon}_{n+1}^p ) \]

The **Closest Point Projection** algorithm minimizes the distance to the yield surface \(f(\sigma, q) = 0\) using a Newton-Raphson iteration.

## 2. Yield Criteria

### Von Mises
The classic J2 flow theory where the equivalent stress is:
\[ f(\pmb{\sigma}) = \sqrt{\frac{3}{2}\pmb{S}:\pmb{S}} \]

### CPB06 (Cazacu-Plunkett-Barlat)
An advanced criterion for orthotropic materials with strength differential effects:
\[ f(\pmb{\sigma}) = B\left[\sum_{i=1}^{3} (|\Sigma_i| - k\Sigma_i)^a \right]^{1/a} \]
Where \(\Sigma_i\) are the eigenvalues of the transformed deviatoric stress tensor \(\pmb{\Sigma} = \mathbb{L} : \pmb{S}\).

## 3. Hardening Laws
We implement several isotropic hardening models:
- **Swift:** \(\sigma_y = K(\varepsilon_0 + \bar{\varepsilon}^p)^n\)
- **Voce:** \(\sigma_y = k \varepsilon_p + q [ 1 - e^{-n \varepsilon_p} ]\)