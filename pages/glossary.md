title: Glossary

# Glossary

Abbreviations used in the names of the declared materials of the test registry
(`tools/fypp/muscle_materials.py`), of the models and of the code.

## Material names

A material of the test registry is named `<solver or kinematics>_<models>`, for example `cp_vm_swift` (closest
point, von Mises, Swift hardening) or `kin_ul_pz_vm_bilinear` (updated Lagrangian kinematics, Perzyna
overstress, von Mises, bilinear hardening).

| Prefix | Meaning |
|--------|---------|
| `cp` | Closest point projection return mapping (`closest_point`, Newton on the plastic multiplier) |
| `cp2` | Monolithic closest point (`closest_point_2`, 7x7 Newton with an own LU solver) |
| `pz` | Closest point with overstress viscoplasticity (`closest_point_perzyna`: Perzyna or Chaboche law) |
| `tm` | Closest point with consistency viscoplasticity (`closest_point_time`: flow stress laws of the `vp` family) |
| `cut` | Cutting plane return mapping (`cutting_plane`, explicit, continuum tangent) |
| `kin` | Finite strain material: the next token is the kinematics, then the solver and the models |

| Token | Model | Family |
|-------|-------|--------|
| `vm` | von Mises | yield |
| `hill48` | Hill (1948) | yield |
| `dp` | extended linear **Drucker-Prager** (not double precision!) | yield |
| `cpb06`, `cpb06ex2` | Cazacu-Plunkett-Barlat 2006 and its two-transformation extension (2008) | yield |
| `yld2004` | Yld2004-18p (Barlat et al. 2005) | yield |
| `bilinear`, `ludwik`, `swift`, `vocemod` | bilinear, Ludwik, Swift, modified Voce | hard |
| `jc`, `jcfull` | Johnson-Cook, full Johnson-Cook | vp |
| `rk`, `mrk` | Rusinek-Klepaczko, modified Rusinek-Klepaczko | vp |
| `nnl` | Nemat-Nasser-Li | vp |
| `va` | Voyiadjis-Almasri | vp |
| `cha` | Chaboche (overstress) | os |
| `maxwell`, `kelvin` | generalized Maxwell, Kelvin-Voigt | visco |
| `ss` | small strain | kin |
| `tl`, `ul` | total, updated Lagrangian | kin |
| `ml`, `sl` | material, spatial logarithmic strain | kin |
| `cr` | corotational (strain rotated by the FE code) | kin |

## Model families (prefix of the generic functions)

| Prefix | Family | Example generic |
|--------|--------|-----------------|
| `elas` | elasticity (linear; frozen viscoelastic) | `elas_stress(prm, strain)` |
| `yield` | yield criteria | `yield_stress_eq(prm, sigma)` |
| `hard` | isotropic hardening laws | `hard_stress(prm, ep)` |
| `vp` | viscoplastic flow stress laws \( \sigma(\epsilon_p, \dot\epsilon_p) \) | `vp_flow_stress(prm, ep, epd, sy)` |
| `os` | overstress viscoplasticity (Perzyna, Chaboche) | `os_flow_rate(prm, sigma_eq, ep, sy)` |
| `visco` | viscoelasticity (generalized Maxwell, Kelvin-Voigt) | `visco_freeze(prm, coefs, state_n)` |
| `kin` | kinematics (small strain, finite strain measures) | `call kin_update(prm, F, state, status)` |

## Other terms

| Term | Meaning |
|------|---------|
| `sp`, `dp` | single and double precision kinds (`muscle_kinds`); a suffix of every type and specific procedure |
| `2O`, `4O`, `sym` | second and fourth order tensors, symmetric storage (`ten_3D2Osym_dp`, `ten_3D4O3sym_dp`) |
| `hsv` | history state variables: the flat array of the state of a point in an FE code (`MAT_<NAME>_NSTATE`) |
| `props` | flat array of material parameters of an FE code (`mat_init_<name>`, `MAT_<NAME>_NPROPS`) |
| SoA | structure of arrays: the batch API takes `eps(npts, 6)`, one row per point |
| kernel | a `pure` procedure that runs per integration point, also on a GPU (`@:DEVICE()`) |
| `ep`, `epd` | equivalent plastic strain and its rate |
| golden | reference values stored in `tests/reference/` from a validated version, compared within the tolerance of their category ([conventions](conventions.html#tests)) |
| FE | finite element (code): Abaqus, ANSYS, LS-DYNA, ... |
