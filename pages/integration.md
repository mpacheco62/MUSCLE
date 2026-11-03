title: FEA Integration

# Integration with FEA Solvers

MUSCLE ships no user subroutine for any FE code (D-018.17). Each declared material
`muscle_mat_<name>` (generated from the material registry, see the module documentation and the
[Developer Guide](contributing.html)) provides
the pieces an adapter needs, in single (`_sp`) and double (`_dp`) precision:

| Piece | Purpose |
|-------|---------|
| `mat_<name>_<k>` | Material type (parameters of its models, `iter_nw`, `tol`). |
| `mat_init_<name>(mat, props, status)`, `MAT_<NAME>_NPROPS` | Material from a flat array of parameters (order: see *Material parameters*). |
| `MAT_<NAME>_NSTATE` | Number of state variables of a point (the stress is not part of them). |
| `muscle_point_<name>(...)` | Point API on flat arrays (no derived tensor types). |
| `mat_step(mat, dt, coefs, status)` | Only materials with a viscoelastic (frozen) elasticity: coefficients of a time step. |
| `mat_rotate_hsv_<name>(hsv, drot)` | Only materials with a strain input: rotates the tensorial parts of the state. |
| `muscle_update_<name>(...)` | Batch API over SoA arrays (OpenMP / GPU), not used by the examples below. |

The point API takes the **total** strain at \( t_{n+1} \) (or the deformation gradients at \( t_n \)
and \( t_{n+1} \) for a finite strain material) and the converged state, and returns the stress,
the new state and, optionally, the tangent:

```fortran
! small strain or corotational material
call muscle_point_<name>(mat, strain, dt, hsv_n, hsv_np1, stress, status, iters [, ddsdde])
! frozen viscoelastic elasticity: coefficients of the step first (once per step)
call mat_step(mat, dt, coefs, status)
call muscle_point_<name>(mat, strain, dt, coefs, hsv_n, hsv_np1, stress, status, iters [, ddsdde])
! finite strain material (total/updated Lagrangian, logarithmic, small strain kinematics on F)
call muscle_point_<name>(mat, F_n, F_np1, dt, hsv_n, hsv_np1, stress, status, iters [, ddsdde])
```

`hsv_n` and `hsv_np1` must be different arrays (copy the FE history array to a local one first).
`status /= STATUS_OK` (module `muscle_status`) means the outputs must not be used: ask the FE code
for a smaller increment (`PNEWDT` in Abaqus, `keycut` in ANSYS).

## Material parameters

`mat_init_<name>(mat, props, status)` builds the material from the first `MAT_<NAME>_NPROPS` values of
`props`. More values are ignored. `STATUS_INVALID_PARAMETER` is returned for fewer values, for a value
outside the domain that `mat_init_<name>` checks (`tol` or `iter_nw` <= 0, Young modulus <= 0, Poisson
ratio outside (-1, 1/2), Drucker-Prager hardening mode outside 1..3), and by the initializers of the
models that validate their inputs (e.g. a Maxwell relaxation time \( \tau \le 0 \)). The order is
generated per material from the model metadata (D-018.17e):

- the inputs of the model of each component, in the order of the components of the material type
  (elasticity, or the viscoelastic model of a frozen viscoelastic elasticity; yield criterion;
  hardening, overstress or viscoplastic law and the hardening law it composes; a kinematics has no
  parameters). A model with an initializer takes its inputs (linear elasticity: Young modulus and
  Poisson ratio; Drucker-Prager: friction angle, K, hardening mode; CPB06/Yld2004: the coefficients),
  the others their parameter components (including those with a default value, e.g. `epmin`);
- Maxwell: bulk, shear_inf, then `visco_branches` pairs (shear, tau);
- an integer input (Drucker-Prager hardening mode: 1 tension, 2 compression, 3 shear) is given as a
  real and rounded with `nint`;
- last, the solver settings `tol` (relative tolerance of the Newton iterations; the default of the
  library is `NEWTON_TOL_dp` = 1e-10, `NEWTON_TOL_sp` = 1e-5) and `iter_nw` (maximum number of
  iterations, e.g. 200), both > 0.

The table of a material is printed in Markdown by

```bash
python3 tools/fypp/muscle_meta.py props <name> [registry.py]
```

e.g. for `cp_vm_swift`: `young`, `poisson`, `k`, `e0`, `n`, `tol`, `iter_nw`
(`MAT_CP_VM_SWIFT_NPROPS` = 7).

## Conventions and component order

MUSCLE stores symmetric tensors in **tensorial** Voigt order `(11, 22, 33, 12, 23, 13)`, *without*
the factor 2 of the shear strains. The tangent `ddsdde(i, j)` is
\( \partial\sigma_i/\partial\varepsilon_j \) for \( j = 1..3 \) and
\( \partial\sigma_i/\partial\gamma_j \) with the engineering shear \( \gamma_j = 2\varepsilon_j \) for
\( j = 4..6 \), the convention of the FE codes below. Deformation gradients are 9 reals in
column-major order (`reshape(F, [9])` of a `F(3, 3)` array).

| Code | Stress/strain order | Shear strain | MUSCLE component of code component `i` |
|------|---------------------|--------------|-----------------------------------------|
| MUSCLE | 11, 22, 33, 12, 23, 13 | tensorial \( \varepsilon_{12} \) | — |
| Abaqus (`STRESS`, `STRAN`, `DDSDDE`) | 11, 22, 33, 12, **13, 23** | engineering \( \gamma_{12} \) | `[1, 2, 3, 4, 6, 5]` |
| ANSYS (`stress`, `Strain`, `dsdePl`) | 11, 22, 33, 12, 23, 13 | engineering \( \gamma_{12} \) | `[1, 2, 3, 4, 5, 6]` |
| LS-DYNA (`sig`, `eps`) | xx, yy, zz, xy, yz, zx | engineering \( \gamma_{xy} \) (increments) | `[1, 2, 3, 4, 5, 6]` |

So an adapter permutes the components (only Abaqus: swap 5 and 6, in rows *and* columns of the
tangent) and divides the shear strains of the FE code by 2. Stresses and tangents need no scaling.

Precision: Abaqus/Standard and ANSYS call the user routine in double precision (`_dp`); the single
precision build of LS-DYNA can use the `_sp` material.

Abaqus/Explicit does not call a UMAT but a **VUMAT**: the table and the UMAT skeleton below do not apply
to it. A VUMAT receives a block of `nblock` points (arrays `(nblock, ...)`, the layout of the batch API
`muscle_update_<name>`), strain **increments** with *tensorial* shears in the order 11, 22, 33, 12,
23, 31 (= MUSCLE order, no factor 2), stresses in the corotational frame, and `stateOld`/`stateNew`
that it does not rotate (the rotation increment follows from the polar decompositions of
`defgradOld`/`defgradNew`). An adapter keeps the total strain in the state, rotates the tensors of the
state (`mat_rotate_hsv_<name>` per point) and calls the batch API on the block; `_sp` for the single
precision build. Not written nor tested.

## Rotations and state variables

For a small strain or corotational material the FE code passes the strain (and stress) already
rotated to the corotational frame (`kin='corotational'` in the registry, D-018.15f). The plastic
strain and the tensors of a viscoelastic state stored in the history array must be rotated too:

- **Abaqus** does not rotate `STATEV`: call `mat_rotate_hsv_<name>(hsv, DROT)` before the point API.
- **ANSYS**: *open question*, it is not confirmed whether `ustatev` is rotated by ANSYS with
  large deflections (`NLGEOM,ON`). Check it with a rigid rotation test before relying on either
  choice. USERMAT receives no rotation increment: if ANSYS does not rotate `ustatev`, compute it from
  the polar decompositions \( \mathbf{F} = \mathbf{R}\mathbf{U} \) of `defGrad_t` and `defGrad`,
  \( \Delta\mathbf{R} = \mathbf{R}_{n+1}\mathbf{R}_n^T \).
- **LS-DYNA** rotates `sig` (objective stress update) but not `hsv`; `umat41` receives no rotation
  increment. The total strain an adapter accumulates in `hsv` (below) is also a tensor in the unrotated
  frame, like the plastic strain. With large rotations use a finite strain material instead
  (deformation gradient, `IHYPER=1` on `*MAT_USER_DEFINED_MATERIAL_MODELS`): `umat41` then receives
  only the current F (in `hsv`), so the adapter must also store F_n in `hsv` (9 more history
  variables) to call the point API with `F_n` and `F_np1`.

A finite strain material (input F) keeps its state objective by itself: do not rotate its history.

## Abaqus (UMAT)

Skeleton, not compiled nor tested against Abaqus. `<name>` is a material of the registry with a
strain input; `*DEPVAR` must be at least `MAT_<NAME>_NSTATE`, `*USER MATERIAL, CONSTANTS=` equal to
`MAT_<NAME>_NPROPS`.

```fortran
subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, &
                time, dtime, temp, dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, &
                nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, jstep, kinc)
    use muscle_kinds, only : dp
    use muscle_status, only : STATUS_OK
    use muscle_mat_<name>, only : mat_<name>_dp, mat_init_<name>, muscle_point_<name>, &
        mat_rotate_hsv_<name>, MAT_<NAME>_NSTATE
    implicit none
    ! ... declarations of the UMAT arguments (double precision, ntens = 6) ...
    integer, parameter :: PERM(6) = [1, 2, 3, 4, 6, 5]   ! Abaqus <-> MUSCLE (its own inverse)
    real(dp), parameter :: SHEAR(6) = [1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp, 0.5_dp]
    type(mat_<name>_dp) :: mat
    real(dp) :: strain(6), sig(6), tangent(6, 6), hsv_n(MAT_<NAME>_NSTATE), hsv_np1(MAT_<NAME>_NSTATE)
    integer :: status, iters

    ! The material could be built once per analysis and cached (thread safety is up to the adapter)
    call mat_init_<name>(mat, props(1:nprops), status)
    if (status /= STATUS_OK) call xit()
    ! Total strain at t_n+1 in the corotational frame, tensorial shears, MUSCLE order
    strain = (stran(PERM) + dstran(PERM))*SHEAR
    hsv_n = statev(1:MAT_<NAME>_NSTATE)
    call mat_rotate_hsv_<name>(hsv_n, drot)       ! Abaqus does not rotate STATEV
    call muscle_point_<name>(mat, strain, dtime, hsv_n, hsv_np1, sig, status, iters, tangent)
    if (status /= STATUS_OK) then
        pnewdt = 0.5_dp                            ! ask for a smaller increment
        return
    end if
    stress(1:6) = sig(PERM)
    ddsdde(1:6, 1:6) = tangent(PERM, PERM)
    statev(1:MAT_<NAME>_NSTATE) = hsv_np1
end subroutine umat
```

With a finite strain material call
`muscle_point_<name>(mat, reshape(dfgrd0, [9]), reshape(dfgrd1, [9]), dtime, ...)` instead and do
not rotate `statev`. With a viscoelastic elasticity, `call mat_step(mat, dtime, coefs, status)`
before the point API and pass `coefs` after `dtime`.

## ANSYS (USERMAT, 3D)

Skeleton, not compiled nor tested against ANSYS. Same component order as MUSCLE; engineering shear
strains.

```fortran
subroutine usermat3d(matId, elemId, kDomIntPt, kLayer, kSectPt, ldstep, isubst, keycut, nDirect, &
                     nShear, ncomp, nStatev, nProp, Time, dTime, Temp, dTemp, stress, ustatev, dsdePl, &
                     sedEl, sedPl, epseq, Strain, dStrain, epsPl, prop, coords, var0, defGrad_t, defGrad, &
                     tsstif, epsZZ, cutFactor, var1, var2, var3, var4, var5, var6, var7)
    use muscle_kinds, only : dp
    use muscle_status, only : STATUS_OK
    use muscle_mat_<name>, only : mat_<name>_dp, mat_init_<name>, muscle_point_<name>, MAT_<NAME>_NSTATE
    implicit none
    ! ... declarations of the USERMAT arguments (ncomp = 6) ...
    real(dp), parameter :: SHEAR(6) = [1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp, 0.5_dp]
    type(mat_<name>_dp) :: mat
    real(dp) :: strain6(6), hsv_n(MAT_<NAME>_NSTATE), hsv_np1(MAT_<NAME>_NSTATE)
    integer :: status, iters

    call mat_init_<name>(mat, prop(1:nProp), status)
    if (status /= STATUS_OK) then
        keycut = 1
        return
    end if
    strain6 = (Strain(1:6) + dStrain(1:6))*SHEAR
    hsv_n = ustatev(1:MAT_<NAME>_NSTATE)
    ! OPEN QUESTION: whether ANSYS rotates ustatev with NLGEOM,ON is not confirmed. If it does not,
    ! rotate hsv_n here with mat_rotate_hsv_<name>(hsv_n, drot): USERMAT gives no rotation increment,
    ! drot = R_np1 R_n^T from the polar decompositions of defGrad_t and defGrad.
    call muscle_point_<name>(mat, strain6, dTime, hsv_n, hsv_np1, stress(1:6), status, iters, dsdePl)
    if (status /= STATUS_OK) then
        keycut = 1                                 ! ask for a bisection
        return
    end if
    ustatev(1:MAT_<NAME>_NSTATE) = hsv_np1
    epseq = hsv_np1(7)                             ! equivalent plastic strain
end subroutine usermat3d
```

## LS-DYNA (umat41)

Skeleton, not compiled nor tested against LS-DYNA. `umat41` receives strain **increments** and no
tangent output (explicit; the implicit tangent goes in `utan41`): the total strain is kept in 6
extra history variables after the state of the material (`NHV` on the material card at least
`MAT_<NAME>_NSTATE + 6`). The material constants `cm(*)` are the `props` of the material; LS-DYNA
also needs the bulk and shear moduli at the positions given on the card (`IBULK`/`IG`): put them after
the `MAT_<NAME>_NPROPS` values (they are ignored by `mat_init_<name>`).

```fortran
subroutine umat41(cm, eps, sig, epsp, hsv, dt1, capa, etype, tt, temper, failel, crv, nnpcrv, cma, &
                  qmat, elsiz, idele, reject)
    use muscle_kinds, only : dp
    use muscle_status, only : STATUS_OK
    use muscle_mat_<name>, only : mat_<name>_dp, mat_init_<name>, muscle_point_<name>, &
        MAT_<NAME>_NSTATE, MAT_<NAME>_NPROPS
    implicit none
    ! ... declarations of the umat41 arguments ...
    real(dp), parameter :: SHEAR(6) = [1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp, 0.5_dp, 0.5_dp]
    integer, parameter :: NS = MAT_<NAME>_NSTATE
    type(mat_<name>_dp) :: mat
    real(dp) :: strain(6), hsv_np1(NS)
    integer :: status, iters

    call mat_init_<name>(mat, cm(1:MAT_<NAME>_NPROPS), status)
    if (status /= STATUS_OK) error stop 'umat41: invalid parameters of mat_init_<name>'
    ! Total strain: accumulated increments (xx, yy, zz, xy, yz, zx = MUSCLE order), tensorial shears;
    ! not rotated by LS-DYNA (see Rotations and state variables)
    hsv(NS + 1:NS + 6) = hsv(NS + 1:NS + 6) + eps(1:6)*SHEAR
    strain = hsv(NS + 1:NS + 6)
    call muscle_point_<name>(mat, strain, dt1, hsv(1:NS), hsv_np1, sig(1:6), status, iters)
    if (status /= STATUS_OK) failel = .true.       ! or another strategy of the analysis
    hsv(1:NS) = hsv_np1
end subroutine umat41
```

LS-DYNA rotates `sig` but not `hsv`; see *Rotations and state variables* above for large rotations.

Code written for the former object-oriented API (solver objects, `Plastic_material_history`): see the
[migration guide](migration.html).
