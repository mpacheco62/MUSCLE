title: FEA Integration

# Integration with FEA Solvers

UMatLib is designed to be easily wrapped by standard user subroutines.

## 1. LS-DYNA (UMAT)
To use UMatLib in LS-DYNA, you should map the `sigma` and `eps` arrays to the `ten_3D2Osym` type.

**Example Wrapper Structure:**
1. Convert LS-DYNA history variables to internal library state.
2. Call `closest_point_solve`.
3. Update the `stress` and `hsv` arrays.

## 2. ANSYS APDL (USERMAT)
For ANSYS, the consistent tangent matrix (Jacobian) must be provided in the `dsde` array. You can obtain it directly from the solver:

```fortran
call solver%tangent(strain, data, tangent_tensor)
! Export tangent_tensor%vals to ANSYS dsde format
```
