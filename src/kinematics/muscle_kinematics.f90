! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kinematics
    !! Module muscle_kinematics
    !! ========================
    !! Umbrella module that re-exports all continuum mechanics kinematic types and strategies.

    use muscle_kinematics_base,          only : Base_kinematics
    use muscle_kin_small_strain,         only : Small_strain_kinematics
    use muscle_kin_finite_base,          only : Base_F_kinematics
    use muscle_kin_total_lagrangian,     only : Total_lagrangian_kinematics
    use muscle_kin_updated_lagrangian,   only : Updated_lagrangian_kinematics
    use muscle_kin_corotational,         only : Corotational_kinematics
    use muscle_kin_material_logarithmic, only : Material_logarithmic_kinematics
    use muscle_kin_spatial_logarithmic,  only : Spatial_logarithmic_kinematics

    implicit none
    public :: Base_kinematics
    public :: Small_strain_kinematics
    public :: Base_F_kinematics
    public :: Total_lagrangian_kinematics
    public :: Updated_lagrangian_kinematics
    public :: Corotational_kinematics
    public :: Material_logarithmic_kinematics
    public :: Spatial_logarithmic_kinematics

end module muscle_kinematics