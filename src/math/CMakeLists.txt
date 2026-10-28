target_sources(muscle PRIVATE
                 muscle_math_operations.f90
                 muscle_math_derivatives.f90
                 muscle_math_linsolvers.f90
                 muscle_math_spectral_derivs.f90
                 )

muscle_add_fypp_sources(muscle
                 muscle_math_inverses.fypp
                 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/muscle_math_lu_inverse.fypp
                 )
