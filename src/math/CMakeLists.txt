muscle_add_fypp_sources(muscle
                 muscle_math_inverses.fypp
                 muscle_math_lu_solve.fypp
                 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/muscle_math_lu_inverse.fypp
                 )

muscle_add_fypp_sources(muscle
                 muscle_math_operations.fypp
                 muscle_math_derivatives.fypp
                 muscle_math_spectral_derivs.fypp
                 muscle_math_linsolvers.fypp
                 )
