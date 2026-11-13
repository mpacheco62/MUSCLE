muscle_add_fypp_sources(muscle
                 muscle_math_fp.fypp
                 muscle_math_inverses.fypp
                 muscle_math_lu_solve.fypp
                 DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/muscle_math_lu_bodies.fypp
                 )

muscle_add_fypp_sources(muscle
                 muscle_math_matrix_functions.fypp
                 muscle_math_spectral_derivs.fypp
                 muscle_math_eigen.fypp
                 muscle_math_linsolvers.fypp
                 )
