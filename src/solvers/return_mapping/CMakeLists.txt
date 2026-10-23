target_sources(muscle PRIVATE
                 muscle_solver_closest_point_2.f90
                 muscle_solver_closest_point_time.f90
                 muscle_solver_cutting_plane.f90
                 )
muscle_add_source(muscle_solver_closest_point.f90)
