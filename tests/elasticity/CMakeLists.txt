# Kernel API of the elastic family (D-017.12) in every kind and dimension: the linear model through
# the generics of muscle_elas_laws against Hooke law (dp reference), 2D vs 3D embedding of
# every function, the generics against the dp specific functions (bit for bit), and pure @:DEVICE()
# point procedures
add_muscle_fypp_test(elas_kinds test_muscle_elas_laws_kinds.fypp)

# Compile-fail check (D-017.10b): every component of the kernel type elas_linear_<k> is derived
# (private), so it can only be built with elas_linear_init. Passes when the compiler rejects
# linear_private_params.f90 with its "PRIVATE component" error; any other error (a missing module, a
# fatal error) fails it. Only with gfortran and with ifx on Windows (the message of other compilers is
# not known here; the file name itself contains "private", so a loose pattern would always pass: ifx is
# matched by its error number #8289, measured with ifx 2026.2 on Windows).
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(_muscle_syntax_only -fsyntax-only)
    set(_muscle_include_flag -I)
    set(_muscle_pass_regex "PRIVATE component")
    set(_muscle_fail_regex "Cannot open module|Fatal Error")
elseif(WIN32 AND CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    set(_muscle_syntax_only /nologo /syntax-only)
    set(_muscle_include_flag /I)
    set(_muscle_pass_regex "error #8289")
    set(_muscle_fail_regex "error #7002|catastrophic error")
endif()
if(_muscle_syntax_only)
    add_test(NAME elasticity-linear_private_params
             COMMAND ${CMAKE_Fortran_COMPILER} ${_muscle_syntax_only}
                     ${_muscle_include_flag}$<TARGET_PROPERTY:muscle,Fortran_MODULE_DIRECTORY>
                     ${CMAKE_CURRENT_SOURCE_DIR}/linear_private_params.f90)
    set_tests_properties(elasticity-linear_private_params PROPERTIES
                         PASS_REGULAR_EXPRESSION "${_muscle_pass_regex}"
                         FAIL_REGULAR_EXPRESSION "${_muscle_fail_regex}")
endif()
