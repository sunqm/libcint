add_executable(cartesian_c c_call_cartesian.c        )
add_executable(spheric_c   c_call_spheric.c          )
add_executable(spinor_c    c_call_spinor.c           )
add_executable(cartesian_f fortran_call_cartesian.F90)
add_executable(spheric_f   fortran_call_spheric.F90  )
add_executable(spinor_f    fortran_call_spinor.F90   )
add_executable(time_c2h6   time_c2h6.c               )
add_executable(time_c60    time_c60.c                )

target_link_libraries(cartesian_c cint m)
target_link_libraries(spheric_c   cint m)
target_link_libraries(spinor_c    cint m)
target_link_libraries(cartesian_f cint m)
target_link_libraries(spheric_f   cint m)
target_link_libraries(spinor_f    cint m)
target_link_libraries(time_c2h6   cint m)
target_link_libraries(time_c60    cint m)

set_target_properties(time_c2h6 PROPERTIES
  COMPILE_FLAGS ${OpenMP_C_FLAGS}
  LINK_FLAGS ${OpenMP_C_FLAGS})
set_target_properties(time_c60 PROPERTIES
  COMPILE_FLAGS ${OpenMP_C_FLAGS}
  LINK_FLAGS ${OpenMP_C_FLAGS})

