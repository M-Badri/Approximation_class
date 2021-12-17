add_executable(example_1_2d simple_2d_linear_interpolation_example.f90)

#add_subdirectory(evolution)
#add_subdirectory(initial)
#add_subdirectory(io)
#add_subdirectory(parser)

target_link_libraries(example_1_2d
  PRIVATE
    app_class
  )
