

add_library(app_class "")

#target_sources(parser
#  PUBLIC
#    ${CMAKE_CURRENT_LIST_DIR}/parser.f90
#  )


target_sources(app_class
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/mod_feqparse.f90
    ${CMAKE_CURRENT_LIST_DIR}/mod_aux_string_functions.f90
    ${CMAKE_CURRENT_LIST_DIR}/mod_math_lib.f90
    ${CMAKE_CURRENT_LIST_DIR}/mod_ac_linear.f90
    ${CMAKE_CURRENT_LIST_DIR}/mod_aux_procs_for_tests.f90
    ${CMAKE_CURRENT_LIST_DIR}/mod_aux_precision.f90
  )




#target_link_libraries(main.bin
#        PRIVATE utils
#        PRIVATE allocators)

#add_subdirectory(allocators)
