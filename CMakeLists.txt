cmake_minimum_required(VERSION 3.5)



### options
option(WARNINGS "ON - display compiler warnings, OFF - do not display compiler warnings" ON)
option(OPENMP "ON - enables adding of OpenMP targets, OFF - disables adding of OpenMP targets" OFF)
# option(PFUNIT "ON - using pFUnit for unit tests, OFF - skip building unit tests" OFF)

# defines project name and used language
project(Approximation_class Fortran)








# expose options to the user
#option(USE_LIBRARY "Compile sources into a library" OFF)

#message(STATUS "Compile sources into a library? ${USE_LIBRARY}")

#include(CMakeDependentOption)

# second option depends on the value of the first
#cmake_dependent_option(
#  MAKE_STATIC_LIBRARY "Compile sources into a static library" OFF
#  "USE_LIBRARY" ON
#  )

# third option depends on the value of the first
#cmake_dependent_option(
#  MAKE_SHARED_LIBRARY "Compile sources into a shared library" ON
#  "USE_LIBRARY" ON
#  )

#set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)






# we default to Release build type
#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
#endif()

#message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#message(STATUS "C flags, Debug configuration: ${CMAKE_C_FLAGS_DEBUG}")
#message(STATUS "C flags, Release configuration: ${CMAKE_C_FLAGS_RELEASE}")
#message(STATUS "C flags, Release configuration with Debug info: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
#message(STATUS "C flags, minimal Release configuration: ${CMAKE_C_FLAGS_MINSIZEREL}")

#message(STATUS "C++ flags, Debug configuration: ${CMAKE_CXX_FLAGS_DEBUG}")
#message(STATUS "C++ flags, Release configuration: ${CMAKE_CXX_FLAGS_RELEASE}")
#message(STATUS "C++ flags, Release configuration with Debug info: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
#message(STATUS "C++ flags, minimal Release configuration: ${CMAKE_CXX_FLAGS_MINSIZEREL}")







#message("C++ compiler flags: ${CMAKE_CXX_FLAGS}")

#list(APPEND flags "-fPIC" "-Wall")
#if(NOT WIN32)
#  list(APPEND flags "-Wextra" "-Wpedantic")
#endif()

#add_library(geometry
#  STATIC
#    geometry_circle.cpp
#    geometry_circle.hpp
#    geometry_polygon.cpp
#    geometry_polygon.hpp
#    geometry_rhombus.cpp
#    geometry_rhombus.hpp
#    geometry_square.cpp
#    geometry_square.hpp
#  )

#target_compile_options(geometry
#  PRIVATE
#    ${flags}
#  )

#add_executable(compute-areas compute-areas.cpp)

#target_compile_options(compute-areas
#  PRIVATE
#    "-fPIC"
#  )

#target_link_libraries(compute-areas geometry)




# compactly let the preprocessor know about the compiler vendor
#target_compile_definitions(hello-world
#  PUBLIC "IS_${CMAKE_Fortran_COMPILER_ID}_FORTRAN_COMPILER"
#)



#target_include_directories(hello-embedded-python
#  PRIVATE
#    ${PYTHON_INCLUDE_DIRS}
#  )

#target_link_libraries(hello-embedded-python
#  PRIVATE
#    ${PYTHON_LIBRARIES}
#  )
  
 
 
  
  
############# cmake -D CMAKE_Fortran_COMPILER=gfortran -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ ..
  
#include(GNUInstallDirs)
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# defines targets and sources
#add_subdirectory(src)



#if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall"
#    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
#    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST")
#endif()

#if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
#    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
#    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
#    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native")
#endif()














#if(USE_LIBRARY)
#  message(STATUS "Compile sources into a STATIC library? ${MAKE_STATIC_LIBRARY}")
#  message(STATUS "Compile sources into a SHARED library? ${MAKE_SHARED_LIBRARY}")

#  if(MAKE_SHARED_LIBRARY)
#    add_library(message SHARED ${_sources})

#    add_executable(hello-world hello-world.cpp)

#    target_link_libraries(hello-world message)
#  endif()

#  if(MAKE_STATIC_LIBRARY)
#    add_library(message STATIC ${_sources})

#    add_executable(hello-world hello-world.cpp)

#    target_link_libraries(hello-world message)
#  endif()
#else()
#  add_executable(hello-world hello-world.cpp ${_sources})
#endif()









# introduce a toggle for using a library
#set(USE_LIBRARY OFF)

#message(STATUS "Compile sources into a library? ${USE_LIBRARY}")

# BUILD_SHARED_LIBS is a global flag offered by CMake
# to toggle the behavior of add_library
#set(BUILD_SHARED_LIBS OFF)

# list sources
#list(APPEND _sources Message.hpp Message.cpp)

#if(USE_LIBRARY)
  # add_library will create a static library
  # since BUILD_SHARED_LIBS is OFF
  #add_library(message ${_sources})

  #add_executable(hello-world hello-world.cpp)

  #target_link_libraries(hello-world message)
#else()
#  add_executable(hello-world hello-world.cpp ${_sources})
# endif()














# sets compiler flags
include(cmake/utils/setFlags.cmake)

# sets path for executables
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

# adds custom cmake functions
include(cmake/utils/fortranFunctions.cmake)

# adding pFUnit support
#include(cmake/utils/pFUnit.cmake)

# adds subdirectories
#add_subdirectory(utils)
add_subdirectory(src)
add_subdirectory(examples)

