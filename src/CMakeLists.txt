cmake_minimum_required(VERSION 3.10)

option(DIMENSION "specify the number of dimensions")

enable_language(CXX) 

# add the executable
add_executable(opar 
  common.cpp 
  currents.cpp
  fields.cpp
  functions.cpp
  globals.cpp
  opar.cpp
  particle_boundary.cpp
  particle_diagnostic.cpp
  particle_exchange.cpp
  particles.cpp
  random.cpp
  species.cpp
  util.cpp
)

set(HDF5_PREFER_PARALLEL TRUE)
set(MPI_CXX_COMPILER TRUE)
set(MPI_CXX_SKIP_MPICXX TRUE)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS})
target_link_libraries(opar "-lschnek" ${MPI_C_LIBRARIES} ${HDF5_LIBRARIES})

set(HAVE_MPI ${MPI_FOUND})
set(HAVE_HDF5 ${HDF5_FOUND})

if (DIMENSION STREQUAL "1")
  set(ONE_DIMENSIONAL 1)
elseif (DIMENSION STREQUAL "3")
  set(THREE_DIMENSIONAL 1)
else()
  set(TWO_DIMENSIONAL 1)
endif()


configure_file(config.hpp.in ../src/config.hpp @ONLY)
