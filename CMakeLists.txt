cmake_minimum_required(VERSION 3.10)

# set the project name
project(opar)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Boost REQUIRED)
find_package(Schnek REQUIRED)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# add the executable
add_executable(opar1d
    src/common.cpp
    src/util.cpp
    src/species.cpp
    src/random.cpp
    src/opar.cpp
    src/globals.cpp
    src/currents.cpp
    src/particle_diagnostic.cpp
    src/particle_boundary.cpp
    src/particles.cpp
    src/functions.cpp
    src/particle_exchange.cpp
    src/fields.cpp
)

add_executable(opar2d
    src/common.cpp
    src/util.cpp
    src/species.cpp
    src/random.cpp
    src/opar.cpp
    src/globals.cpp
    src/currents.cpp
    src/particle_diagnostic.cpp
    src/particle_boundary.cpp
    src/particles.cpp
    src/functions.cpp
    src/particle_exchange.cpp
    src/fields.cpp
)

add_executable(opar3d
    src/common.cpp
    src/util.cpp
    src/species.cpp
    src/random.cpp
    src/opar.cpp
    src/globals.cpp
    src/currents.cpp
    src/particle_diagnostic.cpp
    src/particle_boundary.cpp
    src/particles.cpp
    src/functions.cpp
    src/particle_exchange.cpp
    src/fields.cpp
)

target_compile_definitions(opar1d PRIVATE HUERTO_ONE_DIM ONE_DIMENSIONAL)
target_compile_definitions(opar2d PRIVATE HUERTO_TWO_DIM TWO_DIMENSIONAL)
target_compile_definitions(opar3d PRIVATE HUERTO_THREE_DIM THREE_DIMENSIONAL)

function(setoptions target)
    target_include_directories(${target} PUBLIC ${MPI_INCLUDE_PATH})
    target_include_directories(${target} PUBLIC ${HDF5_INCLUDE_DIRS})

    target_link_libraries(${target} ${MPI_C_LIBRARIES})
    target_link_libraries(${target} ${HDF5_LIBRARIES})
    target_link_libraries(${target} ${Schnek_LIBRARIES})
    target_link_libraries(${target} schnek)
endfunction()

setoptions(opar1d)
setoptions(opar2d)
setoptions(opar3d)

message("MPI ${MPI_INCLUDE_PATH} | ${MPI_C_LIBRARIES}")
message("HDF5 ${HDF5_INCLUDE_DIRS} | ${HDF5_LIBRARIES}")
