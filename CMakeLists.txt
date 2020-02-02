#CMakeLists.txt file

cmake_minimum_required(VERSION 3.9...3.14)
project(SPARSE VERSION 0.1 LANGUAGES CXX)

find_package(MPI)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    message([NOTICE] "MPI was found")
endif(MPI_FOUND)


# Set MPICXX as the CXX compiler

SET(CMAKE_CXX_COMPILER mpicxx)

add_subdirectory(matrix)
add_subdirectory(apps)
