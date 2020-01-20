#CMakeLists.txt file

cmake_minimum_required(VERSION 3.9...3.14)
project(SPARSE VERSION 0.1 LANGUAGES CXX)

find_package(MPI)
if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)

add_subdirectory(matrix)
add_subdirectory(apps)
