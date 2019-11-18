#CMakeLists.txt file

cmake_minimum_required(VERSION 3.9...3.14)
project(headermatrix ${VerifyCXX})

set(HEADERS
    Matrix.hpp)

set(SRC 
    Matrix.hpp)

add_library(headermatrix ${HEADERS})