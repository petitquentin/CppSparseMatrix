#CMakeLists.txt file

cmake_minimum_required(VERSION 3.9...3.14)
#project(CppSparseMatrix)

set(SRCS 
    matrix.cpp)

add_library(matrix OBJECT ${SRCS})
#add_subdirectory(headermatrix)
target_include_directories(matrix PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(CppSparseMatrix ${SRCS})

target_link_libraries(CppSparseMatrix libheadermatrix)