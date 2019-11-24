set(HAEDER_LIST "${SPARSE_SOURCE_DIR}/include/matrix/matrix.hpp")

add_library(matrix_library matrix.cpp ${HEADER_LIST})

target_include_directories(matrix_library PRIVATE ../include)

target_compile_features(matrix_library PUBLIC cxx_std_11)


