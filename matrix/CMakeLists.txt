set(HAEDER_LIST "${SPARSE_SOURCE_DIR}/include/matrix/matrix.hpp")
set(FILES_LIST "matrix.cpp" "COO.cpp" "CSR.cpp" "ELL.cpp" "SGP.cpp" "function.cpp" "function_spmvs.cpp")

add_library(matrix_library ${FILES_LIST} ${HEADER_LIST})

target_include_directories(matrix_library PRIVATE ../include)

target_compile_features(matrix_library PUBLIC cxx_std_11)


