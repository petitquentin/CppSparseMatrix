set(HAEDER_LIST "${SPARSE_SOURCE_DIR}/include/matrix/matrix.hpp")
set(FILES_LIST "matrix.cpp" "COO.cpp" "CSR.cpp" "ELL.cpp" "SGP.cpp" "puissancef_time.cpp" "function_spmv.cpp" "puissance_time.cpp" "function_spmvf.cpp" "function_spmv_time.cpp" "function_spmvf_time.cpp" "function_spmvs.cpp" "pageRank.cpp" "pageRankf.cpp" "pageRank_time.cpp" "pageRankf_time.cpp" "ELLf.cpp" "CSRf.cpp" "COOf.cpp" "puissancef.cpp" "puissance.cpp")

add_library(matrix_library ${FILES_LIST} ${HEADER_LIST})

target_include_directories(matrix_library PRIVATE ../include)

target_compile_features(matrix_library PUBLIC cxx_std_11)


