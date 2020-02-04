file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/pde900.mtx
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/test.mtx
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/well1033.mtx
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)


add_executable(app main.cpp)
target_compile_features(app PRIVATE cxx_std_17)

target_link_libraries(app PRIVATE matrix_library)
target_include_directories(app PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(app PUBLIC data)
install(TARGETS app DESTINATION bin)

add_executable(testmpi testMPI.cpp)
target_compile_features(testmpi PRIVATE cxx_std_17)

target_link_libraries(testmpi PUBLIC MPI::MPI_CXX)
target_link_libraries(testmpi PRIVATE matrix_library)
target_include_directories(testmpi PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testmpi PUBLIC data)
install(TARGETS testmpi DESTINATION bin)

add_executable(testmpi2 testMPI2.cpp)
target_compile_features(testmpi2 PRIVATE cxx_std_17)

target_link_libraries(testmpi2 PUBLIC MPI::MPI_CXX)
target_link_libraries(testmpi2 PRIVATE matrix_library)
target_include_directories(testmpi2 PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testmpi2 PUBLIC data)
install(TARGETS testmpi2 DESTINATION bin)

