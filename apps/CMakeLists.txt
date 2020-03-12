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

add_executable(testmpis testMPI2s.cpp)
target_compile_features(testmpis PRIVATE cxx_std_17)

target_link_libraries(testmpis PUBLIC MPI::MPI_CXX)
target_link_libraries(testmpis PRIVATE matrix_library)
target_include_directories(testmpis PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testmpis PUBLIC data)
install(TARGETS testmpis DESTINATION bin)

add_executable(testmpi2rf testMPI2rf.cpp)
target_compile_features(testmpi2rf PRIVATE cxx_std_17)

target_link_libraries(testmpi2rf PUBLIC MPI::MPI_CXX)
target_link_libraries(testmpi2rf PRIVATE matrix_library)
target_include_directories(testmpi2rf PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testmpi2rf PUBLIC data)
install(TARGETS testmpi2rf DESTINATION bin)

add_executable(testnorm testNorm.cpp)
target_compile_features(testnorm PRIVATE cxx_std_17)

target_link_libraries(testnorm PUBLIC MPI::MPI_CXX)
target_link_libraries(testnorm PRIVATE matrix_library)
target_include_directories(testnorm PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testnorm PUBLIC data)
install(TARGETS testnorm DESTINATION bin)

add_executable(testpuissance testPuissance.cpp)
target_compile_features(testpuissance PRIVATE cxx_std_17)

target_link_libraries(testpuissance PUBLIC MPI::MPI_CXX)
target_link_libraries(testpuissance PRIVATE matrix_library)
target_include_directories(testpuissance PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testpuissance PUBLIC data)
install(TARGETS testpuissance DESTINATION bin)

add_executable(testpagerank testPageRank.cpp)
target_compile_features(testpagerank PRIVATE cxx_std_17)

target_link_libraries(testpagerank PUBLIC MPI::MPI_CXX)
target_link_libraries(testpagerank PRIVATE matrix_library)
target_include_directories(testpagerank PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testpagerank PUBLIC data)
install(TARGETS testpagerank DESTINATION bin)

add_executable(sequentiel testSequencial.cpp)
target_compile_features(sequentiel PRIVATE cxx_std_17)

target_link_libraries(sequentiel PUBLIC MPI::MPI_CXX)
target_link_libraries(sequentiel PRIVATE matrix_library)
target_include_directories(sequentiel PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(sequentiel PUBLIC data)
install(TARGETS sequentiel DESTINATION bin)

add_executable(testmpi2time testMPI2time.cpp)
target_compile_features(testmpi2time PRIVATE cxx_std_17)

target_link_libraries(testmpi2time PUBLIC MPI::MPI_CXX)
target_link_libraries(testmpi2time PRIVATE matrix_library)
target_include_directories(testmpi2time PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testmpi2time PUBLIC data)
install(TARGETS testmpi2time DESTINATION bin)

add_executable(testpuissancetime testPuissanceTime.cpp)
target_compile_features(testpuissancetime PRIVATE cxx_std_17)

target_link_libraries(testpuissancetime PUBLIC MPI::MPI_CXX)
target_link_libraries(testpuissancetime PRIVATE matrix_library)
target_include_directories(testpuissancetime PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testpuissancetime PUBLIC data)
install(TARGETS testpuissancetime DESTINATION bin)

add_executable(testpageranktime testPageRankTime.cpp)
target_compile_features(testpageranktime PRIVATE cxx_std_17)

target_link_libraries(testpageranktime PUBLIC MPI::MPI_CXX)
target_link_libraries(testpageranktime PRIVATE matrix_library)
target_include_directories(testpageranktime PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(testpageranktime PUBLIC data)
install(TARGETS testpageranktime DESTINATION bin)