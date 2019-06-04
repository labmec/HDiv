
# Configure the path.h.in file.
# This file can only use the @VARIABLE@ syntax in the file
configure_file(path.h.in ${PROJECT_SOURCE_DIR}/HDiv_Benchmarks_II/path.h @ONLY)


file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(HDivBenchmarks_II HDivBenchmarks_II.cpp ${sources} ${header})


target_link_libraries(HDivBenchmarks_II pz ${CBLAS_LIB} ${LAPACK_LIB}) 

