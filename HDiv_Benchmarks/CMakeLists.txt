file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Hdiv3DBenchmarks ${sources} ${header})

target_link_libraries(Hdiv3DBenchmarks pz ${CBLAS_LIB} ${LAPACK_LIB}) 
