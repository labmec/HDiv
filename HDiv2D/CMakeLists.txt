
set(HDiv_SRCS ErrorOnFaces.cpp HDiv.cpp main_ZigZagMesh.cpp ErrorOnFaces.h  )

add_executable(HDiv ${HDiv_SRCS} )

target_link_libraries(HDiv pz)




set(sourcesMixed TPZCompelDiscScaled.h TPZCompelDiscScaled.cpp meshgen.h meshgen.cpp pzmixedelasmat.h pzmixedelasmat.cpp)

set(Mymain main.cpp )

add_executable(MixedElasticitySquare ${Mymain} ${sourcesMixed})

target_link_libraries(MixedElasticitySquare pz)

set(MymainG mainGirkmann.cpp ${sourcesMixed} )

add_executable(MixedElasticityGirkmann ${MymainG})

target_link_libraries(MixedElasticityGirkmann pz)

configure_file(Girkmann.msh Girkmann.msh COPYONLY)


