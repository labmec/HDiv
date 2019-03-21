//
//  TMultiphysicsMesh.h
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#ifndef TMultiphysicsMesh_h
#define TMultiphysicsMesh_h

#include <stdio.h>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzbuildmultiphysicsmesh.h"

class TMultiphysicsMesh : public TPZCompMesh {
    
    /// Vector of active physics
    TPZVec<int> m_active_approx_spaces;
    
    /// Vector of computational meshes
    TPZVec<TPZCompMesh * > m_mesh_vector;
    
public:
    
    /// Default constructor
    TMultiphysicsMesh();
    
    /// Constructor based on TPZGeoMesh pointer and vector of meshes
    TMultiphysicsMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > & mesh_vector);
    
    /// Copy constructor
    TMultiphysicsMesh(const TMultiphysicsMesh &other);
    
    /// Assignement constructor
    TMultiphysicsMesh & operator=(const TMultiphysicsMesh &other);
    
    /// Automatic builder for the computational mesh structure
    void AutoBuild();
  
    /// Set active approximation spaces
    void SetActiveApproxSpaces(TPZVec<int> & active_approx_spaces);
};

#endif /* TMultiphysicsMesh_h */
