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
    TPZVec<int> m_active_physics;
    
    /// Vector of inert physics
    TPZVec<int> m_inert_physics;
    
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
  
    void SetActivePhysics(TPZVec<int> & active_physics);
    
    void SetInertPhysics(TPZVec<int> & m_inert_physics);
};

#endif /* TMultiphysicsMesh_h */
