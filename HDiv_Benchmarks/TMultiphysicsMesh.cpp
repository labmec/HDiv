//
//  TMultiphysicsMesh.cpp
//  pz
//
//  Created by Omar Durán on 3/21/19.
//

#include "TMultiphysicsMesh.h"


TMultiphysicsMesh::TMultiphysicsMesh(){
    
    m_active_physics.Resize(0);
    m_inert_physics.Resize(0);
    m_mesh_vector.Resize(0);
    
}

TMultiphysicsMesh::TMultiphysicsMesh(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh * > & mesh_vector) : TPZCompMesh(gmesh){
    m_mesh_vector = mesh_vector;
}


void TMultiphysicsMesh::AutoBuild(){
    
    TPZCompMesh::AutoBuild();
    
    if (m_mesh_vector.size() == 0) {
        DebugStop();
    }
    
    TPZBuildMultiphysicsMesh::AddElements(m_mesh_vector, this);
    TPZBuildMultiphysicsMesh::AddConnects(m_mesh_vector, this);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(m_mesh_vector, this);
}

void TMultiphysicsMesh::SetActivePhysics(TPZVec<int> & active_physics){
    m_active_physics = active_physics;
}

void TMultiphysicsMesh::SetInertPhysics(TPZVec<int> & m_inert_physics){
    m_inert_physics = m_inert_physics;
}

TMultiphysicsMesh::TMultiphysicsMesh(const TMultiphysicsMesh &other) : TPZCompMesh(other) {
    
}

TMultiphysicsMesh & TMultiphysicsMesh::operator=(const TMultiphysicsMesh &other){
    
    if (this != & other) // prevent self-assignment
    {
        TPZCompMesh::operator=(other);
        m_active_physics    = other.m_active_physics;
        m_inert_physics     = other.m_inert_physics;
        m_mesh_vector       = other.m_mesh_vector;
    }
    return *this;
}
