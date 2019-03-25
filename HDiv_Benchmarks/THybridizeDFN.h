//
//  THybridizeDFN.h
//  HDiv
//
//  Created by Omar Durán on 3/23/19.
//

#ifndef THybridizeDFN_h
#define THybridizeDFN_h

#include <stdio.h>
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzstack.h"
#include "pzintel.h"
#include "TFracture.h"
#include "TPZLagrangeMultiplier.h"
#include "pzl2projection.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"

class THybridizeDFN : public TPZHybridizeHDiv {
    
public: /// turn to private when is complete and implement access methods
    
    /// Set of fractures identifiers
    std::set<int> m_fracture_ids;
    
    /// List of fracture characteristics
    TPZStack<TFracture> m_fracture_data;
    
    /// stands for base geometry
    TPZGeoMesh * m_geometry;
    
    /// stands for geometry dimension
    int m_geometry_dim;
    
public:
    
    /// Default constructor
    THybridizeDFN();
    
    /// Default destructor
    ~THybridizeDFN();
    
    /// Method that duplicate and dissociate the connect belonging to the right computational element side
    std::tuple<int, int> DissociateConnects(int flux_trace_id, int lagrange_mult_id, const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> & mesh_vec);
    
    /// Create interface multiphysics elements
    void CreateInterfaceElements(int interface_id, TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
    /// Set fracture data
    void SetFractureData(TPZStack<TFracture> & fracture_data);
    
    /// Set geometry dimension
    void SetDimension(int dimension);
    
    void ComputeAndInsertMaterials(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    void ApplyHibridizationOnInternalFaces(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id);
    
    TPZCompMesh * DuplicateMultiphysicsCMeshMaterials(TPZCompMesh * cmesh);
    
    TPZManVector<int,5> & ExtractActiveApproxSpaces(TPZCompMesh * cmesh);
    
    TPZManVector<TPZCompMesh *, 3> & GetMeshVector(TPZCompMesh * cmesh);
    
    void CleanUpMultiphysicsCMesh(TPZCompMesh * cmesh);
    
    void InsertMaterialsForHibridization(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    /// Construct a lagrange multiplier approximation space over the target dimension elements
    void Hybridize(TPZCompMesh * cmesh, int target_dim);
    
};

#endif /* THybridizeDFN_h */
