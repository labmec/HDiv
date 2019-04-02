//
//  THybridizeDFN.h
//  HDiv
//
//  Created by Omar Durán on 3/23/19.
//

#ifndef THybridizeDFN_h
#define THybridizeDFN_h

#include <stdio.h>
#include <tuple>

#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzstack.h"
#include "pzintel.h"
#include "TFracture.h"
#include "TPZLagrangeMultiplier.h"
#include "pzl2projection.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "TPZMixedDarcyFlow.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"


/// This class is dedicated for conformal geometrical partitions and mixed meshes.
class THybridizeDFN : public TPZHybridizeHDiv {
    
public: /// turn to private when is complete and implement access methods
    
    /// Set of boundary material ids - boundary condition type - boundary data associated to 1D elements
    std::set<std::tuple<int,int,REAL>> m_bc_ids_1d;
    
    /// Set of boundary material ids - boundary condition type - boundary data associated to 0D elements
    std::set<std::tuple<int,int,REAL>> m_bc_ids_0d;
    
    /// List of fracture characteristics
    TPZStack<TFracture> m_fracture_data;
    
    /// Set of fractures identifiers
    std::set<int> m_fracture_ids;
    
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
    void CreateInterfaceElements(int target_dim, int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec);
    
    /// Set fracture data
    void SetFractureData(TPZStack<TFracture> & fracture_data);
    
    /// Set geometry dimension
    void SetDimension(int dimension);
    
    void LoadReferencesByDimension(TPZCompMesh * flux_cmesh, int dim);
    
    void ComputeMaterialIds(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id, int shift = 0);
    
    void InsertMaterials(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    void ApplyHibridizationOnInternalFaces(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id);
    
    TPZCompMesh * DuplicateMultiphysicsCMeshMaterials(TPZCompMesh * cmesh);
    
    TPZManVector<int,5> & ExtractActiveApproxSpaces(TPZCompMesh * cmesh);
    
    TPZManVector<TPZCompMesh *, 3> & GetMeshVector(TPZCompMesh * cmesh);
    
    void CleanUpMultiphysicsCMesh(TPZCompMesh * cmesh);
    
    void InsertMaterialsForHibridization(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    void InsertMaterialsForMixedOperatorOnFractures(int target_dim, TPZCompMesh * cmesh);
    
    void BuildMixedOperatorOnFractures(int p_order, int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    /// Construct a lagrange multiplier approximation space over the target dimension elements
    TPZCompMesh * Hybridize(TPZCompMesh * cmesh, int target_dim);
    
    TPZStack<TPZCompElSide> AnalyzeSide(int target_dim, TPZGeoEl * gel, int side);
    
    std::pair<int,int> HybridizeInSide(const TPZVec<TPZCompElSide> &candidates, TPZCompMesh * flux_cmesh, int & flux_trace_mat_id, int & lagrange_mat_id);
    
    int CreateBCGeometricalElement(const TPZCompElSide & cel_side, TPZCompMesh * flux_cmesh,int & bc_impervious_id);
    
    void ClassifyCompelSides(int target_dim, TPZCompMesh * flux_cmesh, TPZStack<std::pair<int, int>> & gel_index_and_order_lagrange_mult, int impervious_bc_id, int & flux_trace_id, int & lagrange_id);
    
    void CreareLagrangeMultiplierSpace(TPZCompMesh * p_cmesh, TPZStack<std::pair<int, int>> & gel_index_and_order_stack);
    
    void BuildMultiphysicsCMesh(int dim, TPZCompMesh * hybrid_cmesh, TPZManVector<int,5> & approx_spaces, TPZManVector<TPZCompMesh *, 3> mesh_vec);
    
};

#endif /* THybridizeDFN_h */
