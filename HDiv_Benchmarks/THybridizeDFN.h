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
    
private:

    /// Set of boundary material ids - boundary condition type - boundary data associated to 2D elements
    std::vector<std::tuple<int,int,REAL>> m_bc_ids_2d;
    
    /// Set of boundary material ids - boundary material ids of fractures intersections 1d
    std::map<int,int> m_bc_ids_1d;
    
    /// Set of boundary material ids - boundary material ids of fractures intersections 0d
    std::map<int,int> m_bc_ids_0d;
    
    /// List of fracture characteristics
    TPZStack<TFracture> m_fracture_data;
    
    /// Set of fractures identifiers
    std::set<int> m_fracture_ids;
    
    /// stands for base geometry
    TPZGeoMesh * m_geometry;
    
public:
    
    /// Default constructor
    THybridizeDFN();
    
    /// Default destructor
    ~THybridizeDFN();
    
    /// Create interface multiphysics elements
    void CreateInterfaceElements(int target_dim, int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec);
    
    /// Set fracture data
    void SetFractureData(TPZStack<TFracture> & fracture_data);
    
    void LoadReferencesByDimension(TPZCompMesh * flux_cmesh, int dim);
    
    void ComputeMaterialIds(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id, int shift = 0);
    
    void InsertMaterials(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    
    TPZCompMesh * DuplicateMultiphysicsCMeshMaterials(TPZCompMesh * cmesh);
    
    TPZManVector<int,5> & ExtractActiveApproxSpaces(TPZCompMesh * cmesh);
    
    TPZManVector<TPZCompMesh *, 3> & GetMeshVector(TPZCompMesh * cmesh);
    
    void CleanUpMultiphysicsCMesh(TPZCompMesh * cmesh);
    
    void InsertMaterialsForHibridization(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    void InsertMaterialsForMixedOperatorOnFractures(int target_dim, TPZCompMesh * cmesh);
    
    void BuildMixedOperatorOnFractures(int p_order, int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id);
    
    /// Construct a lagrange multiplier approximation spaces for a DFN
    TPZCompMesh * Hybridize(TPZCompMesh * cmesh);
    
    TPZStack<TPZCompElSide> AnalyzeSide(int target_dim, TPZGeoEl * gel, int side);
    
    std::pair<int,int> HybridizeInSide(const TPZVec<TPZCompElSide> &candidates, TPZCompMesh * flux_cmesh, int & flux_trace_mat_id, int & lagrange_mat_id);
    
    int CreateBCGeometricalElement(const TPZCompElSide & cel_side, TPZCompMesh * flux_cmesh,int & bc_impervious_id);
    
    void ClassifyCompelSides(int target_dim, TPZCompMesh * flux_cmesh, TPZStack<std::pair<int, int>> & gel_index_and_order_lagrange_mult, int & flux_trace_id, int & lagrange_id);
    
    void CreareLagrangeMultiplierSpace(TPZCompMesh * p_cmesh, TPZStack<std::pair<int, int>> & gel_index_and_order_stack);
    
    void BuildMultiphysicsCMesh(int dim, TPZCompMesh * hybrid_cmesh, TPZManVector<int,5> & approx_spaces, TPZManVector<TPZCompMesh *, 3> mesh_vec);
    
    void CreateFractureBCGeoElements(int target_dim, TPZGeoMesh * gmesh, std::set<int> bc_indexes, std::set<int> bc_frac_indexes);
    
    
    /// Set the set of boundary material ids - boundary condition type - boundary data associated to 2D elements
    void SetReservoirBoundaryData(std::vector<std::tuple<int,int,REAL>> & bc_ids_2d){
        m_bc_ids_2d = bc_ids_2d;
    }
    
    /// Set the set of boundary material ids - boundary material ids of fractures intersections 1d
    void SetMapReservoirBCToDFNBC1DIds( std::map<int,int> & bc_ids_1d){
        m_bc_ids_1d = bc_ids_1d;
    }
    
    /// Set the set of boundary material ids - boundary material ids of fractures intersections 0d
    void SetMapReservoirBCToDFNBC0DIds( std::map<int,int> & bc_ids_0d){
        m_bc_ids_0d = bc_ids_0d;
    }
    
    
};

#endif /* THybridizeDFN_h */
