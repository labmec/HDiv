

#ifndef CreateDFN_h
#define CreateDFN_h

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
#include "TPZNormalDarcyFlow.h"


/// This class is dedicated for conformal geometrical partitions and mixed meshes.
class CreateDFN : public TPZHybridizeHDiv {
    
private:
    
    int fHDivWrapMatid =0;
    int fLagrangeInterface=0;
    int fInterfaceMatid=0;
    int fflux_resistivity_id =0;
    
    /// Set of boundary material ids - boundary condition type - boundary data associated to 2D elements
    // first : material id
    // second : boundary condition type
    // third : value of the boundary condition
    std::vector<std::tuple<int,int,REAL>> m_bc_ids_2d;
    
    /// Set of boundary material ids - boundary material ids of fractures intersections 1d
    // first : material id of the 3d boundary condition
    // second : material id of the 2d fracture boundary condition
    std::map<int,int> m_bc_ids_1d;
    
    /// Set of boundary material ids - boundary material ids of fractures intersections 0d
    // first : material id of the 3d boundary condition
    // second : material id of the 1d fracture boundary condition
    std::map<int,int> m_bc_ids_0d;
    public:
    /// List of fracture characteristics
    TPZStack<TFracture> m_fracture_data;
    
    /// Set of fractures identifiers
    // first : dimension of the fracture
    // second : set of fracture ids
    std::map<int, std::set<int>> m_fracture_ids;
    
    /// stands for base geometry
    TPZGeoMesh * m_geometry;
    

    
    /// Default constructor
    CreateDFN();
    
    /// Default destructor
    ~CreateDFN();
    
    /// Create interface multiphysics elements
    void CreateInterfaceElements(int target_dim, int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec);
    
    /// Set fracture data
    void SetFractureData(TPZStack<TFracture> & fracture_data);
    
    void LoadReferencesByDimension(TPZCompMesh * flux_cmesh, int dim);
    
    void ComputeMaterialIds(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & flux_resistivity_id, int & mp_nterface_id, int shift = 0);
    
    /// insert the given material ids in the fluxmesh.
    // insert the impervious material to seal the fractures
    // put the lagrange matid in the pressure mesh
    void InsertMaterials(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & flux_resistivity_id);
    
    
    TPZCompMesh * DuplicateMultiphysicsCMeshMaterials(TPZCompMesh * cmesh);
    
    TPZVec<int> & ExtractActiveApproxSpaces(TPZCompMesh * cmesh);
    
    TPZVec<TPZCompMesh *> & GetMeshVector(TPZCompMesh * cmesh);
    
    void CleanUpMultiphysicsCMesh(TPZCompMesh * cmesh);
    
    void InsertMaterialsForHibridization(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & flux_resistivity_id, int & mp_nterface_id);
    
    void InsertMaterialsForMixedOperatorOnFractures(int target_dim, TPZCompMesh * cmesh);
    
    /// switch the reference of the pressure computational element to the geometric element of the neighbouring fracture element
    // it will also change the material id of HDivBound element to resistivity material id
    // inserts the material objects 0f boundary conditions of fractures in the flux mesh
    void BuildFracturesAproxSpace(int p_order, int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & flux_resistivity_id, int & mp_nterface_id);
    
    void BuildMultiphysicsCMesh(int dim, TPZCompMesh * hybrid_cmesh, TPZVec<int> & approx_spaces, TPZManVector<TPZCompMesh *, 3> mesh_vec);
    
    
    /// create boundary elements of dimension target_dim-1 which lay on the boundary of the three dimensional mesh
    // target_dim : consider fracture elements of dimension target_dim
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
    
    /// group and condense the elements
    static void GroupElements(TPZMultiphysicsCompMesh *cmesh);
    
    TPZCompMesh *CreateDFNCmesh(TPZMultiphysicsCompMesh *initial_mesh);
    
    void SetPeriferalMaterialIds(int HDivWrapMatid, int LagrangeInterface, int InterfaceMatid, int flux_resistivity_id){
        fHDivWrapMatid = HDivWrapMatid;
        fLagrangeInterface =LagrangeInterface;
        fInterfaceMatid =InterfaceMatid;
        fflux_resistivity_id =flux_resistivity_id;
    }
    
    void * Hybridize(TPZMultiphysicsCompMesh *multiphysics, bool group_elements, double Lagrange_term_multiplier);
    
    void CheckODElements(TPZCompMesh *dfn_hybrid_cmesh);
    void CreateAllInterfaceElements(int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec);
    
};

#endif /* THybridizeDFN_h */
