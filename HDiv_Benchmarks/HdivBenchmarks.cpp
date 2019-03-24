
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZRefPatternTools.h"

#include "tpzhierarquicalgrid.h"
#include "TPZReadGIDGrid.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"
#include "tpzquadraticquad.h"
#include "tpzarc3d.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "pzfunction.h"
#include "tpzchangeel.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"

#include "TPZPrimalPoisson.h"
#include "TPZDualPoisson.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZCompMeshTools.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"

#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckmesh.h"
#include "TPZGmshReader.h"
#include "pzintel.h"

#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZLagrangeMultiplier.h"


#include "THybridzeDFN.h"
#include "pzl2projection.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif


struct SimulationCase {
    bool            IsMHMQ;
    bool            IsHybrid;
    bool            UsePardisoQ;
    bool            UseFrontalQ;
    bool            UseGmshMeshQ;
    bool            NonAffineQ;
    int             elemen_type;
    int             n_h_levels;
    int             n_p_levels;
    int             n_acc_terms;
    int             int_order;
    int             n_threads;
    int             perturbation_type;
    std::string     mesh_type;
    std::string     domain_type;
    std::string     conv_summary;
    std::string     dump_folder;
    TPZStack<int>   omega_ids;
    TPZStack<int>   gamma_ids;
    TPZStack<REAL>   permeabilities;
    TPZStack<REAL>   type;
    TPZStack<REAL>   vals;
    
    SimulationCase() : IsMHMQ(false), UsePardisoQ(true), IsHybrid(false),UseFrontalQ(false), UseGmshMeshQ(false), NonAffineQ(false), elemen_type(0), n_h_levels(0), n_p_levels(1), n_acc_terms(0), int_order(1), n_threads(0),perturbation_type(0), mesh_type(""), domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),gamma_ids(), permeabilities(), type(), vals()
    {
        
    }
    
    SimulationCase(const SimulationCase &copy) :IsHybrid(copy.IsHybrid) ,IsMHMQ(copy.IsMHMQ), UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
    UseGmshMeshQ(copy.UseGmshMeshQ), NonAffineQ(copy.NonAffineQ), elemen_type(copy.elemen_type), n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), n_acc_terms(copy.n_acc_terms), int_order(copy.int_order),n_threads(copy.n_threads),perturbation_type(copy.perturbation_type), mesh_type(copy.mesh_type), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
    dump_folder(copy.dump_folder), omega_ids(copy.omega_ids), gamma_ids(copy.gamma_ids),
    permeabilities(copy.permeabilities),
    type(copy.type),
    vals(copy.vals)
    {
        
    }
    
    SimulationCase &operator=(const SimulationCase &copy)
    {
        
        IsMHMQ = copy.IsMHMQ;
        IsHybrid = copy.IsHybrid;
        UsePardisoQ = copy.UsePardisoQ;
        UseFrontalQ = copy.UseFrontalQ;
        UseGmshMeshQ = copy.UseGmshMeshQ;
        NonAffineQ = copy.NonAffineQ;
        elemen_type = copy.elemen_type;
        n_h_levels = copy.n_h_levels;
        n_p_levels = copy.n_p_levels;
        n_acc_terms = copy.n_acc_terms;
        int_order = copy.int_order;
        n_threads = copy.n_threads;
        perturbation_type = copy.perturbation_type;
        mesh_type = copy.mesh_type;
        domain_type = copy.domain_type;
        conv_summary = copy.conv_summary;
        dump_folder = copy.dump_folder;
        omega_ids = copy.omega_ids;
        gamma_ids = copy.gamma_ids;
        permeabilities=copy.permeabilities;
        type=copy.type;
        vals=copy.vals;
        return *this;
    }
};

void InsertFrac(TPZGeoMesh *gmesh, TPZFMatrix<REAL> corners, int matids );
TPZGeoMesh * case2mesh();
TPZGeoMesh * PrettyCubemesh();
TPZCompMesh * CMeshMixed(TPZGeoMesh * geometry, int p, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);
TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, int p, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);
TPZCompMesh * FluxMesh(TPZGeoMesh * gmesh, int order, SimulationCase sim);
TPZCompMesh * PressureMesh(TPZGeoMesh * gmesh, int order,SimulationCase sim);
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
void forcing(const TPZVec<REAL> &p, TPZVec<STATE> &f);
void SeparateConnectsByFracId(TPZCompMesh * mixed_cmesh,int fracid);
void InsertFractureMaterial(TPZCompMesh *cmesh);
void FractureTest();

/// Executes case 1
void Case_1();

/// Executes case 2
void Case_2();

/// Executes cube
void Pretty_cube();

int main(){
    Pretty_cube();
//    Case_1();
//     Case_2();
//    FractureTest();
}


/// Executes cube
void Pretty_cube(){
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.omega_ids.push_back(1);
    sim.permeabilities.push_back(1.0);
    
    sim.gamma_ids.push_back(-1);
    sim.gamma_ids.push_back(-2);
    sim.gamma_ids.push_back(-3);
    sim.gamma_ids.push_back(-4);
    sim.gamma_ids.push_back(-5);
    sim.gamma_ids.push_back(-6);
    
    
    //    D = 0;
    //    N = 1;
    sim.type.push_back(0);
    sim.type.push_back(1);
    sim.type.push_back(0);
    sim.type.push_back(1);
    sim.type.push_back(1);
    sim.type.push_back(1);
    
    sim.vals.push_back(2.0);
    sim.vals.push_back(0.0);
    sim.vals.push_back(1.0);
    sim.vals.push_back(0.0);
    sim.vals.push_back(0.0);
    sim.vals.push_back(0.0);
    
    
    TPZGmshReader Geometry;
    TPZGeoMesh *gmesh = PrettyCubemesh();
    std::ofstream file("geometry_cube.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    
    std::ofstream file_txt("geometry_cube_base.txt");
    gmesh->Print(file_txt);
    
    int p_order = 1;
    TPZVec<TPZCompMesh *> meshvec;
    TPZCompMesh *cmixedmesh = NULL;
    cmixedmesh = MPCMeshMixed(gmesh, p_order, sim, meshvec);
    std::ofstream filemixed("mixedMesh.txt");
    cmixedmesh->Print(filemixed);
    
    TPZCompMesh *cmeshm =NULL;
    if(sim.IsHybrid){
        TPZCompMesh * cmesh_m_Hybrid;
        
        TPZHybridizeHDiv hybridizer;
        
        bool old_version_Q = false;
        
        if(old_version_Q)
        {
            TPZManVector<TPZCompMesh*, 2> meshvector_Hybrid(2);
            std::ofstream file_txt("geometry_cube_base_before_Hybridize.txt");
            cmixedmesh->Reference()->Print(file_txt);
            tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(cmixedmesh, meshvec, true, -1.);
            cmesh_m_Hybrid->InitializeBlock();
            
            std::ofstream file_geo_hybrid("geometry_cube_hybrid_old.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(cmesh_m_Hybrid->Reference(), file_geo_hybrid);
            
            std::ofstream file_hybrid_mixed("Hybrid_mixed_cmesh_old.txt");
            cmesh_m_Hybrid->Print(file_hybrid_mixed);
            {
                std::ofstream file_hybrid_mixed_q("Hybrid_mixed_cmesh_old_q.txt");
                meshvector_Hybrid[0]->Print(file_hybrid_mixed_q);
                std::ofstream file_hybrid_mixed_p("Hybrid_mixed_cmesh_old_p.txt");
                meshvector_Hybrid[1]->Print(file_hybrid_mixed_p);
            }
            return;
        }
        
        
        THybridzeDFN dfn_hybridzer;
        int target_dim = 3;
        TPZStack<int> fracture_ids;
        fracture_ids.Push(100);
        dfn_hybridzer.Hybridize(cmixedmesh,target_dim,fracture_ids);
        
        std::set<int> m_fracture_ids;
        m_fracture_ids.insert(100);
        {
            int material_id = 1; /// Material ids being hybridized
            TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmixedmesh);
            if (!mp_cmesh) {
                DebugStop();
            }
            
//            int n_meshes = mp_cmesh->MeshVector().size();
//            TPZManVector<TPZCompMesh *, 2> dfn_mixed_mesh_vec(n_meshes, 0);
//            for (int i = 0; i <  n_meshes; i++) {
//                dfn_mixed_mesh_vec[i] =  mp_cmesh->MeshVector()[i]->Clone();
//            }
            
            TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = mp_cmesh->MeshVector();
            
            int lagrange_mult_id; /// compute a proper LM id
            int flux_trace_id; /// compute a proper flux_trace id
            int truly_lagrange_mult_id;
            int n_state = 1;
            
            /// Computes maximum material identifier
            {
                int maxMatId = std::numeric_limits<int>::min();
                for (auto &mesh : dfn_mixed_mesh_vec) {
                    for (auto &mat : mesh->MaterialVec()) {
                        maxMatId = std::max(maxMatId, mat.first);
                    }
                }
                if (maxMatId == std::numeric_limits<int>::min()) {
                    maxMatId = 0;
                }
                flux_trace_id = maxMatId + 1;
                lagrange_mult_id = maxMatId + 2;
                truly_lagrange_mult_id = maxMatId + 3;
            }
            
            /// Insert LM and flux_trace materials
            {
                TPZFNMatrix<1, STATE> xk(n_state, n_state, 0.), xb(n_state, n_state, 0.), xc(n_state, n_state, 0.), xf(n_state, 1, 0.);
                TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
                if (!flux_cmesh) {
                    DebugStop();
                }
                if (!flux_cmesh->FindMaterial(flux_trace_id)) {

                    if (target_dim == 2) {
                        auto flux_trace_mat = new TPZMat1dLin(flux_trace_id);
                        flux_trace_mat->SetMaterial(xk, xb, xc, xf);
                        flux_cmesh->InsertMaterialObject(flux_trace_mat);
                    }else if(target_dim == 3){
                        auto flux_trace_mat = new TPZMat2dLin(flux_trace_id);
                        flux_trace_mat->SetMaterial(xk, xc, xf);
                        flux_cmesh->InsertMaterialObject(flux_trace_mat);
                    }
                }

                TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
                if (!pressure_cmesh) {
                    DebugStop();
                }
                if (!pressure_cmesh->FindMaterial(lagrange_mult_id)) {

                    if (target_dim == 2) {
                        auto lagrange_mult_mat = new TPZMat1dLin(lagrange_mult_id);
                        lagrange_mult_mat->SetMaterial(xk, xb, xc, xf);
                        pressure_cmesh->InsertMaterialObject(lagrange_mult_mat);
                    }else if(target_dim == 3){
                        auto lagrange_mult_mat = new TPZMat2dLin(lagrange_mult_id);
                        lagrange_mult_mat->SetMaterial(xk, xc, xf);
                        pressure_cmesh->InsertMaterialObject(lagrange_mult_mat);
                    }
                }
                
            }
            
            /// HybridizeInternalFaces
            {
                TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
                TPZGeoMesh * gmesh = flux_cmesh->Reference();
                int dim = gmesh->Dimension();
                gmesh->ResetReference();
                flux_cmesh->LoadReferences();
                int64_t nel = flux_cmesh->NElements();
                std::list<std::tuple<int64_t, int> > pressures;
                for (int64_t el = 0; el < nel; el++) {
                    TPZCompEl *cel = flux_cmesh->Element(el);
                    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                    if (!intel || (intel->Reference()->Dimension() != dim && material_id != cel->Material()->Id())) {
                        continue;
                    }
                    
                    // loop over the side of dimension dim-1
                    TPZGeoEl *gel = intel->Reference();
                    for (int side = gel->NCornerNodes(); side < gel->NSides() - 1; side++) {
                        if (gel->SideDimension(side) != dim - 1) {
                            continue;
                        }

                        TPZCompElSide celside(intel, side);
                        TPZCompElSide neighcomp = hybridizer.RightElement(intel, side); //it makes sense
                        if (neighcomp) {
                           pressures.push_back(dfn_hybridzer.DissociateConnects(flux_trace_id,lagrange_mult_id,celside, neighcomp, dfn_mixed_mesh_vec));
                        }
                    }
                }
                flux_cmesh->InitializeBlock();
                flux_cmesh->ComputeNodElCon();
                
                TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
                gmesh->ResetReference();
                pressure_cmesh->SetDimModel(gmesh->Dimension()-1);
                for (auto pindex : pressures) {
                    int64_t elindex;
                    int order;
                    std::tie(elindex, order) = pindex;
                    TPZGeoEl *gel = gmesh->Element(elindex);
                    int64_t celindex;
                    TPZCompEl *cel = pressure_cmesh->ApproxSpace().CreateCompEl(gel, *pressure_cmesh, celindex);
                    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
                    TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
                    if (intel){
                        intel->PRefine(order);
                        //            intel->SetSideOrder(gel->NSides() - 1, order);
                    } else if (intelDisc) {
                        intelDisc->SetDegree(order);
                        intelDisc->SetTrueUseQsiEta();
                    } else {
                        DebugStop();
                    }
                    int n_connects = cel->NConnects();
                    for (int i = 0; i < n_connects; ++i) {
                        cel->Connect(i).SetLagrangeMultiplier(2);
                    }
                    gel->ResetReference();
                }
                pressure_cmesh->InitializeBlock();
                pressure_cmesh->SetDimModel(gmesh->Dimension());
                
            }
            
            /// ResconstrucMultiphysicsCMesh
            
            TPZGeoMesh *gmesh = mp_cmesh->Reference();
            TPZMultiphysicsCompMesh * dfn_hybrid_cmesh = new TPZMultiphysicsCompMesh(gmesh);
            TPZManVector<int,5> & active_approx_spaces = mp_cmesh->GetActiveApproximationSpaces();
            {
                
                /// copy already inserted materials
                {
                    cmixedmesh->CopyMaterials(*dfn_hybrid_cmesh);
                }
                
                /// Cleanup mp_cmesh
                {
                    cmixedmesh->CleanUp();
                    cmixedmesh->SetReference(NULL);
                }
                
                /// Insert fractures material objects
                {
                    
                    TPZFNMatrix<1, STATE> xk(n_state, n_state, 0.), xb(n_state, n_state, 0.), xc(n_state, n_state, 0.), xf(n_state, 1, 0.);
                    
                    if (!dfn_hybrid_cmesh) {
                        DebugStop();
                    }
                    if (!dfn_hybrid_cmesh->FindMaterial(flux_trace_id)) {
                        
                        if (target_dim == 2) {
                            auto flux_trace_mat = new TPZMat1dLin(flux_trace_id);
                            flux_trace_mat->SetMaterial(xk, xb, xc, xf);
                            dfn_hybrid_cmesh->InsertMaterialObject(flux_trace_mat);
                        }else if(target_dim == 3){
                            auto flux_trace_mat = new TPZMat2dLin(flux_trace_id);
                            flux_trace_mat->SetMaterial(xk, xc, xf);
                            dfn_hybrid_cmesh->InsertMaterialObject(flux_trace_mat);
                        }
                    }
                    
                    if (!dfn_hybrid_cmesh->FindMaterial(lagrange_mult_id)) {
                        
                        if (target_dim == 2) {
                            auto lagrange_mult_mat = new TPZMat1dLin(lagrange_mult_id);
                            lagrange_mult_mat->SetMaterial(xk, xb, xc, xf);
                            dfn_hybrid_cmesh->InsertMaterialObject(lagrange_mult_mat);
                        }else if(target_dim == 3){
                            auto lagrange_mult_mat = new TPZMat2dLin(lagrange_mult_id);
                            lagrange_mult_mat->SetMaterial(xk, xc, xf);
                            dfn_hybrid_cmesh->InsertMaterialObject(lagrange_mult_mat);
                        }
                    }
                    
                    if (!dfn_hybrid_cmesh->FindMaterial(truly_lagrange_mult_id)) {
                        auto lagrange_multiplier = new TPZLagrangeMultiplier(truly_lagrange_mult_id, target_dim - 1, n_state);
                        lagrange_multiplier->SetMultiplier(1.0);
                        dfn_hybrid_cmesh->InsertMaterialObject(lagrange_multiplier);
                    }
                }
                dfn_hybrid_cmesh->SetDimModel(target_dim);
                dfn_hybrid_cmesh->BuildMultiphysicsSpace(active_approx_spaces, dfn_mixed_mesh_vec);
                
                dfn_hybridzer.CreateInterfaceElements(truly_lagrange_mult_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
                
                dfn_hybrid_cmesh->InitializeBlock();
                cmeshm=dfn_hybrid_cmesh;
                
                
                /// Change lagrange multipliers identifiers to fractures identifiers
                {
//                    TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
                    dfn_hybrid_cmesh->Reference()->BuildConnectivity();
                    dfn_hybrid_cmesh->Reference()->ResetReference();
                    dfn_hybrid_cmesh->LoadReferences();
                    TPZGeoMesh  * gmesh = dfn_hybrid_cmesh->Reference();
                    int64_t n_cel = dfn_hybrid_cmesh->NElements();
                    for (auto gel : gmesh->ElementVec()) {
                
                        if (!gel) {
                            DebugStop();
                        }
                        
                        int gel_mat_id = gel->MaterialId();
                        std::set<int>::iterator it = m_fracture_ids.find(gel_mat_id);
                        bool is_not_fracture_Q = it == m_fracture_ids.end();
                        if (is_not_fracture_Q) {
                            continue;
                        }
                        
                        TPZCompEl * cel = gel->Reference();
                        if (!cel) {
                            DebugStop();
                        }
                        
                        
                        if (gel->MaterialId() == lagrange_mult_id) {
                            int side = gel->NSides() - 1;
                            TPZGeoElSide gel_side(gel,side);
                            TPZStack<TPZGeoElSide> gelsides;
                            gel_side.AllNeighbours(gelsides);
                            
                            for (auto iside : gelsides) {
                                if(iside.Element()->MaterialId() == material_id){
                                    int side = iside.Side();
                                    TPZGeoEl * gel = iside.Element();
                                    TPZGeoElSide gel_side(gel,side);
                                    TPZStack<TPZGeoElSide> vol_gelsides;
                                    gel_side.AllNeighbours(vol_gelsides);
                                    int aka = 0;
                                }
                            }
                            
                        }
                        
                        
                    }
                    
                }
                
                

                
                {
                    std::ofstream file_hybrid_mixed_q("Hybrid_mixed_cmesh_q.txt");
                    dfn_hybrid_cmesh->MeshVector()[0]->Print(file_hybrid_mixed_q);
                    std::ofstream file_hybrid_mixed_p("Hybrid_mixed_cmesh_p.txt");
                    dfn_hybrid_cmesh->MeshVector()[1]->Print(file_hybrid_mixed_p);
                }
                
            }

        }

    }
    else{
        cmeshm=cmixedmesh;
        
    }
    
    std::ofstream file_geo_hybrid("geometry_cube_hybrid.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(cmeshm->Reference(), file_geo_hybrid);
    
    std::ofstream file_geo_hybrid_txt("geometry_cube_hybrid.txt");
    cmeshm->Reference()->Print(file_geo_hybrid_txt);
    
    std::ofstream file_hybrid_mixed("Hybrid_mixed_cmesh.txt");
    cmeshm->Print(file_hybrid_mixed);
    
    TPZAnalysis *an = CreateAnalysis(cmeshm, sim);
    
    std::cout << "Assembly neq = " << cmeshm->NEquations() << std::endl;
    an->Assemble();
    
    std::cout << "Solution of the system" << std::endl;
    an->Solve();
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    
    int div = 0;
    std::string fileresult("cube.vtk");
    an->DefineGraphMesh(3,scalnames,vecnames,fileresult);
    an->PostProcess(div,3);
    
}


void Case_1(){
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=false;
    sim.omega_ids.push_back(1);
    sim.omega_ids.push_back(2);
    sim.permeabilities.push_back(1.0);
    sim.permeabilities.push_back(0.1);
    
    sim.gamma_ids.push_back(4);
    sim.gamma_ids.push_back(5);
    sim.gamma_ids.push_back(6);
    
    //    d = 0;
    //    n = 1;
    sim.type.push_back(1);
    sim.type.push_back(0);
    sim.type.push_back(0);
    
    sim.vals.push_back(0.0);
    sim.vals.push_back(4.0);
    sim.vals.push_back(1.0);
    
    
    
    TPZGmshReader Geometry;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    Geometry.SetFormatVersion("3");
    gmesh = Geometry.GeometricGmshMesh("case_1_1k.msh");
    Geometry.PrintPartitionSummary(std::cout);
    std::ofstream file("geometry_case_1.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);

    TPZVec<TPZCompMesh *> meshvec;
    TPZCompMesh *cmixedmesh = NULL;
    cmixedmesh = MPCMeshMixed(gmesh, 1, sim, meshvec);
    std::ofstream filemixed("mixedMesh.txt");
    cmixedmesh->Print(filemixed);
    
    TPZCompMesh *cmeshm =NULL;
    if(sim.IsHybrid){
        TPZCompMesh * cmesh_m_Hybrid;
        TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
        TPZHybridizeHDiv hybridizer;
        tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(cmixedmesh, meshvec, true, -1.);
        cmesh_m_Hybrid->InitializeBlock();
        cmeshm=cmesh_m_Hybrid;
    }
    else{
        cmeshm=cmixedmesh;
        
    }
    
    TPZAnalysis *an = CreateAnalysis(cmeshm, sim);
    
    std::cout << "Assembly neq = " << cmeshm->NEquations() << std::endl;
    an->Assemble();
    
    std::cout << "Solution of the system" << std::endl;
    an->Solve();
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    
    int div = 0;
    std::string fileresult("case_1_1k.vtk");
    an->DefineGraphMesh(3,scalnames,vecnames,fileresult);
    an->PostProcess(div,3);
    
}

void Case_2(){
    
    SimulationCase sim;
    
    sim.UsePardisoQ=true;
    sim.IsHybrid=false; /// For now testing without hybrid mesh.
    sim.omega_ids.push_back(1);
    sim.omega_ids.push_back(2);
    sim.permeabilities.push_back(1.0);
    sim.permeabilities.push_back(0.1);
    
    sim.gamma_ids.push_back(-1);
    sim.gamma_ids.push_back(-2);
    sim.gamma_ids.push_back(-3);
    sim.type.push_back(1);
    sim.type.push_back(0);
    sim.type.push_back(1);
    sim.vals.push_back(0.0);
    sim.vals.push_back(1.0);
    sim.vals.push_back(-1.0);
    
    TPZGeoMesh *gmesh = case2mesh();
    std::ofstream file("geometry_case_2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file, true);
    
    TPZVec<TPZCompMesh *> meshvec;
    TPZCompMesh *cmixedmesh = NULL;
    cmixedmesh = MPCMeshMixed(gmesh, 1, sim, meshvec);
    std::ofstream filemixed("mixedMesh.txt");
    cmixedmesh->Print(filemixed);
    
    TPZCompMesh *cmeshm =NULL;
    if(sim.IsHybrid){
        TPZCompMesh * cmesh_m_Hybrid;
        TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
        TPZHybridizeHDiv hybridizer;
        tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(cmixedmesh, meshvec, true, -1.);
        cmesh_m_Hybrid->InitializeBlock();
        cmeshm=cmesh_m_Hybrid;
    }
    else{
        cmeshm=cmixedmesh;
    }
    
    TPZAnalysis *an = CreateAnalysis(cmeshm, sim);
    std::cout << "Assembly neq = " << cmeshm->NEquations() << std::endl;
    an->Assemble();
    
    std::cout << "Solution of the system" << std::endl;
    an->Solve();
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    
    int div = 0;
    std::string fileresult("case_2.vtk");
    an->DefineGraphMesh(3,scalnames,vecnames,fileresult);
    an->PostProcess(div,3);
}

void InsertFrac(TPZGeoMesh *gmesh, TPZFMatrix<REAL> corners, int matid){
    
    //Set Frac_1
    TPZManVector<REAL,3> co(3,0.0);
    co[0] = corners(0,0);
    co[1] = corners(1,0);
    co[2] = corners(2,0);
    
    TPZGeoNode * node = gmesh->FindNode(co);
    int index_1 = gmesh->NodeIndex(node);
    
    co[0] = corners(0,1);
    co[1] = corners(1,1);
    co[2] = corners(2,1);
    
    node = gmesh->FindNode(co);
    int index_2 = gmesh->NodeIndex(node);
    
    co[0] = corners(0,2);
    co[1] = corners(1,2);
    co[2] = corners(2,2);
    
    node = gmesh->FindNode(co);
    int index_3 = gmesh->NodeIndex(node);
    
    co[0] = corners(0,3);
    co[1] = corners(1,3);
    co[2] = corners(2,3);
    
    node = gmesh->FindNode(co);
    int index_4 = gmesh->NodeIndex(node);
    TPZVec<int64_t> cords(4);
    
    cords[0]=index_1;
    cords[1]=index_2;
    cords[2]=index_3;
    cords[3]=index_4;
    
    int64_t Nels = gmesh->NElements();
    gmesh->CreateGeoElement(EQuadrilateral, cords, matid, Nels);
    
    gmesh->BuildConnectivity();
    
    TPZGeoEl * frac_gel = gmesh->Element(Nels);
    TPZVec<TPZGeoEl *> sons;
    frac_gel->Divide(sons);
}
TPZGeoMesh * case2mesh(){
    // Creating the Geo mesh
    int dimel=1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.0);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,dimel);
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    gengrid.SetRefpatternElements(true);
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -1);
    gengrid.SetBC(gmesh, 6, -1);
    gengrid.SetBC(gmesh, 7, -1);
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = 3;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    //
    
    TPZExtendGridDimension extend(gmesh,1.0/8);
    extend.SetElType(1);
    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(8,-1,-1);
    gmesh3d->BuildConnectivity();
    
    
    int nels = gmesh3d->NElements();
    for (int iel =0; iel<nels; iel++){
        TPZGeoEl * gel = gmesh3d->Element(iel);
        if(!gel){continue;}
        if (gel->Dimension() == 1){continue;}
        TPZManVector<REAL,3> qsi(gel->Dimension()), xCenter(3,0.);
        gel->CenterPoint(gel->NSides()-1, qsi);
        gel->X(qsi,xCenter);
        
        if (gel->Dimension()==3) {
            
            if (xCenter[0]>0.5 && xCenter[1]<0.5){
                gel->SetMaterialId(2);
            }
            
            if ((xCenter[0]>0.625 && xCenter[0]<0.75 ) && (xCenter[1]>0.5 && xCenter[1]<0.625)&&(xCenter[2]>0.5 && xCenter[2]<0.75)){
                gel->SetMaterialId(2);
            }
            
            if (xCenter[0]>0.75 && (xCenter[1]>0.5 && xCenter[1]<0.75)&&(xCenter[2]>0.5)){
                gel->SetMaterialId(2);
            }
        }
        if (gel->Dimension()==2) {
            
            
            if (xCenter[0]>0.875 && xCenter[1]>0.875 && xCenter[2]>0.875 ){
                gel->SetMaterialId(-2);
            }
            
            if (xCenter[0]<0.25 && xCenter[1]<0.25 && xCenter[2]<0.25 ){
                gel->SetMaterialId(-3);
            }
        }
    }
    
    TPZFMatrix<REAL> frac1(3,4,0.0);
    
    frac1(0,0)=0.5;
    frac1(1,0)=0.0;
    frac1(2,0)=0.0;
    
    frac1(0,1)=0.5;
    frac1(1,1)=0.0;
    frac1(2,1)=1.0;
    
    frac1(0,2)=0.5;
    frac1(1,2)=1.0;
    frac1(2,2)=1.0;
    
    frac1(0,3)=0.5;
    frac1(1,3)=1.0;
    frac1(2,3)=0.0;
    
    TPZFMatrix<REAL> frac2(3,4,0.0);
    
    frac2(0,0) = 1.0;
    frac2(1,0) = 0.5;
    frac2(2,0) = 0.0;
    
    frac2(0,1) = 1.0;
    frac2(1,1) = 0.5;
    frac2(2,1) = 1.0;
    
    frac2(0,2) = 0.0;
    frac2(1,2) = 0.5;
    frac2(2,2) = 1.0;
    
    frac2(0,3) = 0.0;
    frac2(1,3) = 0.5;
    frac2(2,3) = 0.0;
    
    //frac3
    TPZFMatrix<REAL> frac3(3,4,0.0);
    frac3(0,0) = 0.75;
    frac3(1,0) = 0.5;
    frac3(2,0) = 0.5;
    
    frac3(0,1) = 0.75;
    frac3(1,1) = 0.5;
    frac3(2,1) = 1.0;
    
    frac3(0,2) = 0.75;
    frac3(1,2) = 1.0;
    frac3(2,2) = 1.0;
    
    frac3(0,3) = 0.75;
    frac3(1,3) = 1.0;
    frac3(2,3) = 0.5;
    
    
    //Set Frac_4
    TPZFMatrix<REAL> frac4(3,4,0.0);
    frac4(0,0) = 1.0;
    frac4(1,0) = 0.75;
    frac4(2,0) = 0.5;
    
    frac4(0,1) = 1.0;
    frac4(1,1) = 0.75;
    frac4(2,1) = 1.0;
    
    frac4(0,2) = 0.5;
    frac4(1,2)= 0.75;
    frac4(2,2) = 1.0;
    
    frac4(0,3)= 0.50;
    frac4(1,3) = 0.75;
    frac4(2,3) = 0.50;
    
    
    //Set Frac_5
    TPZFMatrix<REAL> frac5(3,4,0.0);
    frac5(0,0) = 0.675;
    frac5(1,0) = 0.5;
    frac5(2,0) = 0.5;
    
    frac5(0,1) = 0.675;
    frac5(1,1) = 0.5;
    frac5(2,1) = 0.75;
    
    frac5(0,2) = 0.675;
    frac5(1,2) = 0.75;
    frac5(2,2) = 0.75;
    
    frac5(0,3) = 0.675;
    frac5(1,3) = 0.75;
    frac5(2,3) = 0.5;
    
    
    //Set Frac_6
    TPZFMatrix<REAL> frac6(3,4,0.0);
    frac6(0,0) = 0.75;
    frac6(1,0) = 0.675;
    frac6(2,0) = 0.5;
    
    frac6(0,1) = 0.75;
    frac6(1,1) = 0.675;
    frac6(2,1) = 0.75;
    
    frac6(0,2) = 0.5;
    frac6(1,2) = 0.675;
    frac6(2,2) = 0.75;
    
    frac6(0,3) = 0.50;
    frac6(1,3) = 0.675;
    frac6(2,3)= 0.50;
    
    TPZFMatrix<REAL> frac7(3,4,0.0);
    
    frac7(0,0) = 0.0;
    frac7(1,0) = 0.0;
    frac7(2,0) = 0.5;
    
    frac7(0,1) = 1.0;
    frac7(1,1) = 0.0;
    frac7(2,1) = 0.5;
    
    frac7(0,2) = 1.0;
    frac7(1,2) = 1.0;
    frac7(2,2) = 0.5;
    
    frac7(0,3) = 0.0;
    frac7(1,3) = 1.0;
    frac7(2,3)= 0.5;
    
    
    TPZFMatrix<REAL> frac8(3,4,0.0);
    frac8(0,0) = 0.75;
    frac8(1,0) = 0.5;
    frac8(2,0) = 0.675;
    
    frac8(0,1) = 0.75;
    frac8(1,1) = 0.75;
    frac8(2,1) = 0.675;
    
    frac8(0,2) = 0.5;
    frac8(1,2) = 0.75;
    frac8(2,2) = 0.675;
    
    frac8(0,3) = 0.50;
    frac8(1,3) = 0.5;
    frac8(2,3)= 0.675;
    
    //Set Frac_4
    TPZFMatrix<REAL> frac9(3,4,0.0);
    frac9(0,0) = 1.0;
    frac9(1,0) = 0.5;
    frac9(2,0) = 0.75;
    
    frac9(0,1) = 1.0;
    frac9(1,1) = 1.0;
    frac9(2,1) = 0.75;
    
    frac9(0,2) = 0.5;
    frac9(1,2)= 1.0;
    frac9(2,2) = 0.75;
    
    frac9(0,3)= 0.50;
    frac9(1,3) = 0.50;
    frac9(2,3) = 0.75;
    
    //  InsertFrac(mesh, corners, matid)
    InsertFrac(gmesh3d,frac1,3);
    InsertFrac(gmesh3d,frac2,3);
    InsertFrac(gmesh3d,frac3,3);
    InsertFrac(gmesh3d,frac4,3);
    InsertFrac(gmesh3d,frac5,3);
    InsertFrac(gmesh3d,frac6,3);
    InsertFrac(gmesh3d,frac7,3);
    InsertFrac(gmesh3d,frac8,3);
    InsertFrac(gmesh3d,frac9,3);
    return gmesh3d;
}

TPZGeoMesh * PrettyCubemesh(){
    
    // Creating the Geo mesh
    int dimel=1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.0);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,dimel);
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    gengrid.SetRefpatternElements(true);
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    
    TPZVec<TPZGeoEl *> sons;
    const int nref = 1;
    for (int iref = 0; iref < nref; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if(!gel->HasSubElement())
            {
                gel->Divide(sons);
            }
        }
    }
    
    TPZExtendGridDimension extend(gmesh,1.0/2);
    extend.SetElType(1);
    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(2,-5,-6);
    gmesh3d->BuildConnectivity();
    
    TPZFMatrix<REAL> frac1(3,4,0.0);
    
    frac1(0,0)=0.5;
    frac1(1,0)=0.0;
    frac1(2,0)=0.0;
    
    frac1(0,1)=0.5;
    frac1(1,1)=0.0;
    frac1(2,1)=1.0;
    
    frac1(0,2)=0.5;
    frac1(1,2)=1.0;
    frac1(2,2)=1.0;
    
    frac1(0,3)=0.5;
    frac1(1,3)=1.0;
    frac1(2,3)=0.0;
    
    TPZFMatrix<REAL> frac2(3,4,0.0);
    
    frac2(0,0) = 1.0;
    frac2(1,0) = 0.5;
    frac2(2,0) = 0.0;
    
    frac2(0,1) = 1.0;
    frac2(1,1) = 0.5;
    frac2(2,1) = 1.0;
    
    frac2(0,2) = 0.0;
    frac2(1,2) = 0.5;
    frac2(2,2) = 1.0;
    
    frac2(0,3) = 0.0;
    frac2(1,3) = 0.5;
    frac2(2,3) = 0.0;
    
    //frac3
    TPZFMatrix<REAL> frac3(3,4,0.0);
    frac3(0,0) = 0.75;
    frac3(1,0) = 0.5;
    frac3(2,0) = 0.5;
    
    frac3(0,1) = 0.75;
    frac3(1,1) = 0.5;
    frac3(2,1) = 1.0;
    
    frac3(0,2) = 0.75;
    frac3(1,2) = 1.0;
    frac3(2,2) = 1.0;
    
    frac3(0,3) = 0.75;
    frac3(1,3) = 1.0;
    frac3(2,3) = 0.5;
    
    
    //Set Frac_4
    TPZFMatrix<REAL> frac4(3,4,0.0);
    frac4(0,0) = 1.0;
    frac4(1,0) = 0.75;
    frac4(2,0) = 0.5;
    
    frac4(0,1) = 1.0;
    frac4(1,1) = 0.75;
    frac4(2,1) = 1.0;
    
    frac4(0,2) = 0.5;
    frac4(1,2)= 0.75;
    frac4(2,2) = 1.0;
    
    frac4(0,3)= 0.50;
    frac4(1,3) = 0.75;
    frac4(2,3) = 0.50;
    
    
    //Set Frac_5
    TPZFMatrix<REAL> frac5(3,4,0.0);
    frac5(0,0) = 0.675;
    frac5(1,0) = 0.5;
    frac5(2,0) = 0.5;
    
    frac5(0,1) = 0.675;
    frac5(1,1) = 0.5;
    frac5(2,1) = 0.75;
    
    frac5(0,2) = 0.675;
    frac5(1,2) = 0.75;
    frac5(2,2) = 0.75;
    
    frac5(0,3) = 0.675;
    frac5(1,3) = 0.75;
    frac5(2,3) = 0.5;
    
    
    //Set Frac_6
    TPZFMatrix<REAL> frac6(3,4,0.0);
    frac6(0,0) = 0.75;
    frac6(1,0) = 0.675;
    frac6(2,0) = 0.5;
    
    frac6(0,1) = 0.75;
    frac6(1,1) = 0.675;
    frac6(2,1) = 0.75;
    
    frac6(0,2) = 0.5;
    frac6(1,2) = 0.675;
    frac6(2,2) = 0.75;
    
    frac6(0,3) = 0.50;
    frac6(1,3) = 0.675;
    frac6(2,3)= 0.50;
    
    TPZFMatrix<REAL> frac7(3,4,0.0);
    
    frac7(0,0) = 0.0;
    frac7(1,0) = 0.0;
    frac7(2,0) = 0.5;
    
    frac7(0,1) = 1.0;
    frac7(1,1) = 0.0;
    frac7(2,1) = 0.5;
    
    frac7(0,2) = 1.0;
    frac7(1,2) = 1.0;
    frac7(2,2) = 0.5;
    
    frac7(0,3) = 0.0;
    frac7(1,3) = 1.0;
    frac7(2,3)= 0.5;
    
    
    TPZFMatrix<REAL> frac8(3,4,0.0);
    frac8(0,0) = 0.75;
    frac8(1,0) = 0.5;
    frac8(2,0) = 0.675;
    
    frac8(0,1) = 0.75;
    frac8(1,1) = 0.75;
    frac8(2,1) = 0.675;
    
    frac8(0,2) = 0.5;
    frac8(1,2) = 0.75;
    frac8(2,2) = 0.675;
    
    frac8(0,3) = 0.50;
    frac8(1,3) = 0.5;
    frac8(2,3)= 0.675;
    
    //Set Frac_4
    TPZFMatrix<REAL> frac9(3,4,0.0);
    frac9(0,0) = 1.0;
    frac9(1,0) = 0.5;
    frac9(2,0) = 0.75;
    
    frac9(0,1) = 1.0;
    frac9(1,1) = 1.0;
    frac9(2,1) = 0.75;
    
    frac9(0,2) = 0.5;
    frac9(1,2)= 1.0;
    frac9(2,2) = 0.75;
    
    frac9(0,3)= 0.50;
    frac9(1,3) = 0.50;
    frac9(2,3) = 0.75;
    
    //  InsertFrac(mesh, corners, matid)
    InsertFrac(gmesh3d,frac1,100);
    InsertFrac(gmesh3d,frac2,100);
//    InsertFrac(gmesh3d,frac3,3);
//    InsertFrac(gmesh3d,frac4,3);
//    InsertFrac(gmesh3d,frac5,3);
//    InsertFrac(gmesh3d,frac6,3);
    InsertFrac(gmesh3d,frac7,100);
//    InsertFrac(gmesh3d,frac8,3);
//    InsertFrac(gmesh3d,frac9,3);
    gmesh->BuildConnectivity();
    return gmesh3d;
    
}

TPZCompMesh * FluxMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    int nbound = sim_data.gamma_ids.size();
   
    
    
    
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedPoisson * volume = new TPZMixedPoisson(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(forcing, 5);
        cmesh->InsertMaterialObject(volume);
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "q_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}
TPZCompMesh * PressureMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);

    for (int ivol=0; ivol < nvols; ivol++) {
        TPZMixedPoisson * volume = new TPZMixedPoisson(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
    }
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "p_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}
void forcing(const TPZVec<REAL> &p, TPZVec<STATE> &f){
    REAL x = p[0];
    REAL y = p[0];
    REAL z = p[0];
    f[0]=0.0*x*y*z;
}

TPZCompMesh * CMeshMixed(TPZGeoMesh * geometry, int order, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    int nbound= sim_data.gamma_ids.size();
    int dim = 3;
    if (nvols<1) {
        std::cout<<"Error: Omega is not defined."<<std::endl;
        DebugStop();
    }
    if (nbound<1) {
        std::cout<<"Error: Gamma is not defined."<<std::endl;
        DebugStop();
    }
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedPoisson * volume = new TPZMixedPoisson(sim_data.omega_ids[ivol],dim);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    
    meshvec.resize(2);
    meshvec[0] = FluxMesh(geometry, order, sim_data);
    meshvec[1] = PressureMesh(geometry, order, sim_data);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvec, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmesh);
    
    std::cout << "Created multi physics mesh\n";
    if (sim_data.IsMHMQ) {
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    else{
        TPZCompMeshTools::GroupElements(cmesh);
        std::cout << "Created grouped elements\n";
        bool keepmatrix = false;
        bool keeponelagrangian = true;
        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
        std::cout << "Created condensed elements\n";
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    
#ifdef PZDEBUG2
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Dual_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    return cmesh;
}

TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, int order, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    int nbound= sim_data.gamma_ids.size();
    if (nvols<1) {
        std::cout<<"Error: Omega is not defined."<<std::endl;
        DebugStop();
    }
    if (nbound<1) {
        std::cout<<"Error: Gamma is not defined."<<std::endl;
        DebugStop();
    }
    
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry);
    
    TPZFNMatrix<9,STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedPoisson * volume = new TPZMixedPoisson(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    
    TPZManVector<TPZCompMesh * ,2> mesh_vec(2);
    mesh_vec[0] = FluxMesh(geometry, order, sim_data);
    mesh_vec[1] = PressureMesh(geometry, order, sim_data);
    TPZManVector<int,5> active_approx_spaces(2); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);

    std::cout << "Created multi physics mesh\n";
    if (sim_data.IsMHMQ) {
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    else{
        TPZCompMeshTools::GroupElements(cmesh);
        std::cout << "Created grouped elements\n";
        bool keepmatrix = false;
        bool keeponelagrangian = true;
        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
        std::cout << "Created condensed elements\n";
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Dual_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    meshvec = mesh_vec;
    return cmesh;
}

TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
    if (sim_data.UsePardisoQ) {
        
        TPZSymetricSpStructMatrix matrix(cmesh);
        //        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    
    if (sim_data.UseFrontalQ) {
        
        TPZParFrontStructMatrix<TPZFrontSym<STATE> > matrix(cmesh);
        matrix.SetDecomposeType(ELDLt);
        matrix.SetNumThreads(sim_data.n_threads);
        
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        analysis->SetSolver(step);
        
        return analysis;
    }
    else{
        
        TPZSkylineStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ECholesky);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
    }
    
    return analysis;
    
}
void FractureTest(){
    int dimel=1;
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.0);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,1);
    nelx[0]=2;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(EQuadrilateral);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    gengrid.SetRefpatternElements(true);
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    TPZManVector<REAL,3> co(3,0.0);
    co[0] = 0.5;
    co[1] = 1.0;
    
    
    TPZGeoNode * node = gmesh->FindNode(co);
    int index_1 = gmesh->NodeIndex(node);
    
    co[0] = 0.5;
    co[1] = 0.0;
    
    node = gmesh->FindNode(co);
    int index_2 = gmesh->NodeIndex(node);
    
    TPZVec<int64_t> cords(2);
    
    cords[0]=index_1;
    cords[1]=index_2;
    int64_t Nels = gmesh->NElements();
    gmesh->CreateGeoElement(EOned, cords, 2, Nels);
    gmesh->BuildConnectivity();
    std::ofstream file("test.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    
    
    //Flux Mesh
    SimulationCase casetest;
    
    casetest.UsePardisoQ=true;
    casetest.IsHybrid=true;
    casetest.omega_ids.push_back(1);
    casetest.permeabilities.push_back(1.0);
    
    
    casetest.gamma_ids.push_back(-1);
    casetest.gamma_ids.push_back(-2);
    casetest.gamma_ids.push_back(-3);
    casetest.gamma_ids.push_back(-4);
    
    //    d = 0;
    //    n = 1;
    casetest.type.push_back(1);
    casetest.type.push_back(0);
    casetest.type.push_back(1);
    casetest.type.push_back(0);
    
    casetest.vals.push_back(0.0);
    casetest.vals.push_back(100.0);
    casetest.vals.push_back(0.0);
    casetest.vals.push_back(10.0);
    
    TPZCompMesh *cmeshq = FluxMesh(gmesh, 1, casetest);
    std::ofstream filecomptxt("cmeshq.txt");
    cmeshq->Print(filecomptxt);
    SeparateConnectsByFracId(cmeshq, 2);
    std::ofstream filecomptxt2("cmeshqq.txt");
    cmeshq->Print(filecomptxt2);
    
    //
    TPZCompMesh *cmeshp = PressureMesh(gmesh, 1, casetest);
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmeshq;
    fmeshvec[1]=cmeshp;
    
    TPZCompMesh *cmixedmesh = CMeshMixed(gmesh, 1, casetest, fmeshvec);
    std::ofstream filemixed("mixedMesh.txt");
    cmixedmesh->Print(filemixed);
    
    
    
    
    TPZCompMesh *cmeshm =NULL;
    if(casetest.IsHybrid){
        //        TPZCompMesh * cmesh_m_Hybrid;
        //        TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
        //        TPZHybridizeHDiv hybridizer;
        //        tie(cmesh_m_Hybrid, meshvector_Hybrid) = hybridizer.Hybridize(cmixedmesh, meshvec, true, -1.);
        //        cmesh_m_Hybrid->InitializeBlock();
        //        cmeshm=cmesh_m_Hybrid;
        cmeshm=cmixedmesh;
    }
    else{
        cmeshm=cmixedmesh;
        
    }
    
    
    
    //
    
    int ncompel = fmeshvec[0]->NElements();
    int ngeoel = gmesh->NElements();
    
    
    //  SeparateConnectsByNeighborhood(cmeshm);
    
    
    
    int index= fmeshvec[0]->Element(1)->ConnectIndex(0);
    std::cout<<"connect index: "<<index<<std::endl;
    //
    
    TPZAnalysis *an = CreateAnalysis(cmeshm, casetest);
    
    std::cout << "Assembly neq = " << cmeshm->NEquations() << std::endl;
    an->Assemble();
    
    std::cout << "Solution of the system" << std::endl;
    an->Solve();
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("Flux");
    scalnames.Push("Pressure");
    scalnames.Push("Permeability");
    
    int div = 0;
    int dim=gmesh->Dimension();
    std::string fileresult("casetest.vtk");
    an->DefineGraphMesh(dim,scalnames,vecnames,fileresult);
    an->PostProcess(div,dim);
    
    
    
    
}
void SeparateConnectsByFracId(TPZCompMesh * cmesh,int fracid){

#ifdef PZDEBUG
    if(!cmesh){
        DebugStop();
    }
#endif
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    gmesh->ResetReference();
    cmesh->LoadReferences();
    cmesh->ComputeNodElCon();
    int64_t nelg = gmesh->NElements();
    int64_t nel = cmesh->NElements();
    for (int64_t iel=0; iel<nelg; iel++) {
        TPZGeoEl *gel = gmesh->Element(iel);
        std::cout<<"mat id: "<<gel->MaterialId()<<std::endl;
        if (gel->MaterialId()==fracid) {
            //            TPZCompEl *cel = gel->Reference();
            //            if (!gel) {continue;}
            TPZGeoElSide gelside = gel->Neighbour(2);
            TPZStack<TPZGeoElSide> allneigh;
            gelside.ComputeNeighbours(allneigh);
            TPZCompEl *cel = allneigh[0].Element()->Reference();
            //desconecta 2d
            int side = allneigh[0].Side()-4;
            TPZConnect &c = cel->Connect(side);
            int64_t cindex = cmesh->AllocateNewConnect(c);
            TPZConnect &newc = cmesh->ConnectVec()[cindex];
            newc = c;
            c.DecrementElConnected();
            newc.DecrementElConnected();
            cel->SetConnectIndex(side, cindex);
            InsertFractureMaterial(cmesh);
            
            
//            //crea nuevo connect side 2
//            TPZCompEl *cel1d = gel->Reference();
//            int64_t cindex2 = cmesh->AllocateNewConnect(c);
//            TPZConnect &newc2 = cmesh->ConnectVec()[cindex];
//            newc2 = c;
//            c.DecrementElConnected();
//            newc2.DecrementElConnected();
//            cel1d->SetConnectIndex(0, cindex2);

            
        }
    }
    
    cmesh->ExpandSolution();
    std::ofstream fileprint("print.txt");
    cmesh->Print(fileprint);
    
}
void InsertFractureMaterial(TPZCompMesh *cmesh){
    //
    REAL perm_0 = 10;
    REAL conv=0;
    TPZVec<REAL> convdir(2,0.0);
    TPZMatPoisson3d *surface = new TPZMatPoisson3d(2,1);
    surface->SetParameters(perm_0, conv, convdir);
    cmesh->InsertMaterialObject(surface);
    cmesh->LoadReferences();
    cmesh->AutoBuild();
    //
}
