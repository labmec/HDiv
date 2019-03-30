
#include "path.h"
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
#include "pzmatmixedpoisson3d.h"
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

#include "TFracture.h"
#include "THybridizeDFN.h"
#include "pzl2projection.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "TPZMixedDarcyFlow.h"


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
        /// check for self-assignment
        if(&copy == this){
            return *this;
        }
        
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
TPZVec<REAL> MidPoint(TPZVec<REAL> & x_i, TPZVec<REAL> & x_e);
void FractureTest();

/// Executes case 1
void Case_1();

/// Executes case 2
void Case_2();

/// Executes cube
void Pretty_cube();

int main(){
    
    std::string source_dir = SOURCE_DIR;
#ifdef LOG4CXX
    std::string log_file = source_dir;
    log_file += "/dfn.cfg";
    InitializePZLOG(log_file);
//    InitializePZLOG();
#endif
    

    
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
    sim.type.push_back(1);
    sim.type.push_back(1);
    sim.type.push_back(1);
    sim.type.push_back(1);
    sim.type.push_back(0);
    sim.type.push_back(0);
    
    sim.vals.push_back(0.0);
    sim.vals.push_back(0.0);
    sim.vals.push_back(0.0);
    sim.vals.push_back(0.0);
    sim.vals.push_back(2.0);
    sim.vals.push_back(1.0);
    
    /// Defining DFN data
    
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id               = 100;
    fracture.m_kappa_normal     = 0.001;
    fracture.m_kappa_tangential = 0.001;
    fracture.m_d_opening        = 1.0e-2;
    fracture_data.push_back(fracture);
    
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
            cmeshm = cmesh_m_Hybrid;
        }else{
            int dimension = 3;
            THybridizeDFN dfn_hybridzer;
            dfn_hybridzer.SetFractureData(fracture_data);
            dfn_hybridzer.SetDimension(dimension);
            
            /// step 1 apply process on dimension 3 entities
            int target_dim = 3;
//            cmeshm = dfn_hybridzer.Hybridize(cmixedmesh,target_dim);
            cmeshm = dfn_hybridzer.Hybridize_II(cmixedmesh,target_dim);
        }

    }
    else{
        cmeshm=cmixedmesh;
        
    }
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmeshm);
    TPZManVector<TPZCompMesh * > mesh_vec = mp_cmesh->MeshVector();

    {
        std::ofstream file_hybrid_mixed_q("Hybrid_mixed_cmesh_q.txt");
        mesh_vec[0]->ComputeNodElCon();
        mesh_vec[0]->Print(file_hybrid_mixed_q);
        
        std::ofstream file_hybrid_mixed_p("Hybrid_mixed_cmesh_p.txt");
        mesh_vec[1]->ComputeNodElCon();
        mesh_vec[1]->Print(file_hybrid_mixed_p);
        
        std::ofstream file_hybrid_mixed("Hybrid_mixed_cmesh.txt");
        cmeshm->ComputeNodElCon();
        cmeshm->Print(file_hybrid_mixed);
    }

    
    TPZAnalysis *an = CreateAnalysis(cmeshm, sim);
    
    std::cout << "Assembly neq = " << cmeshm->NEquations() << std::endl;
    an->Assemble();
    
//    an->Rhs().Print("r = ",std::cout,EMathematicaInput);
    
    std::cout << "Solution of the system" << std::endl;
    an->Solve();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mesh_vec, cmeshm);
    
    std::ofstream file_geo_hybrid("geometry_cube_hybrid.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(cmeshm->Reference(), file_geo_hybrid);
    
    std::ofstream file_geo_hybrid_txt("geometry_cube_hybrid.txt");
    cmeshm->Reference()->Print(file_geo_hybrid_txt);

    std::ofstream file_hybrid_mixed_q("Hybrid_mixed_cmesh_q.txt");
    mesh_vec[0]->Print(file_hybrid_mixed_q);
    
    std::ofstream file_hybrid_mixed_p("Hybrid_mixed_cmesh_p.txt");
    mesh_vec[1]->Print(file_hybrid_mixed_p);
    
    std::ofstream file_hybrid_mixed("Hybrid_mixed_cmesh.txt");
    cmeshm->Print(file_hybrid_mixed);
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    vecnames.Push("kappa");
    scalnames.Push("p");
    
    int div = 0;
    std::string file_reservoir("cube.vtk");
    an->DefineGraphMesh(3,scalnames,vecnames,file_reservoir);
    an->PostProcess(div,3);
    
    { /// fracture postprocessor
        TPZStack<std::string,10> scalnames, vecnames;
        scalnames.Push("state");
        std::string file_frac("fracture.vtk");
        auto material = mesh_vec[1]->FindMaterial(100);
        TPZMixedDarcyFlow * fract_2d = dynamic_cast<TPZMixedDarcyFlow *>(material);
        fract_2d->SetDimension(2);
        TPZAnalysis frac_an(mesh_vec[1],false);
        frac_an.DefineGraphMesh(2,scalnames,vecnames,file_frac);
        frac_an.PostProcess(div,2);
    }

    { /// lagrange postprocessor
        TPZStack<std::string,10> scalnames, vecnames;
        scalnames.Push("state");
        std::string file_frac("lagrange_1d.vtk");
        auto material = mesh_vec[1]->FindMaterial(100);
        TPZMixedDarcyFlow * fract_2d = dynamic_cast<TPZMixedDarcyFlow *>(material);
        fract_2d->SetDimension(1);
        TPZAnalysis frac_an(mesh_vec[1],false);
        frac_an.DefineGraphMesh(1,scalnames,vecnames,file_frac);
        frac_an.PostProcess(div,1);
    }
    
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
    TPZManVector<REAL,3> co1(3,0.0);
    co1[0] = corners(0,0);
    co1[1] = corners(1,0);
    co1[2] = corners(2,0);
    
    TPZGeoNode * node = gmesh->FindNode(co1);
    int index_1 = gmesh->NodeIndex(node);
    
    TPZManVector<REAL,3> co2(3,0.0);
    co2[0] = corners(0,1);
    co2[1] = corners(1,1);
    co2[2] = corners(2,1);
    
    node = gmesh->FindNode(co2);
    int index_2 = gmesh->NodeIndex(node);
    
    TPZManVector<REAL,3> co3(3,0.0);
    co3[0] = corners(0,2);
    co3[1] = corners(1,2);
    co3[2] = corners(2,2);
    
    node = gmesh->FindNode(co3);
    int index_3 = gmesh->NodeIndex(node);
    
    TPZManVector<REAL,3> co4(3,0.0);
    co4[0] = corners(0,3);
    co4[1] = corners(1,3);
    co4[2] = corners(2,3);
    
    node = gmesh->FindNode(co4);
    int index_4 = gmesh->NodeIndex(node);
    
    int64_t n_nodes_base = gmesh->NNodes();
    int64_t n_nodes = n_nodes_base + 5;
    gmesh -> NodeVec().Resize(n_nodes);
    gmesh->SetMaxNodeId(n_nodes-1);
    
    TPZManVector<REAL,3> co5,co6,co7,co8,co9;
    co5 = MidPoint(co1, co2);
    co6 = MidPoint(co2, co3);
    co7 = MidPoint(co3, co4);
    co8 = MidPoint(co4, co1);
    co9 = MidPoint(co8, co6);
    
    int64_t index_5 = n_nodes_base;
    TPZGeoNode(index_5, co5, *gmesh);

    int64_t index_6 = n_nodes_base + 1;
    TPZGeoNode(index_6, co6, *gmesh);
    
    int64_t index_7 = n_nodes_base + 2;
    TPZGeoNode(index_7, co7, *gmesh);
    
    int64_t index_8 = n_nodes_base + 3;
    TPZGeoNode(index_8, co8, *gmesh);
    
    int64_t index_9 = n_nodes_base + 4;
    TPZGeoNode(index_9, co9, *gmesh);
    

    TPZVec<int64_t> node_index(4);
    int64_t el_index;;
    
    node_index[0]=index_1;
    node_index[1]=index_5;
    node_index[2]=index_9;
    node_index[3]=index_8;
    el_index = gmesh->NElements();
    gmesh->CreateGeoElement(EQuadrilateral, node_index, matid, el_index);
    
    node_index[0]=index_5;
    node_index[1]=index_2;
    node_index[2]=index_6;
    node_index[3]=index_9;
    el_index = gmesh->NElements();
    gmesh->CreateGeoElement(EQuadrilateral, node_index, matid, el_index);

    node_index[0]=index_8;
    node_index[1]=index_9;
    node_index[2]=index_7;
    node_index[3]=index_4;
    el_index = gmesh->NElements();
    gmesh->CreateGeoElement(EQuadrilateral, node_index, matid, el_index);

    node_index[0]=index_9;
    node_index[1]=index_6;
    node_index[2]=index_3;
    node_index[3]=index_7;
    el_index = gmesh->NElements();
    gmesh->CreateGeoElement(EQuadrilateral, node_index, matid, el_index);
    
}

TPZVec<REAL> MidPoint(TPZVec<REAL> & x_i, TPZVec<REAL> & x_e) {
    
    if (x_i.size()!=x_e.size()) {
        DebugStop();
    }
    TPZVec<REAL> x(x_i.size());
    for (int i = 0;  i < x_i.size(); i++) {
        x[i] = x_i[i] + 0.5*(x_e[i] - x_i[i]);
    }
    
    return x;
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
    gmesh->BuildConnectivity();
    int n_layers = 2;
    TPZExtendGridDimension extend(gmesh,1.0/n_layers);
    extend.SetElType(1);
    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(n_layers,-5,-6);
    gmesh3d->BuildConnectivity();
    
    int fracture_id= 100;
    int dim = gmesh3d->Dimension();
    /// Insert fractures
    {
        std::map<std::pair<int,int>,int> vol_to_vol_side_indexes;
        std::vector<int> sides = {20,21,22,23,24,25};
        for (auto gel : gmesh3d->ElementVec()) {
            
            if (!gel) continue;
            if (gel->Dimension() != dim) continue;
            
            for (auto side: sides) {
                TPZStack<TPZGeoElSide> all_neigh;
                TPZGeoElSide gelside(gel, side);
                gelside.AllNeighbours(all_neigh);
                std::set<int> vols;
                vols.insert(gel->Index());
                for (auto gel_side : all_neigh) {
                    if (gel_side.Element()->Dimension() == dim) {
                        int vol_index = gel_side.Element()->Index();
                        vols.insert(vol_index);
                    }
                }
                if (vols.size()!=2) {
                    continue;
                }
                
                std::set<int>::iterator it;
                it=vols.begin();
                int left = *it;
                ++it;
                int right = *it;
                
                vol_to_vol_side_indexes.insert(std::make_pair(std::make_pair(left, right),side));
                
            }
        }
        
        for (auto chunk : vol_to_vol_side_indexes) {
            TPZGeoEl * gel_l = gmesh3d->Element(chunk.first.first);
            TPZGeoEl * gel_r = gmesh3d->Element(chunk.first.second);
            TPZGeoElSide gelside(gel_l, chunk.second);
            TPZManVector<REAL,2> qsi(2,0.0);
            TPZManVector<REAL,3> normal(3);
            gelside.Normal(qsi, gel_l, gel_r, normal);
            if (fabs(fabs(normal[2]) - 1.0) <= 1.0e-8) {
                TPZGeoElBC gbc(gelside, fracture_id);
            }
            if (fabs(fabs(normal[1]) - 1.0) <= 1.0e-8) {
                TPZGeoElBC gbc(gelside, fracture_id);
            }
            if (fabs(fabs(normal[0]) - 1.0) <= 1.0e-8) {
                TPZGeoElBC gbc(gelside, fracture_id);
            }
        }
        
    }
    gmesh3d->BuildConnectivity();
    
    bool insert_fractrures_intersection_Q = true;
    if(insert_fractrures_intersection_Q){
        /// Insert fractures intersections
        {
            std::map<std::pair<int,int>,int> vol_to_vol_side_indexes;
            std::vector<int> sides = {4,5,6,7};
            for (auto gel : gmesh3d->ElementVec()) {
                
                if (!gel) continue;
                if (gel->Dimension() != dim-1 || gel->MaterialId() != 100) continue;
                
                for (auto side: sides) {
                    TPZStack<TPZGeoElSide> all_neigh;
                    TPZGeoElSide gelside(gel, side);
                    gelside.AllNeighbours(all_neigh);
                    std::set<int> vols;
                    vols.insert(gel->Index());
                    for (auto gel_side : all_neigh) {
                        if (gel_side.Element()->Dimension() == dim-1 && gel_side.Element()->MaterialId() == 100) {
                            int vol_index = gel_side.Element()->Index();
                            vols.insert(vol_index);
                        }
                    }
                    if (vols.size()!=4) {
                        continue;
                    }
                    
                    std::set<int>::iterator it;
                    it=vols.begin();
                    int left = *it;
                    ++it;
                    int right = *it;
                    
                    vol_to_vol_side_indexes.insert(std::make_pair(std::make_pair(left, right),side));
                    
                }
            }
            
            for (auto chunk : vol_to_vol_side_indexes) {
                TPZGeoEl * gel_l = gmesh3d->Element(chunk.first.first);
                TPZGeoElSide gelside(gel_l, chunk.second);
                TPZGeoElBC gbc(gelside, fracture_id);
            }
            
        }
    }
    
    gmesh3d->BuildConnectivity();
    return gmesh3d;
    
    
//    TPZFMatrix<REAL> frac1(3,4,0.0);
//
//    frac1(0,0)=0.5;
//    frac1(1,0)=0.0;
//    frac1(2,0)=0.0;
//
//    frac1(0,1)=0.5;
//    frac1(1,1)=0.0;
//    frac1(2,1)=1.0;
//
//    frac1(0,2)=0.5;
//    frac1(1,2)=1.0;
//    frac1(2,2)=1.0;
//
//    frac1(0,3)=0.5;
//    frac1(1,3)=1.0;
//    frac1(2,3)=0.0;
//
//    TPZFMatrix<REAL> frac2(3,4,0.0);
//
//    frac2(0,0) = 1.0;
//    frac2(1,0) = 0.5;
//    frac2(2,0) = 0.0;
//
//    frac2(0,1) = 1.0;
//    frac2(1,1) = 0.5;
//    frac2(2,1) = 1.0;
//
//    frac2(0,2) = 0.0;
//    frac2(1,2) = 0.5;
//    frac2(2,2) = 1.0;
//
//    frac2(0,3) = 0.0;
//    frac2(1,3) = 0.5;
//    frac2(2,3) = 0.0;
//
//    //frac3
//    TPZFMatrix<REAL> frac3(3,4,0.0);
//    frac3(0,0) = 0.75;
//    frac3(1,0) = 0.5;
//    frac3(2,0) = 0.5;
//
//    frac3(0,1) = 0.75;
//    frac3(1,1) = 0.5;
//    frac3(2,1) = 1.0;
//
//    frac3(0,2) = 0.75;
//    frac3(1,2) = 1.0;
//    frac3(2,2) = 1.0;
//
//    frac3(0,3) = 0.75;
//    frac3(1,3) = 1.0;
//    frac3(2,3) = 0.5;
//
//
//    //Set Frac_4
//    TPZFMatrix<REAL> frac4(3,4,0.0);
//    frac4(0,0) = 1.0;
//    frac4(1,0) = 0.75;
//    frac4(2,0) = 0.5;
//
//    frac4(0,1) = 1.0;
//    frac4(1,1) = 0.75;
//    frac4(2,1) = 1.0;
//
//    frac4(0,2) = 0.5;
//    frac4(1,2)= 0.75;
//    frac4(2,2) = 1.0;
//
//    frac4(0,3)= 0.50;
//    frac4(1,3) = 0.75;
//    frac4(2,3) = 0.50;
//
//
//    //Set Frac_5
//    TPZFMatrix<REAL> frac5(3,4,0.0);
//    frac5(0,0) = 0.675;
//    frac5(1,0) = 0.5;
//    frac5(2,0) = 0.5;
//
//    frac5(0,1) = 0.675;
//    frac5(1,1) = 0.5;
//    frac5(2,1) = 0.75;
//
//    frac5(0,2) = 0.675;
//    frac5(1,2) = 0.75;
//    frac5(2,2) = 0.75;
//
//    frac5(0,3) = 0.675;
//    frac5(1,3) = 0.75;
//    frac5(2,3) = 0.5;
//
//
//    //Set Frac_6
//    TPZFMatrix<REAL> frac6(3,4,0.0);
//    frac6(0,0) = 0.75;
//    frac6(1,0) = 0.675;
//    frac6(2,0) = 0.5;
//
//    frac6(0,1) = 0.75;
//    frac6(1,1) = 0.675;
//    frac6(2,1) = 0.75;
//
//    frac6(0,2) = 0.5;
//    frac6(1,2) = 0.675;
//    frac6(2,2) = 0.75;
//
//    frac6(0,3) = 0.50;
//    frac6(1,3) = 0.675;
//    frac6(2,3)= 0.50;
//
//    TPZFMatrix<REAL> frac7(3,4,0.0);
//
//    frac7(0,0) = 0.0;
//    frac7(1,0) = 0.0;
//    frac7(2,0) = 0.5;
//
//    frac7(0,1) = 1.0;
//    frac7(1,1) = 0.0;
//    frac7(2,1) = 0.5;
//
//    frac7(0,2) = 1.0;
//    frac7(1,2) = 1.0;
//    frac7(2,2) = 0.5;
//
//    frac7(0,3) = 0.0;
//    frac7(1,3) = 1.0;
//    frac7(2,3)= 0.5;
//
//
//    TPZFMatrix<REAL> frac8(3,4,0.0);
//    frac8(0,0) = 0.75;
//    frac8(1,0) = 0.5;
//    frac8(2,0) = 0.675;
//
//    frac8(0,1) = 0.75;
//    frac8(1,1) = 0.75;
//    frac8(2,1) = 0.675;
//
//    frac8(0,2) = 0.5;
//    frac8(1,2) = 0.75;
//    frac8(2,2) = 0.675;
//
//    frac8(0,3) = 0.50;
//    frac8(1,3) = 0.5;
//    frac8(2,3)= 0.675;
//
//    //Set Frac_4
//    TPZFMatrix<REAL> frac9(3,4,0.0);
//    frac9(0,0) = 1.0;
//    frac9(1,0) = 0.5;
//    frac9(2,0) = 0.75;
//
//    frac9(0,1) = 1.0;
//    frac9(1,1) = 1.0;
//    frac9(2,1) = 0.75;
//
//    frac9(0,2) = 0.5;
//    frac9(1,2)= 1.0;
//    frac9(2,2) = 0.75;
//
//    frac9(0,3)= 0.50;
//    frac9(1,3) = 0.50;
//    frac9(2,3) = 0.75;
//
//    //  InsertFrac(mesh, corners, matid)
//    InsertFrac(gmesh3d,frac1,100);
////    InsertFrac(gmesh3d,frac2,100);
////    InsertFrac(gmesh3d,frac3,3);
////    InsertFrac(gmesh3d,frac4,3);
////    InsertFrac(gmesh3d,frac5,3);
////    InsertFrac(gmesh3d,frac6,3);
////    InsertFrac(gmesh3d,frac7,100);
////    InsertFrac(gmesh3d,frac8,3);
////    InsertFrac(gmesh3d,frac9,3);
//    gmesh3d->BuildConnectivity();
//    return gmesh3d;
    
}

TPZCompMesh * FluxMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    int nbound = sim_data.gamma_ids.size();
   

    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMatMixedPoisson3D * volume = new TPZMatMixedPoisson3D(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(forcing, 5);
        cmesh->InsertMaterialObject(volume);
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
                TPZMaterial * face_1d = volume->CreateBC(volume,1000*sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face_1d);
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
        TPZMatMixedPoisson3D * volume = new TPZMatMixedPoisson3D(sim_data.omega_ids[ivol],dimension);
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
        TPZMatMixedPoisson3D * volume = new TPZMatMixedPoisson3D(sim_data.omega_ids[ivol],dim);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
                TPZMaterial * face_1d = volume->CreateBC(volume,1000*sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face_1d);
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

        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
                TPZMaterial * face_1d = volume->CreateBC(volume,1000*sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face_1d);
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
//        TPZCompMeshTools::GroupElements(cmesh);
//        std::cout << "Created grouped elements\n";
//        bool keepmatrix = false;
//        bool keeponelagrangian = true;
//        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
//        std::cout << "Created condensed elements\n";
//        cmesh->CleanUpUnconnectedNodes();
//        cmesh->ExpandSolution();
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "Dual_cmesh" << ".txt";
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
