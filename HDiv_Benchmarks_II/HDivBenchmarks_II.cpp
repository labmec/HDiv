
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
#include "CreateDFN.h"

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
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSpStructMatrix.h"

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
#include "TPZTracerFlow.h"
#include "TPZNormalDarcyFlow.h"

#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

using namespace std;

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
    TPZStack<int>   omega_dim;
    TPZStack<int>   gamma_ids;
    TPZStack<int>   gamma_dim;
    TPZStack<REAL>   permeabilities;
    TPZStack<REAL>   porosities;
    TPZStack<REAL>   type;
    TPZStack<REAL>   vals;
    REAL            c_inlet;

    SimulationCase() : IsMHMQ(false), UsePardisoQ(true), IsHybrid(false),UseFrontalQ(false), UseGmshMeshQ(false), NonAffineQ(false), elemen_type(0), n_h_levels(0), n_p_levels(1), n_acc_terms(0), int_order(1), n_threads(0),perturbation_type(0), mesh_type(""), domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),omega_dim(),gamma_ids(), gamma_dim(), permeabilities(),porosities(), type(), vals(), c_inlet(0)
    {

    }

    SimulationCase(const SimulationCase &copy) :IsHybrid(copy.IsHybrid) ,IsMHMQ(copy.IsMHMQ), UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
    UseGmshMeshQ(copy.UseGmshMeshQ), NonAffineQ(copy.NonAffineQ), elemen_type(copy.elemen_type), n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), n_acc_terms(copy.n_acc_terms), int_order(copy.int_order),n_threads(copy.n_threads),perturbation_type(copy.perturbation_type), mesh_type(copy.mesh_type), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
    dump_folder(copy.dump_folder), omega_ids(copy.omega_ids), omega_dim(copy.omega_dim), gamma_ids(copy.gamma_ids), gamma_dim(copy.gamma_dim),
    permeabilities(copy.permeabilities),porosities(copy.porosities),
    type(copy.type),
    vals(copy.vals),c_inlet(copy.c_inlet)
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
        omega_dim = copy.omega_dim;
        gamma_ids = copy.gamma_ids;
        gamma_dim = copy.gamma_dim;
        permeabilities=copy.permeabilities;
        porosities = copy.porosities;
        type=copy.type;
        vals=copy.vals;
        c_inlet=copy.c_inlet;
        return *this;
    }
};
TPZCompMesh * FluxMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data);
TPZCompMesh * PressureMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data);
TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, int order, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);

/// Executes cube
void Pretty_cube();
void Pretty_cubeII();

int main(){
    
    std::string source_dir = SOURCE_DIR;
#ifdef LOG4CXX
    std::string log_file = source_dir;
    log_file += "/dfn.cfg";
    InitializePZLOG(log_file);
#endif
    
 Pretty_cubeII();
//  Pretty_cube();
   
    
}
void Pretty_cube(){
    
    int h_level = 0;
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 8;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// not used but inserted
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// C inlet value
    sim.c_inlet = 1.0;
    
    int bc_inlet  = 3;
    int bc_outlet = 4;
    int bc_non_flux = 5;
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(3);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 2.0;
    REAL p_outlet = 1.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    sim.vals.push_back(qn);
    
    /// Defining DFN boundary data (id,bc_type,data)
    std::vector<std::tuple<int,int,REAL>> bc_ids_2d;
    bc_ids_2d.push_back(std::make_tuple(bc_inlet,bc_type_D,p_inlet));
    bc_ids_2d.push_back(std::make_tuple(bc_outlet,bc_type_D,p_outlet));
    bc_ids_2d.push_back(std::make_tuple(bc_non_flux,bc_type_N,qn));
    
    int bc_1d_inlet  = 130;
    int bc_1d_outlet = 140;
    int bc_1d_non_flux = 150;
    
    int bc_0d_inlet  = 230;
    int bc_0d_outlet = 240;
    int bc_0d_non_flux = 250;
    
    std::map<int,int> bc_ids_1d_map;
    bc_ids_1d_map.insert(std::make_pair(bc_inlet,bc_1d_inlet));
    bc_ids_1d_map.insert(std::make_pair(bc_outlet,bc_1d_outlet));
    bc_ids_1d_map.insert(std::make_pair(bc_non_flux,bc_1d_non_flux));
    
    std::map<int,int> bc_ids_0d_map;
    bc_ids_0d_map.insert(std::make_pair(bc_inlet,bc_0d_inlet));
    bc_ids_0d_map.insert(std::make_pair(bc_outlet,bc_0d_outlet));
    bc_ids_0d_map.insert(std::make_pair(bc_non_flux,bc_0d_non_flux));
    
    /// Defining DFN data
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id.insert(6);
    fracture.m_id.insert(7);
    fracture.m_id.insert(8);
    fracture.m_dim              = 2;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 0.5;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(9);
    fracture.m_dim              = 1;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(10);
    fracture.m_dim              = 0;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    
    
    /// Benchmarks Material ID convention
    /// 1 and 2 for 3D matrix
    /// 3,4, and 5 for 2D matrix boundaries (3 -> inlet, 4 -> outlet, 5 -> impervious)
    /// 6 fractures
    /// 7 fractures intersections
    /// 8 crossing intersections
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[3]["RockMatrix_1"] = 1;
    dim_name_and_physical_tag[3]["RockMatrix_2"] = 2;
    dim_name_and_physical_tag[2]["BCInlet"] = 3;
    dim_name_and_physical_tag[2]["BCOutlet"] = 4;
    dim_name_and_physical_tag[2]["BCImpervious"] = 5;
    dim_name_and_physical_tag[2]["Fractures_1"] = 6;
    dim_name_and_physical_tag[2]["Fractures_2"] = 7;
    dim_name_and_physical_tag[2]["Fractures_3"] = 8;
    dim_name_and_physical_tag[1]["FracturesIntersections"] = 9;
    dim_name_and_physical_tag[0]["CrossingIntresections"] = 10;
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
    std::string file_gmsh = source_dir + "/meshes/the_cuttest_cube/cube.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
    
#ifdef PZDEBUG
    std::ofstream file("geometry_cube_base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    std::ofstream file_txt("geometry_cube_base.txt");
    gmesh->Print(file_txt);
#endif
    
    int p_order = 1;
    TPZVec<TPZCompMesh *> meshvec;
    TPZCompMesh *cmixedmesh = NULL;
    cmixedmesh = MPCMeshMixed(gmesh, p_order, sim, meshvec);
#ifdef PZDEBUG
    std::ofstream filemixed("mixed_cmesh.txt");
    cmixedmesh->Print(filemixed);
#endif
    
    TPZCompMesh *cmeshm =NULL;
    THybridizeDFN dfn_hybridzer;
    dfn_hybridzer.SetFractureData(fracture_data);
    
    dfn_hybridzer.SetReservoirBoundaryData(bc_ids_2d);
    dfn_hybridzer.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
    dfn_hybridzer.SetMapReservoirBCToDFNBC0DIds(bc_ids_0d_map);
    cmeshm = dfn_hybridzer.Hybridize(cmixedmesh);

    
    
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmeshm);
    TPZManVector<TPZCompMesh *> vectormesh = mp_cmesh->MeshVector();
    TPZCompMesh *Flux_Mesh = vectormesh[0] ;
    TPZCompMesh *Pressure_Mesh = vectormesh[1] ;
    std::ofstream file3("PRESSURE.vtk");
     std::ofstream file2("FLUX.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(Flux_Mesh, file2);
    TPZVTKGeoMesh::PrintCMeshVTK(Pressure_Mesh, file3);
    
#ifdef PZDEBUG
    {
        std::ofstream file("geometry_cube_after.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
        std::ofstream file_txt("geometry_cube_after.txt");
        gmesh->Print(file_txt);
    }
#endif
    
    
   // TPZCompMesh *s_cmesh = CreateTransportMesh(mp_cmesh, fracture_data);
    
#ifdef PZDEBUG
    {
//        std::ofstream file("transport_mesh.vtk");
//        TPZVTKGeoMesh::PrintCMeshVTK(s_cmesh, file);
//        std::ofstream out("transport_mesh.txt");
//        s_cmesh->Print(out);
    }
#endif
    
//    TPZManVector<TPZCompMesh *,3> meshtrvec(3);
//    meshtrvec[0] = meshvec[0];
//    meshtrvec[1] = meshvec[1];
//    meshtrvec[2] = s_cmesh;
    
   // TPZMultiphysicsCompMesh *cmesh_transport = MPTransportMesh(mp_cmesh, fracture_data, sim, meshtrvec);
  //  TPZAnalysis * tracer_analysis = CreateTransportAnalysis(cmesh_transport, sim);
    
    bool solve_dfn_problem_Q = true;
    if (solve_dfn_problem_Q) {
        TPZManVector<TPZCompMesh * > mesh_vec = mp_cmesh->MeshVector();
        
#ifdef PZDEBUG
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
#endif
        
        std::cout << "Condensing DFN equations." << std::endl;
        std::cout << "DFN neq before condensation = " << mp_cmesh->NEquations() << std::endl;
        dfn_hybridzer.GroupElements(mp_cmesh);
        std::cout << "DFN neq = " << mp_cmesh->NEquations() << std::endl;
        
#ifdef PZDEBUG
        {
            
            std::ofstream file_hybrid_condensed("Hybrid_mixed_condensed.txt");
            mp_cmesh->ComputeNodElCon();
            mp_cmesh->Print(file_hybrid_condensed);
        }
#endif
        
        std::cout << "DFN equations are condensed." << std::endl;
        
        TPZAnalysis *an = CreateAnalysis(mp_cmesh, sim);
        std::cout << "Assembly DFN problem neq = " << mp_cmesh->NEquations() << std::endl;
        an->Assemble();
        std::cout << "Assembly for DFN complete." << std::endl;
        
        std::cout << "Solving DFN problem." << std::endl;
        an->Solve();
        std::cout << "DFN problem solved." << std::endl;
        mp_cmesh->LoadSolutionFromMultiPhysics();
        
#ifdef PZDEBUG
        std::ofstream file_geo_hybrid("geometry_cube_hybrid.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(cmeshm->Reference(), file_geo_hybrid);
#endif
        
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        scalnames.Push("p");
        
        int div = 0;
        std::set<int> mat_id_3D;
        mat_id_3D.insert(1);
        mat_id_3D.insert(2);
        std::string file_reservoir("cube.vtk");
        an->DefineGraphMesh(3,mat_id_3D,scalnames,vecnames,file_reservoir);
        an->PostProcess(div,3);
        
        std::set<int> mat_id_2D;
        mat_id_2D.insert(6);
        mat_id_2D.insert(7);
        mat_id_2D.insert(8);
        std::string file_frac("fracture.vtk");
        an->DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
        an->PostProcess(div,2);
        
        std::set<int> mat_id_1D;
        mat_id_1D.insert(7);
         mat_id_1D.insert(8);
         mat_id_1D.insert(9);
        std::string file_frac_intersections("fracture_intersections.vtk");
        an->DefineGraphMesh(1,mat_id_1D,scalnames,vecnames,file_frac_intersections);
        an->PostProcess(div,1);
    }
    
    return;
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
            }
        }
    }
    
    cmesh->SetDimModel(dimension);
    TPZManVector<TPZCompMesh * ,2> mesh_vec(2);
    /// generate the 3D flux mesh (only)
    mesh_vec[0] = FluxMesh(geometry, order, sim_data);
    mesh_vec[1] = PressureMesh(geometry, order, sim_data);
    TPZManVector<int,5> active_approx_spaces(2); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,mesh_vec);
    
    
    
    if (sim_data.IsMHMQ) {
        cmesh->CleanUpUnconnectedNodes();
        cmesh->ExpandSolution();
    }
    //    else{
    //        TPZCompMeshTools::GroupElements(cmesh);
    //        std::cout << "Created grouped elements\n";
    //        bool keepmatrix = false;
    //        bool keeponelagrangian = true;
    //        TPZCompMeshTools::CreatedCondensedElements(cmesh, keeponelagrangian, keepmatrix);
    //        std::cout << "Created condensed elements\n";
    //        cmesh->CleanUpUnconnectedNodes();
    //        cmesh->ExpandSolution();
    //    }
    
    std::cout << "Created multi-physics DFN mesh\n";
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name  << "Dual_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    meshvec = mesh_vec;
    return cmesh;
}
TPZCompMesh * FluxMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    int dimension = geometry->Dimension();
    int nvols = sim_data.omega_ids.size();
    int nbound = sim_data.gamma_ids.size();
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(sim_data.omega_ids[ivol],sim_data.omega_dim[ivol]);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        
        cmesh->InsertMaterialObject(volume);
        if (ivol==0) {
            for (int ibound=0; ibound<nbound; ibound++) {
                if(sim_data.gamma_dim[ibound] != sim_data.omega_dim[ivol]) continue;
                val2(0,0)=sim_data.vals[ibound];
                int condType=sim_data.type[ibound];
                TPZMaterial * face = volume->CreateBC(volume,sim_data.gamma_ids[ibound],condType,val1,val2);
                cmesh->InsertMaterialObject(face);
            }
        }
    }
    
//    for (int dim=3; dim>0; dim--) {
//        cmesh->Reference()->ResetReference();
//        cmesh->SetDimModel(dim);
//        cmesh->SetDefaultOrder(order);
//        std::set<int> matids;
//        for (int ivol=0; ivol<nvols; ivol++) {
//            if(sim_data.omega_dim[ivol]==dim) matids.insert(sim_data.omega_ids[ivol]);
//        }
//        for (int ibound=0; ibound<nbound; ibound++) {
//            if(sim_data.gamma_dim[ibound]==dim) matids.insert(sim_data.gamma_ids[ibound]);
//        }
//        cmesh->SetAllCreateFunctionsHDiv();
//        cmesh->AutoBuild(matids);
//    }
    
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
    cmesh->SetDimModel(3);
    cmesh->AutoBuild();
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "q_cmesh_raw" << ".txt";
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
    
    std::set<int> matids;
    for (int ivol=0; ivol < nvols; ivol++) {
        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        cmesh->InsertMaterialObject(volume);
        if (sim_data.omega_dim[ivol] == dimension) {
            matids.insert(sim_data.omega_ids[ivol]);
        }
    }
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    
    cmesh->AutoBuild(matids);
    cmesh->InitializeBlock();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name << "p_cmesh_raw" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
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
void Pretty_cubeII(){
    
    //Lectura de la malha
    TPZManVector<std::map<std::string,int>,4> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[3]["RockMatrix_1"] = 1;
    dim_name_and_physical_tag[3]["RockMatrix_2"] = 2;
    dim_name_and_physical_tag[2]["BCInlet"] = 3;
    dim_name_and_physical_tag[2]["BCOutlet"] = 4;
    dim_name_and_physical_tag[2]["BCImpervious"] = 5;
    dim_name_and_physical_tag[2]["Fractures_1"] = 26;
    dim_name_and_physical_tag[2]["Fractures_2"] = 26;
    dim_name_and_physical_tag[2]["Fractures_3"] = 26;
    dim_name_and_physical_tag[1]["FracturesIntersections"] = 27;
    dim_name_and_physical_tag[0]["CrossingIntresections"] = 28;
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
    std::string file_gmsh = source_dir + "/meshes/the_cuttest_cube/cube.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
    
#ifdef PZDEBUG
    std::ofstream file("geometry_cube_base.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    std::ofstream file_txt("geometry_cube_base.txt");
    gmesh->Print(file_txt);
    gmesh->SetDimension(3);
#endif
    
    //Configuracion de Materiales
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 1;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// not used but inserted
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(1.0);
    
    /// C inlet value
    sim.c_inlet = 1.0;
    
    int bc_inlet  = 3;
    int bc_outlet = 4;
    int bc_non_flux = 5;
    
    sim.gamma_ids.push_back(bc_inlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_outlet);
    sim.gamma_dim.push_back(3);
    sim.gamma_ids.push_back(bc_non_flux);
    sim.gamma_dim.push_back(3);
    
    int bc_type_D = 0;    //    D = 0;
    int bc_type_N = 1;    //    N = 1;
    REAL p_inlet  = 2.0;
    REAL p_outlet = 1.0;
    REAL qn       = 0.0;
    
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_D);
    sim.type.push_back(bc_type_N);
    
    sim.vals.push_back(p_inlet);
    sim.vals.push_back(p_outlet);
    sim.vals.push_back(qn);
    TPZVec<TPZCompMesh *> meshvec;
    
    
    
    
    /// Defining DFN data
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id.insert(26);
    fracture.m_dim              = 2;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 0.5;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(27);
    fracture.m_dim              = 1;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
    fracture.m_id.clear();
    fracture.m_id.insert(28);
    fracture.m_dim              = 0;
    fracture.m_kappa_normal     = 1.0e20;
    fracture.m_kappa_tangential = 1.0;
    fracture.m_d_opening        = 1.0;
    fracture.m_porosity         = 1.0;
    fracture_data.push_back(fracture);
 
    
    
    //
    /// Defining DFN boundary data (id,bc_type,data)
    std::vector<std::tuple<int,int,REAL>> bc_ids_2d;
    bc_ids_2d.push_back(std::make_tuple(bc_inlet,bc_type_D,p_inlet));
    bc_ids_2d.push_back(std::make_tuple(bc_outlet,bc_type_D,p_outlet));
    bc_ids_2d.push_back(std::make_tuple(bc_non_flux,bc_type_N,qn));
    
    int bc_1d_inlet  = 130;
    int bc_1d_outlet = 140;
    int bc_1d_non_flux = 150;
    
    int bc_0d_inlet  = 230;
    int bc_0d_outlet = 240;
    int bc_0d_non_flux = 250;
    
    std::map<int,int> bc_ids_1d_map;
    bc_ids_1d_map.insert(std::make_pair(bc_inlet,bc_1d_inlet));
    bc_ids_1d_map.insert(std::make_pair(bc_outlet,bc_1d_outlet));
    bc_ids_1d_map.insert(std::make_pair(bc_non_flux,bc_1d_non_flux));
    
    std::map<int,int> bc_ids_0d_map;
    bc_ids_0d_map.insert(std::make_pair(bc_inlet,bc_0d_inlet));
    bc_ids_0d_map.insert(std::make_pair(bc_outlet,bc_0d_outlet));
    bc_ids_0d_map.insert(std::make_pair(bc_non_flux,bc_0d_non_flux));
    //
    CreateDFN dfn;
    int flux_trace_id=10;
    int lagrange_id=11;
    int flux_resistivity_id=13;
    int mp_nterface_id =12;
    
    dfn.SetFractureData(fracture_data);
    dfn.SetReservoirBoundaryData(bc_ids_2d);
    dfn.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
    dfn.SetMapReservoirBCToDFNBC0DIds(bc_ids_0d_map);
    
    
    TPZMultiphysicsCompMesh *mp_initial_cmesh = MPCMeshMixed(gmesh, 1, sim, meshvec);
    TPZCompMesh * cmesh_m_Hybrid;
    TPZManVector<TPZCompMesh*, 3> meshvector_Hybrid(3);
    TPZHybridizeHDiv hybridizer;
    hybridizer.SetPeriferalMaterialIds(10, 11, 12);
    TPZManVector<TPZCompMesh *> mp_mesh_vector = mp_initial_cmesh->MeshVector();
    TPZCompMesh * q_mesh = mp_mesh_vector[0];
    TPZCompMesh *p_mesh = mp_mesh_vector[1];
    cmesh_m_Hybrid = hybridizer.Hybridize(mp_initial_cmesh, false, 1.);

    dfn.BuildMixedOperatorOnFractures(1, 2, mp_initial_cmesh, flux_trace_id, lagrange_id, flux_resistivity_id, mp_nterface_id);
    dfn.LoadReferencesByDimension(q_mesh, 2);
    
    std::ofstream fileq("flujo1.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(q_mesh, fileq);
    q_mesh->Reference()->SetDimension(2);
    hybridizer.Hybridize(mp_initial_cmesh, false,1);
    
    dfn.BuildMixedOperatorOnFractures(1, 1, mp_initial_cmesh, flux_trace_id, lagrange_id, flux_resistivity_id, mp_nterface_id);
    dfn.LoadReferencesByDimension(q_mesh, 1);
    
    std::ofstream fileq2("flujo2.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(q_mesh, fileq2);
    
   
//    TPZStack<std::pair<int, int> > gel_index_and_order_lagrange_mult ;
//    dfn.ClassifyCompelSides(1, q_mesh, gel_index_and_order_lagrange_mult, flux_trace_id, lagrange_id, flux_resistivity_id);
q_mesh->Reference()->SetDimension(1);
    hybridizer.Hybridize(mp_initial_cmesh, false,1);
  
   
    
    std::ofstream fileq3("flujo3.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(q_mesh, fileq3);
    
    //
    /// Computational multiphysics mesh reconstruction
    TPZVec<int> & active_approx_spaces = dfn.ExtractActiveApproxSpaces(mp_initial_cmesh);
    TPZCompMesh * dfn_hybrid_cmesh = dfn.DuplicateMultiphysicsCMeshMaterials(mp_initial_cmesh);
    dfn.CleanUpMultiphysicsCMesh(mp_initial_cmesh);
    
    int target_dim=3;
    int fractures_dim =2;
    int fractures_intersections_dim =1;
    dfn.InsertMaterialsForHibridization(target_dim, dfn_hybrid_cmesh, flux_trace_id, lagrange_id, flux_resistivity_id, mp_nterface_id);
    dfn.InsertMaterialsForMixedOperatorOnFractures(fractures_dim,dfn_hybrid_cmesh);
    dfn.InsertMaterialsForMixedOperatorOnFractures(fractures_intersections_dim,dfn_hybrid_cmesh);
    dfn.InsertMaterialsForMixedOperatorOnFractures(0,dfn_hybrid_cmesh);
    
   
    TPZManVector<TPZMaterial *> zeroDimMat(dfn.m_fracture_ids[0].size());
    int count = 0;
    for (auto matid:dfn.m_fracture_ids[0]) {
        zeroDimMat[count++] = dfn_hybrid_cmesh->FindMaterial(matid);
        dfn_hybrid_cmesh->MaterialVec().erase(matid);
    }
    dfn.BuildMultiphysicsCMesh(target_dim,dfn_hybrid_cmesh,active_approx_spaces,mp_mesh_vector);
    for (auto matptr:zeroDimMat) {
        dfn_hybrid_cmesh->InsertMaterialObject(matptr);
    }
//    std::set<int> MaterialIDs;
//    MaterialIDs.insert(-1942);
//    mp_mesh_vector[0]->AutoBuild(MaterialIDs);
    dfn.CreateInterfaceElements(target_dim, mp_nterface_id, dfn_hybrid_cmesh, mp_mesh_vector);
    dfn.CreateInterfaceElements(fractures_dim, mp_nterface_id, dfn_hybrid_cmesh, mp_mesh_vector);
    dfn.CreateInterfaceElements(fractures_intersections_dim, mp_nterface_id, dfn_hybrid_cmesh, mp_mesh_vector);
    dfn_hybrid_cmesh->ComputeNodElCon();
    
    if(1){
        int64_t nel = dfn_hybrid_cmesh->NElements();
        for(int64_t el=0; el<nel; el++)
        {
            TPZCompEl *cel = dfn_hybrid_cmesh->Element(el);
            TPZMultiphysicsElement *mpcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
            if(!mpcel) continue;
            if (mpcel->Dimension()==0) {
                std::cout<<" aqui!: "<<std::endl;
            }
            if(!mpcel->Element(0) && mpcel->Element(1))
            {
                
                TPZConnect &c = mpcel->Connect(0);
                int nelcon = c.NElConnected();
                if(nelcon <= 3)
                {
                    continue;
                }
                std::cout << "this element has a pressure but no flux\n";
                
                TPZGeoEl *gel = cel->Reference();
                // find a neighbour with matid = 8
                TPZGeoEl *gel8 = 0;
                int matid_target = 0;
                TPZGeoElSide gelside(gel,gel->NSides()-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(dfn.m_fracture_ids[0].find(neighbour.Element()->MaterialId()) != dfn.m_fracture_ids[0].end())
                    {
                        gel8 = neighbour.Element();
                        matid_target = neighbour.Element()->MaterialId();
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                /// if there is no element 8 found, issue a warning because the mesher should have created one
                if(!gel8)
                {
                    std::cout << "no element found at the intersection of fractures\n";
                    int matid = *dfn.m_fracture_ids[0].begin();
                    std::cout << "changin matid to " << matid << std::endl;
                    gel->SetMaterialId(matid);
                }
                else
                {
                    gel8->SetMaterialId(0);
                    gel->SetMaterialId(matid_target);
                }
            }
        }
    }
    
    dfn_hybrid_cmesh->InitializeBlock();
    //
    
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(dfn_hybrid_cmesh);
    
    bool solve_dfn_problem_Q = true;
    if (solve_dfn_problem_Q) {
        TPZManVector<TPZCompMesh * > mesh_vec = mp_cmesh->MeshVector();
        
        std::ofstream printpress("pres.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(mesh_vec[1], printpress);

        
        std::cout << "Condensing DFN equations." << std::endl;
        std::cout << "DFN neq before condensation = " << mp_cmesh->NEquations() << std::endl;
      //  dfn.GroupElements(mp_cmesh);
        std::cout << "DFN neq = " << mp_cmesh->NEquations() << std::endl;
        
#ifdef PZDEBUG
        {
            
            std::ofstream file_hybrid_condensed("Hybrid_mixed_condensed.txt");
            mp_cmesh->ComputeNodElCon();
            mp_cmesh->Print(file_hybrid_condensed);
        }
#endif
        
        std::cout << "DFN equations are condensed." << std::endl;
        
        TPZAnalysis *an = CreateAnalysis(mp_cmesh, sim);
        std::cout << "Assembly DFN problem neq = " << mp_cmesh->NEquations() << std::endl;
        an->Assemble();
        std::cout << "Assembly for DFN complete." << std::endl;
        
        std::cout << "Solving DFN problem." << std::endl;
        an->Solve();
        std::cout << "DFN problem solved." << std::endl;
        mp_cmesh->LoadSolutionFromMultiPhysics();
        
        
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        scalnames.Push("p");
        
        int div = 0;
        std::set<int> mat_id_3D;
        mat_id_3D.insert(1);
        mat_id_3D.insert(2);
        std::string file_reservoir("cube.vtk");
        an->DefineGraphMesh(3,mat_id_3D,scalnames,vecnames,file_reservoir);
        an->PostProcess(div,3);
        
        std::set<int> mat_id_2D;
        mat_id_2D.insert(26);
        std::string file_frac("fracture.vtk");
        an->DefineGraphMesh(2,mat_id_2D,scalnames,vecnames,file_frac);
        an->PostProcess(div,2);
        
        std::set<int> mat_id_1D;
        mat_id_1D.insert(27);
        std::string file_frac_intersections("fracture_intersections.vtk");
        an->DefineGraphMesh(1,mat_id_1D,scalnames,vecnames,file_frac_intersections);
        an->PostProcess(div,1);
    }
    
    return;
    
}
