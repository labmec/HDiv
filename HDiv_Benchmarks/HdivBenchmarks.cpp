
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
    
    SimulationCase() : IsMHMQ(false), UsePardisoQ(true), IsHybrid(false),UseFrontalQ(false), UseGmshMeshQ(false), NonAffineQ(false), elemen_type(0), n_h_levels(0), n_p_levels(1), n_acc_terms(0), int_order(1), n_threads(0),perturbation_type(0), mesh_type(""), domain_type(""),conv_summary(""),dump_folder(""),omega_ids(),omega_dim(),gamma_ids(), gamma_dim(), permeabilities(),porosities(), type(), vals()
    {
        
    }
    
    SimulationCase(const SimulationCase &copy) :IsHybrid(copy.IsHybrid) ,IsMHMQ(copy.IsMHMQ), UsePardisoQ(copy.UsePardisoQ), UseFrontalQ(copy.UseFrontalQ),
    UseGmshMeshQ(copy.UseGmshMeshQ), NonAffineQ(copy.NonAffineQ), elemen_type(copy.elemen_type), n_h_levels(copy.n_h_levels), n_p_levels(copy.n_p_levels), n_acc_terms(copy.n_acc_terms), int_order(copy.int_order),n_threads(copy.n_threads),perturbation_type(copy.perturbation_type), mesh_type(copy.mesh_type), domain_type(copy.domain_type), conv_summary(copy.conv_summary),
    dump_folder(copy.dump_folder), omega_ids(copy.omega_ids), omega_dim(copy.omega_dim), gamma_ids(copy.gamma_ids), gamma_dim(copy.gamma_dim),
    permeabilities(copy.permeabilities),porosities(copy.porosities),
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
        omega_dim = copy.omega_dim;
        gamma_ids = copy.gamma_ids;
        gamma_dim = copy.gamma_dim;
        permeabilities=copy.permeabilities;
        porosities = copy.porosities;
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
TPZCompMesh * SaturationMesh(TPZVec<TPZCompMesh *> &meshvec, int order, SimulationCase sim);
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
TPZAnalysis * CreateTransportAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
void forcing(const TPZVec<REAL> &p, TPZVec<STATE> &f);
void SeparateConnectsByFracId(TPZCompMesh * mixed_cmesh,int fracid);
void InsertFractureMaterial(TPZCompMesh *cmesh);
TPZVec<REAL> MidPoint(TPZVec<REAL> & x_i, TPZVec<REAL> & x_e);
void AdjustMaterialIdBoundary(TPZMultiphysicsCompMesh *cmesh);
TPZCompMesh *CreateTransportMesh(TPZMultiphysicsCompMesh *cmesh);

void InsertTransportInterfaceElements(TPZMultiphysicsCompMesh *cmesh);
TPZMultiphysicsCompMesh * MPTransportMesh(TPZMultiphysicsCompMesh * mixed, TPZStack<TFracture> & fracture_data ,SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);
void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC);

void UniformRefinement(TPZGeoMesh * geometry, int h_level);

void InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes);

void TimeFoward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt);

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
#endif
    

    Pretty_cube();

//    Case_1();

//     Case_2();

}

/// Executes cube
void Pretty_cube(){
    
    int h_level = 0;
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.n_threads = 0;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(0.25);
    
    /// not used but inserted
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    sim.porosities.push_back(0.25);
    
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
    fracture.m_id               = 6;
    fracture.m_dim              = 2;
    fracture.m_kappa_normal     = 0.001;
    fracture.m_kappa_tangential = 0.001;
    fracture.m_d_opening        = 1.0e-2;
    fracture.m_porosity         = 0.25;
    fracture_data.push_back(fracture);
    fracture.m_id               = 7;
    fracture.m_dim              = 1;
    fracture.m_kappa_normal     = 0.001;
    fracture.m_kappa_tangential = 0.001;
    fracture.m_d_opening        = 1.0e-2;
    fracture.m_porosity         = 0.25;
    fracture_data.push_back(fracture);
    fracture.m_id               = 8;
    fracture.m_dim              = 0;
    fracture.m_kappa_normal     = 0.00;
    fracture.m_kappa_tangential = 0.00;
    fracture.m_d_opening        = 1.0e-2;
    fracture.m_porosity         = 0.25;
    fracture_data.push_back(fracture);

    
    /// Benchmarks Material ID convention
    /// 1 and 2 for 3D matrix
    /// 3,4, and 5 for 2D matrix boundaries (3 -> inlet, 4 -> outlet, 5 -> impervious)
    /// 6 fractures
    /// 7 fractures intersections
    /// 8 crossing intersections
    TPZManVector<std::map<std::string,int>,5> dim_name_and_physical_tag(4); // From 0D to 3D
    dim_name_and_physical_tag[3]["RockMatrix_1"] = 1;
    dim_name_and_physical_tag[3]["RockMatrix_2"] = 2;
    dim_name_and_physical_tag[2]["BCInlet"] = 3;
    dim_name_and_physical_tag[2]["BCOutlet"] = 4;
    dim_name_and_physical_tag[2]["BCImpervious"] = 5;
    dim_name_and_physical_tag[2]["Fractures"] = 6;
    dim_name_and_physical_tag[1]["FracturesIntersections"] = 7;
    dim_name_and_physical_tag[0]["CrossingIntresections"] = 8;
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
    std::string file_gmsh = source_dir + "/meshes/the_cuttest_cube/cube.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    
    Geometry.SetDimNamePhysical(dim_name_and_physical_tag);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
    UniformRefinement(gmesh, h_level);
    
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
    if(sim.IsHybrid){
        
        THybridizeDFN dfn_hybridzer;
        dfn_hybridzer.SetFractureData(fracture_data);
        
        dfn_hybridzer.SetReservoirBoundaryData(bc_ids_2d);
        dfn_hybridzer.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
        dfn_hybridzer.SetMapReservoirBCToDFNBC0DIds(bc_ids_0d_map);
        cmeshm = dfn_hybridzer.Hybridize(cmixedmesh);

    }
    else{
        cmeshm=cmixedmesh;
    }

    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmeshm);

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
        
        TPZAnalysis *an = CreateAnalysis(mp_cmesh, sim);
        std::cout << "Assembly neq = " << mp_cmesh->NEquations() << std::endl;
        an->Assemble();
        
        std::cout << "Solution of the system" << std::endl;
        an->Solve();
        
        mp_cmesh->LoadSolutionFromMultiPhysics();
        
#ifdef PZDEBUG
        std::ofstream file_geo_hybrid("geometry_cube_hybrid.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(cmeshm->Reference(), file_geo_hybrid);
#endif
        
        TPZStack<std::string,10> scalnames, vecnames;
        vecnames.Push("q");
        vecnames.Push("kappa");
        scalnames.Push("p");
        
        int div = 0;
        std::string file_reservoir("cube.vtk");
        an->DefineGraphMesh(3,scalnames,vecnames,file_reservoir);
        an->PostProcess(div,3);
        
#ifdef PZDEBUG
        { /// fracture postprocessor
            TPZStack<std::string,10> scalnames, vecnames;
            scalnames.Push("state");
            std::string file_frac("fracture.vtk");
            auto material = mesh_vec[1]->FindMaterial(6);
            TPZL2Projection * fract_2d = dynamic_cast<TPZL2Projection *>(material);
            fract_2d->SetDimension(2);
            TPZAnalysis frac_an(mesh_vec[1],false);
            frac_an.DefineGraphMesh(2,scalnames,vecnames,file_frac);
            frac_an.PostProcess(div,2);
        }
        
        { /// lagrange postprocessor
            TPZStack<std::string,10> scalnames, vecnames;
            scalnames.Push("state");
            std::string file_frac("lagrange_1d.vtk");
            auto material = mesh_vec[1]->FindMaterial(7);
            TPZL2Projection * fract_2d = dynamic_cast<TPZL2Projection *>(material);
            fract_2d->SetDimension(1);
            TPZAnalysis frac_an(mesh_vec[1],false);
            frac_an.DefineGraphMesh(1,scalnames,vecnames,file_frac);
            frac_an.PostProcess(div,1);
        }
#endif
    }

    TPZCompMesh *s_cmesh = CreateTransportMesh(mp_cmesh);
    TPZManVector<TPZCompMesh *,3> meshtrvec(3);
    meshtrvec[0] = meshvec[0];
    meshtrvec[1] = meshvec[1];
    meshtrvec[2] = s_cmesh;
    
    TPZMultiphysicsCompMesh *cmesh_transport = MPTransportMesh(mp_cmesh, fracture_data, sim, meshtrvec);
    TPZAnalysis * tracer_analysis = CreateTransportAnalysis(cmesh_transport, sim);
    
    int n_steps = 10;
    REAL dt     = 0.1;
    TimeFoward(tracer_analysis, n_steps, dt);
    
    return;

}

void TimeFoward(TPZAnalysis * tracer_analysis, int & n_steps, REAL & dt){
    
    
    TPZMultiphysicsCompMesh * cmesh_transport = dynamic_cast<TPZMultiphysicsCompMesh *>(tracer_analysis->Mesh());
    
    if (!cmesh_transport) {
        DebugStop();
    }
    
    /// Compute mass matrix M.
    TPZAutoPointer<TPZMatrix<STATE> > M;
    {
        
        
        bool mass_matrix_Q = true;
        std::set<int> volumetric_mat_ids = {1,2,6,7,8};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                DebugStop();
            }
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing Mass Matrix." << std::endl;
        tracer_analysis->Assemble();
        M = tracer_analysis->Solver().Matrix()->Clone();
    }
    
    {
        bool mass_matrix_Q = false;
        std::set<int> volumetric_mat_ids = {1,2,6,7,8};
        
        for (auto mat_id: volumetric_mat_ids) {
            TPZMaterial * mat = cmesh_transport->FindMaterial(mat_id);
            TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
            if (!volume) {
                DebugStop();
            }
            volume->SetTimeStep(dt);
            volume->SetMassMatrixAssembly(mass_matrix_Q);
        }
        
        std::cout << "Computing transport operator K = M + T " << std::endl;
        tracer_analysis->Assemble();
    }
    
    /// Time evolution
    std::string file_reservoir("transport_cube.vtk");
    
    {
        int div = 0;
        int64_t n_eq = tracer_analysis->Mesh()->NEquations();
        TPZFMatrix<REAL> s_n(n_eq,1,0.0);
        TPZFMatrix<REAL> last_state_mass;
        TPZFMatrix<REAL> ds,s_np1;
        
        for (int i = 0; i < n_steps; i++) {
            M->Multiply(s_n, last_state_mass);
            
            tracer_analysis->Rhs() -= last_state_mass;
            tracer_analysis->Rhs() *= -1.0;
            
            tracer_analysis->Solve(); /// (LU decomposition)
            ds = tracer_analysis->Solution();
            s_np1 = s_n + ds;
            tracer_analysis->LoadSolution(s_np1);
            cmesh_transport->LoadSolutionFromMultiPhysics();
            tracer_analysis->AssembleResidual();
            
            /// postprocess ...
            TPZStack<std::string,10> scalnames, vecnames;
            scalnames.Push("Sw");
            scalnames.Push("So");
            
            std::map<int,int> volumetric_ids;
            volumetric_ids.insert(std::make_pair(1, 3));
            
            std::map<int,int> fracture_ids;
            fracture_ids.insert(std::make_pair(6, 2));
            
            std::map<int,int> fracture_intersections_ids;
            fracture_intersections_ids.insert(std::make_pair(7, 1));
            
            for (auto data: volumetric_ids) {
                TPZMaterial * mat = cmesh_transport->FindMaterial(data.first);
                TPZTracerFlow * volume = dynamic_cast<TPZTracerFlow * >(mat);
                if (!volume) {
                    DebugStop();
                }
                volume->SetDimension(data.second);
            }
            int dim = 3;
            tracer_analysis->DefineGraphMesh(dim,scalnames,vecnames,file_reservoir);
            tracer_analysis->PostProcess(div,dim);
            
            // configuring next time step
            s_n = s_np1;
        }
        
        
    }
    
}


void Case_1(){
    
    int h_level = 0;
    
    SimulationCase sim;
    sim.UsePardisoQ=true;
    sim.IsHybrid=true;
    sim.omega_ids.push_back(1);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    
    sim.omega_ids.push_back(2);
    sim.omega_dim.push_back(3);
    sim.permeabilities.push_back(1.0);
    
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
    
    
    int bc_1d_non_flux = -1942;
    std::map<int,int> bc_ids_1d_map;
    bc_ids_1d_map.insert(std::make_pair(bc_non_flux,bc_1d_non_flux));

    
    /// Defining DFN data
    TPZStack<TFracture> fracture_data;
    TFracture fracture;
    fracture.m_id               = 6;
    fracture.m_dim              = 2;
    fracture.m_kappa_normal     = 0.001;
    fracture.m_kappa_tangential = 1000.0;
    fracture.m_d_opening        = 1.0e-2;
    fracture_data.push_back(fracture);
   
    
    TPZGmshReader Geometry;
    std::string source_dir = SOURCE_DIR;
    std::string file_gmsh = source_dir + "/meshes/Case_1/case_1.msh";
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    std::string version("4.1");
    Geometry.SetFormatVersion(version);
    gmesh = Geometry.GeometricGmshMesh(file_gmsh.c_str());
    Geometry.PrintPartitionSummary(std::cout);
    
    UniformRefinement(gmesh, h_level);
    
#ifdef PZDEBUG
    std::ofstream file("geometry_case_1.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    
    std::ofstream file_txt("geometry_case_1_base.txt");
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
    if(sim.IsHybrid){
        
        THybridizeDFN dfn_hybridzer;
        dfn_hybridzer.SetFractureData(fracture_data);
        
        dfn_hybridzer.SetReservoirBoundaryData(bc_ids_2d);
        dfn_hybridzer.SetMapReservoirBCToDFNBC1DIds(bc_ids_1d_map);
        cmeshm = dfn_hybridzer.Hybridize(cmixedmesh);
        
    }
    else{
        cmeshm=cmixedmesh;
    }
    
    TPZMultiphysicsCompMesh * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmeshm);
    
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
    
    
    TPZAnalysis *an = CreateAnalysis(cmeshm, sim);
    std::cout << "Assembly neq = " << cmeshm->NEquations() << std::endl;
    an->Assemble();
    
    std::cout << "Solution of the system" << std::endl;
    an->Solve();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mesh_vec, cmeshm);
    
#ifdef PZDEBUG
    std::ofstream file_geo_hybrid("geometry_case_1_hybrid.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(cmeshm->Reference(), file_geo_hybrid);
#endif
    
    TPZStack<std::string,10> scalnames, vecnames;
    vecnames.Push("q");
    vecnames.Push("kappa");
    scalnames.Push("p");
    
    int div = 0;
    std::string file_reservoir("case_1.vtk");
    an->DefineGraphMesh(3,scalnames,vecnames,file_reservoir);
    an->PostProcess(div,3);
    
#ifdef PZDEBUG
    { /// fracture postprocessor
        TPZStack<std::string,10> scalnames, vecnames;
        scalnames.Push("state");
        std::string file_frac("fracture.vtk");
        auto material = mesh_vec[1]->FindMaterial(6);
        TPZL2Projection * fract_2d = dynamic_cast<TPZL2Projection *>(material);
        fract_2d->SetDimension(2);
        TPZAnalysis frac_an(mesh_vec[1],false);
        frac_an.DefineGraphMesh(2,scalnames,vecnames,file_frac);
        frac_an.PostProcess(div,2);
    }
#endif
    
    TPZCompMesh *cmesh_transport = CreateTransportMesh(mp_cmesh);
    TPZManVector<TPZCompMesh *,3> meshtrvec(3);
    meshtrvec[0] = meshvec[0];
    meshtrvec[1] = meshvec[1];
    meshtrvec[2] = cmesh_transport;
    TPZMultiphysicsCompMesh *mp_tr_cmesh = MPTransportMesh(mp_cmesh, fracture_data, sim, meshtrvec);
    DebugStop();
    InsertTransportInterfaceElements(mp_tr_cmesh);
    
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
                TPZGeoElBC gbc(gelside, fracture_id+1);
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

void AdjustMaterialIdBoundary(TPZMultiphysicsCompMesh *cmesh)
{
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++)
    {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = 0;
        if(cel) gel = cel->Reference();
        if(gel->MaterialId() != -1942) continue;
        int dim = gel->Dimension();
        TPZManVector<REAL,3> xi(dim),x(3,0.);
        gel->CenterPoint(gel->NSides()-1, xi);
        gel->X(xi, x);
        if(fabs(x[2]) < 1.e-6)
        {
            std::cout << "Changing matid of element at x = " << x << " to -5\n";
            gel->SetMaterialId(-5);
        }
        if(fabs(x[2]-1.) < 1.e-6)
        {
            std::cout << "Changing matid of element at x = " << x << " to -6\n";
            gel->SetMaterialId(-6);
        }
        
    }
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
    
    for (int dim=3; dim>0; dim--) {
        cmesh->Reference()->ResetReference();
        cmesh->SetDimModel(dim);
        cmesh->SetDefaultOrder(order);
        std::set<int> matids;
        for (int ivol=0; ivol<nvols; ivol++) {
            if(sim_data.omega_dim[ivol]==dim) matids.insert(sim_data.omega_ids[ivol]);
        }
        for (int ibound=0; ibound<nbound; ibound++) {
            if(sim_data.gamma_dim[ibound]==dim) matids.insert(sim_data.gamma_ids[ibound]);
        }
        cmesh->SetAllCreateFunctionsHDiv();
        cmesh->AutoBuild(matids);
    }
    cmesh->InitializeBlock();
    
#ifdef PZDEBUG2
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
    
#ifdef PZDEBUG2
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

        TPZMixedDarcyFlow * volume = new TPZMixedDarcyFlow(sim_data.omega_ids[ivol],dim);

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
    
#ifdef PZDEBUG2
    std::stringstream file_name;
    file_name  << "Dual_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    meshvec = mesh_vec;
    return cmesh;
}

TPZMultiphysicsCompMesh * MPTransportMesh(TPZMultiphysicsCompMesh * mixed, TPZStack<TFracture> & fracture_data ,SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec){
    
    
    TPZGeoMesh *geometry = mixed->Reference();
    int dimension = geometry->Dimension();
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry);

    
    std::set<int> matrix_mat_ids;// = {1,2,6,7,8};
    std::set<int> fracture_mat_ids;
    

    /// Inserting matrix materials
    int n_vols = sim_data.omega_ids.size();
    for (int i = 0; i < n_vols; i++) {
        int mat_id = sim_data.omega_ids[i];
        REAL phi = sim_data.porosities[i];
        TPZTracerFlow * volume = new TPZTracerFlow(mat_id,0);
        volume->SetPorosity(phi);
        cmesh->InsertMaterialObject(volume);
    }
    
    /// Inserting fracture materials
    int n_fracs = fracture_data.size();
    for (int i = 0; i < n_fracs; i++) {
        int mat_id = fracture_data[i].m_id;
        REAL phi = fracture_data[i].m_porosity;
        TPZTracerFlow * volume = new TPZTracerFlow(mat_id,0);
        volume->SetPorosity(phi);
        cmesh->InsertMaterialObject(volume);
    }
    
    std::set<int> bc_inlet_mat_ids = {3,5,130,150,230,250};
    std::set<int> bc_outlet_mat_ids = {4,140,240};
    
    TPZMaterial * material = cmesh->FindMaterial(1);
    if (!material) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0);
    TPZFMatrix<STATE> val2(1,1,0);
    
    int typ_inlet = 0; // inlet
    /// Inserting the materials
    val2(0,0) = 1.0;
    for (auto mat_id: bc_inlet_mat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ_inlet, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    
    int typ_outlet = 1; // outlet
    /// Inserting the materials
    val2(0,0) = 0.0;
    for (auto mat_id: bc_outlet_mat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ_outlet, val1, val2);
        cmesh->InsertMaterialObject(bc);
    }
    
    cmesh->SetDimModel(dimension);
    
    TPZManVector<int,5> active_approx_spaces(3); /// 1 stands for an active approximation spaces
    active_approx_spaces[0] = 0;
    active_approx_spaces[1] = 0;
    active_approx_spaces[2] = 1;
    cmesh->BuildMultiphysicsSpace(active_approx_spaces,meshvec);
    
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
    
    InsertTransportInterfaceElements(cmesh);
    
#ifdef PZDEBUG
    std::ofstream transport("transport_cmesh.txt");
    cmesh->Print(transport);
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

TPZAnalysis * CreateTransportAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data){
    
    TPZAnalysis * analysis = new TPZAnalysis(cmesh, true);
    
    if (sim_data.UsePardisoQ) {
        
        TPZSpStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        analysis->SetStructuralMatrix(matrix);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        
        return analysis;
    }else{
        
        TPZSkylineNSymStructMatrix matrix(cmesh);
        matrix.SetNumThreads(sim_data.n_threads);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELU);
        analysis->SetSolver(step);
        analysis->SetStructuralMatrix(matrix);
        return analysis;
        
    }
    
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
    
//    TPZExtendGridDimension extend(gmesh,1.0);
//    extend.SetElType(1);
//    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(1,-5,-6);
//    gmesh3d->BuildConnectivity();
   
    
    std::ofstream file2("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file2);
    
    std::ofstream file3("gmesh.txt");
    gmesh->Print(file3);
    
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
//    casetest.gamma_ids.push_back(-5);
//    casetest.gamma_ids.push_back(-6);
    //    d = 0;
    //    n = 1;
    casetest.type.push_back(1);
    casetest.type.push_back(0);
    casetest.type.push_back(1);
    casetest.type.push_back(0);
//    casetest.type.push_back(0);
//    casetest.type.push_back(0);
    
    casetest.vals.push_back(0.0);
    casetest.vals.push_back(2.0);
    casetest.vals.push_back(0.0);
    casetest.vals.push_back(1.0);
//    casetest.vals.push_back(2.0);
//    casetest.vals.push_back(1.0);
    
    
    TPZCompMesh *cmeshq = FluxMesh(gmesh, 1, casetest);

    
    //
    TPZCompMesh *cmeshp = PressureMesh(gmesh, 1, casetest);
    TPZVec<TPZCompMesh *> fmeshvec(2);
    fmeshvec[0]=cmeshq;
    fmeshvec[1]=cmeshp;
    
    TPZCompMesh *cmixedmesh = MPCMeshMixed(gmesh, 1, casetest, fmeshvec);
    std::ofstream filemixed("mixedMesh.txt");
   // cmixedmesh->Print(filemixed);
    
    
    
    
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
  //  scalnames.Push("Permeability");
    
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


TPZCompMesh *CreateTransportMesh(TPZMultiphysicsCompMesh *cmesh)
{
    TPZCompMesh *q_cmesh = cmesh->MeshVector()[0];
//    TPZCompMesh *p_cmesh = cmesh->MeshVector()[1];
    TPZCompMesh *s_cmesh = new TPZCompMesh(q_cmesh->Reference());
    
    TPZGeoMesh * geometry = q_cmesh->Reference();
    if (!geometry) {
        DebugStop();
    }
    
    geometry->ResetReference();
    q_cmesh->LoadReferences();

    std::set<int> allmat_ids = {1,2,6,7,8};
    std::set<int> bcmat_ids = {3,4,5,130,140,150,230,240,250};
    
    int nstate = 1;
    TPZVec<STATE> sol(1,0.0);
    
    /// Inserting the materials
    for (auto mat_id: allmat_ids) {
        TPZL2Projection * volume = new TPZL2Projection(mat_id,0,nstate, sol);
        s_cmesh->InsertMaterialObject(volume);
    }
    
    TPZMaterial * material = s_cmesh->FindMaterial(1);
    if (!material) {
        DebugStop();
    }
    TPZFMatrix<STATE> val1(1,1,0);
    TPZFMatrix<STATE> val2(1,1,1);
    
    int typ = 0; // inlet
    /// Inserting the materials
    for (auto mat_id: bcmat_ids) {
        TPZMaterial * bc = material->CreateBC(material, mat_id, typ, val1, val2);
        s_cmesh->InsertMaterialObject(bc);
    }
    geometry->ResetReference();

    s_cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    int s_order = 0;
    s_cmesh->SetDefaultOrder(s_order);
    for (auto cel : q_cmesh->ElementVec()) {
        
        if (!cel) {
            DebugStop();
        }
        
        TPZGeoEl * gel = cel->Reference();
        int gel_index = gel->Index();
        int mat_id = gel->MaterialId();
        
        if (allmat_ids.find(mat_id) != allmat_ids.end()) {
            TPZGeoEl * gel = geometry->Element(gel_index);
            CreateTransportElement(s_order,s_cmesh, gel, false);
        }
        if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
            TPZGeoEl * gel = geometry->Element(gel_index);
            CreateTransportElement(s_order,s_cmesh, gel, true);
        }
    }


    
    s_cmesh->SetDimModel(0);
    s_cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    
    /// Create point element
    for (auto cel : cmesh->ElementVec()) {
        if (!cel) {
            continue;
        }
        
        TPZGeoEl * gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        
        if (gel->Dimension() != 0) {
            continue;
        }
        int matid = gel->MaterialId();
        if(bcmat_ids.find(matid) != bcmat_ids.end())
        {
            continue;
        }
        if(allmat_ids.find(matid) == allmat_ids.end())
        {
            continue;
        }
        
        TPZMultiphysicsInterfaceElement * int_face = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (int_face) {
            continue;
        }
        
        int n_connect = cel->NConnects();
        if (n_connect !=1 && gel->MaterialId() == 8) {
            DebugStop();
        }
        else if(n_connect == 1 && gel->MaterialId() != 8)
        {
            std::cout << "I dont understand matid should be 8\n";
        }
        else if (n_connect == 0)
        {
            DebugStop();
        }
        TPZConnect & c = cel->Connect(0);
        
        if (c.NElConnected() > 2 && gel->MaterialId() == 8) {
            CreateTransportElement(s_order,s_cmesh, gel, false);
        }
        
        
    }
    geometry->ResetReference();
    s_cmesh->InitializeBlock();
    return s_cmesh;
}

void InsertTransportInterfaceElements(TPZMultiphysicsCompMesh *cmesh)
{
    
    int transport_matid = 1;
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int64_t nel = cmesh->NElements();
    TPZManVector<std::vector<int>,4> gel_vol_index_to_cel(4);
   
    for (int64_t el = 0; el< nel; el++) {
        
        TPZCompEl *cel = cmesh->Element(el);
        if(!cel) DebugStop();
        TPZMultiphysicsElement *celmp = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if(!celmp) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        
        int gel_dim = gel->Dimension();
        gel_vol_index_to_cel[gel_dim].push_back(el);
        
    }

    InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[3]);
    InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[2]);
    InsertInterfacesBetweenElements(transport_matid, cmesh, gel_vol_index_to_cel[1]);
    
    cmesh->ComputeNodElCon();
}

void InsertInterfacesBetweenElements(int transport_matid, TPZCompMesh * cmesh, std::vector<int> & cel_indexes){
    
    TPZGeoMesh * geometry = cmesh->Reference();
    if (!geometry) {
        DebugStop();
    }
    
    int mesh_dim = geometry->Dimension();
    
    std::set<int> bcmat_ids = {3,4,5,130,140,150,230,240,250};
    bool needs_all_boundaries_Q = true;
    TPZManVector<int64_t,3> left_mesh_indexes(2,0);
    left_mesh_indexes[0] = 0;
    left_mesh_indexes[1] = 2;
    TPZManVector<int64_t,3> right_mesh_indexes(1,0);
    right_mesh_indexes[0] = 2;
    for (auto cel_index: cel_indexes) {
        TPZCompEl *cel = cmesh->Element(cel_index);
        TPZGeoEl *gel = cel->Reference();
        
        int nsides = gel->NSides();
        int geldim = gel->Dimension();
        for (int side = 0; side<nsides; side++) {
            int sidedim = gel->SideDimension(side);
            if(sidedim < geldim-1) continue;
            int mat_id = gel->MaterialId();
            if (bcmat_ids.find(mat_id) != bcmat_ids.end()) {
                continue;
            }
            TPZStack<TPZCompElSide> celstack;
            TPZCompElSide celside(cel,side);
            TPZGeoElSide gelside(gel,side);
            gelside.EqualLevelCompElementList(celstack, 0, 0);
            if(celstack.size() == 0 && sidedim < mesh_dim && needs_all_boundaries_Q)
            {
                DebugStop();
            }
            if(celstack.size() == 0){
                continue;
            }
            
            if(sidedim == geldim-1)
            {
                
                int count = 0;
                int stack_index = 0;
                for (int icel = 0; icel < celstack.size(); icel++) {
                    TPZCompEl * cel_neigh = celstack[icel].Element();
                    TPZMultiphysicsInterfaceElement * mp_interface_cel = dynamic_cast<TPZMultiphysicsInterfaceElement * >(cel_neigh);
                    int gel_neigh_dim = cel_neigh->Dimension();
                    if (!mp_interface_cel) {
                        if (gel_neigh_dim == geldim-1) {
                            count++;
                            stack_index = icel;
                        }
                    }
                }
                
                /// There must be a lower dimensional transport element to equate the fluxes
                if(count != 1){
                    DebugStop();
                }
                
                TPZGeoEl *neighgel = celstack[stack_index].Element()->Reference();
                if(geldim < neighgel->Dimension())
                {
                    // we only create interfaces from lower to higher dimensional elements
                    // neighgel must be a boundary element
                    continue;
                }
                if(geldim == neighgel->Dimension() && gel->Id() > neighgel->Id())
                {
                    // the interface is created from the lower id to the higher id
                    continue;
                }
                {
                    int mat_interface_id;
                    int neigh_mat_id = neighgel->MaterialId();
                    if (bcmat_ids.find(neigh_mat_id) != bcmat_ids.end()) {
                        mat_interface_id = neigh_mat_id;
                    }else{
                        mat_interface_id = transport_matid;
                    }
                    
                    // we need to create the interface
                    TPZGeoElBC gbc(gelside,mat_interface_id);
                    int64_t index;
                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstack[stack_index]);
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                }
            }
        }
    }
    
}

TPZCompMesh * SaturationMesh(TPZMultiphysicsCompMesh *mpcompmesh, int order, SimulationCase sim){
    
    TPZManVector<TPZCompMesh * > mesh_vec = mpcompmesh->MeshVector();
    TPZCompMesh *qMesh = mesh_vec[0];
    TPZCompMesh *pMesh = mesh_vec[1];
    mpcompmesh->MeshVector().Resize(3);
    
    int dim =3;
    pMesh->LoadReferences();
    
    TPZGeoMesh *gmesh = pMesh->Reference();
    TPZCompMesh *smesh = new TPZCompMesh(gmesh);
    int mat_id_s = 1;
   
    
    int nel = pMesh->NElements();
    TPZStack<TPZCompElSide> el_oned_connec;
    for(int iel =0; iel<nel; iel++){
        TPZGeoEl *gel = pMesh->Element(iel)->Reference();
        if (gel->Dimension()==0) {
            TPZGeoElSide gelside(gel,0);
            TPZStack<TPZCompElSide> elsidevec;
            int onlyinterpolated =0;
            int removeduplicates =0;
            gelside.EqualLevelCompElementList(elsidevec, onlyinterpolated,  removeduplicates);
            int nelconnect = elsidevec.size();
            for (int ioned=0; ioned<nelconnect; ioned++) {
                if (nelconnect>2) {
                    int64_t index;
                    gel->SetMaterialId(mat_id_s);
                    smesh->CreateCompEl(gel, index);
                }
            }
        }
    }
    
    qMesh->LoadReferences();
    qMesh->SetDimModel(dim);
    
    int nelq = qMesh->NElements();
    for (int Dims=3; Dims<1; Dims++) {
        for (int iel =0; iel<nelq; iel++) {
            TPZGeoEl *gel = qMesh->Element(iel)->Reference();
            if (!gel->Reference()) {
                continue;
            }
            int dimel = gel->Dimension();
            int nconnecs = gel->Reference()->NConnects();
            if (dimel == Dims && nconnecs ==1) {
                int64_t index;
                gel->SetMaterialId(mat_id_s);
                smesh->CreateCompEl(gel, index);
                
            }
        }
    }
    
    
    smesh->SetDimModel(dim);
    smesh->SetDefaultOrder(order);
    
    smesh->SetAllCreateFunctionsDiscontinuous();
    //smesh->AutoBuild();
    
    smesh->AdjustBoundaryElements();
    smesh->CleanUpUnconnectedNodes();
    
    mpcompmesh->MeshVector()[2]=smesh;
    
    return smesh;
    
   
}

void CreateTransportElement(int p_order, TPZCompMesh *cmesh, TPZGeoEl *gel, bool is_BC){
    int64_t cel_index;
    int dimension = gel->Dimension();
    cmesh->SetDimModel(dimension);
    if (is_BC) {
        cmesh->SetDimModel(dimension+1);
    }
    TPZCompEl * cel = cmesh->ApproxSpace().CreateCompEl(gel, *cmesh, cel_index);
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
    TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
    if (intel){
        intel->PRefine(p_order);
    } else if (intelDisc) {
        intelDisc->SetDegree(p_order);
        intelDisc->SetTrueUseQsiEta();
    } else {
        DebugStop();
    }
    gel->ResetReference();
}

void UniformRefinement(TPZGeoMesh * geometry, int h_level) {
    
    TPZManVector<TPZGeoEl*> sons;
    for(int i=0; i < h_level; i++)
    {
        int64_t nels = geometry->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = geometry->ElementVec()[elem];
            if(!gel || gel->HasSubElement())
                continue;
            gel->Divide(sons);
        }
    }
    geometry->ResetConnectivities();
    geometry->BuildConnectivity();
}
