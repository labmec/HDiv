
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
TPZCompMesh * CMeshMixed(TPZGeoMesh * geometry, int p, SimulationCase sim_data, TPZVec<TPZCompMesh *> &meshvec);
TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, int p, SimulationCase sim_data);
TPZCompMesh * FluxMesh(TPZGeoMesh * gmesh, int order, SimulationCase sim);
TPZCompMesh * PressureMesh(TPZGeoMesh * gmesh, int order,SimulationCase sim);
TPZAnalysis * CreateAnalysis(TPZCompMesh * cmesh, SimulationCase & sim_data);
void forcing(const TPZVec<REAL> &p, TPZVec<STATE> &f);
void UnwrapMesh(TPZCompMesh *cmesh);

/// Executes case 1
void Case_1();

/// Executes case 2
void Case_2();

int main(){
//    Case_1();
    Case_2();
}

void Case_1(){
    
    DebugStop();
}

#define NewMPCmesh_Q

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
#ifdef NewMPCmesh_Q
    cmixedmesh = MPCMeshMixed(gmesh, 1, sim);
#else
    cmixedmesh = CMeshMixed(gmesh, 1, sim, meshvec);
#endif
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

TPZCompMesh * FluxMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    int dimension = 3;
    int dirichlet = 0;
    int neumann = 1;
    int nvolumes = sim_data.omega_ids.size();
    int nboundaries = sim_data.gamma_ids.size();
    
//#ifdef PZDEBUG
//    if (nvolumes != 1) {
//        std::cout << "Error:: unable to compute the given case = " << &sim_data << std::endl;
//        DebugStop();
//    }
//#endif
    
   
    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    
        TPZMixedPoisson * volume1 = new TPZMixedPoisson(sim_data.omega_ids[0],3);
  
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(forcing, 5);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume1->SetForcingFunction(rhs);
    
    volume1->SetPermeability(1.0);
        cmesh->InsertMaterialObject(volume1);
   
    TPZMixedPoisson * volume2 = new TPZMixedPoisson(sim_data.omega_ids[1],3);
    
    volume2->SetForcingFunction(rhs);
    volume1->SetPermeability(1.0);
    cmesh->InsertMaterialObject(volume2);
    
    
    TPZMaterial * face_1 = volume1->CreateBC(volume1,sim_data.gamma_ids[0],neumann,val1,val2);
    
    cmesh->InsertMaterialObject(face_1);
    
    val2(0,0)=1000;
    TPZMaterial * face_2 = volume1->CreateBC(volume1,sim_data.gamma_ids[1],dirichlet,val1,val2);
    
    cmesh->InsertMaterialObject(face_2);
    
    val2(0,0)=-100;
    TPZMaterial * face_3 = volume1->CreateBC(volume1,sim_data.gamma_ids[2],neumann,val1,val2);
    
    cmesh->InsertMaterialObject(face_3);
    
        
    
    
    
    cmesh->SetDimModel(dimension);
    cmesh->SetDefaultOrder(order);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->AutoBuild();
    
        cmesh->SetDefaultOrder(order);
    
        
        cmesh->AutoBuild();
    
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "q_cmesh" << ".txt";
    std::ofstream sout(file_name.str().c_str());
    cmesh->Print(sout);
#endif
    
    return cmesh;
    
}
TPZCompMesh * PressureMesh(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
    
    int dimension = 3;
    int nvolumes = sim_data.omega_ids.size();
    

    
    TPZCompMesh *cmesh = new TPZCompMesh(geometry);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    TPZMaterial * volume = new TPZMatPoisson3d(sim_data.omega_ids[0]);
    TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(forcing, 5);
    rhs_exact->SetPolynomialOrder(sim_data.int_order);
    TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
    volume->SetForcingFunction(rhs);
        
     cmesh->InsertMaterialObject(volume);
  
    TPZMaterial * volume2 = new TPZMatPoisson3d(sim_data.omega_ids[1]);
    
    volume2->SetForcingFunction(rhs);
    cmesh->InsertMaterialObject(volume2);
        

    
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
    f[0]=x*y*z;
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
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(forcing, 5);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
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
    
    
    TPZManVector<TPZCompMesh * ,2> meshvector(2);
    meshvector[0] = FluxMesh(geometry, order, sim_data);
    meshvector[1] = PressureMesh(geometry, order, sim_data);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, cmesh);
    
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
    
    meshvec = meshvector;
    return cmesh;
}

TPZMultiphysicsCompMesh * MPCMeshMixed(TPZGeoMesh * geometry, int order, SimulationCase sim_data){
    
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
    
    TPZManVector<TPZCompMesh * ,2> meshvec(2);
    meshvec[0] = FluxMesh(geometry, order, sim_data);
    meshvec[1] = PressureMesh(geometry, order, sim_data);
    TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(geometry,meshvec);
    
    TPZFMatrix<STATE> val1(dimension,dimension,0.0),val2(dimension,1,0.0);
    
    for (int ivol=0; ivol<nvols; ivol++) {
        TPZMixedPoisson * volume = new TPZMixedPoisson(sim_data.omega_ids[ivol],dimension);
        volume->SetPermeability(sim_data.permeabilities[ivol]);
        
        TPZDummyFunction<STATE> * rhs_exact = new TPZDummyFunction<STATE>(forcing, 5);
        rhs_exact->SetPolynomialOrder(sim_data.int_order);
        TPZAutoPointer<TPZFunction<STATE> > rhs = rhs_exact;
        volume->SetForcingFunction(rhs);
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
    
#ifdef PZDEBUG
    std::stringstream file_name;
    file_name   << sim_data.dump_folder << "/" << "Dual_cmesh" << ".txt";
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

