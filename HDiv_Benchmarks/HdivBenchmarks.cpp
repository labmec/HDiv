//
//  Hdiv3DCurved.cpp
//  Publication about 3D Curved piola mapping for mixed formulation


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


#ifdef USING_BOOST
#include "boost/date_time/posix_time/posix_time.hpp"
#endif

void InsertFrac(TPZGeoMesh *gmesh, TPZFMatrix<REAL> corners, int matids );
TPZGeoMesh * case2mesh();

int main(){
    
    
    TPZGeoMesh *gmesh = case2mesh();
    std::ofstream file("testCase2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    
    TPZGeoMesh *geomesh = new TPZGeoMesh;
    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetfDimensionlessL(s);
    
    geomesh = Geometry.GeometricGmshMesh("case_1_1k.msh");
    
    std::ofstream file2("testCase1.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, file2);
    
    
    
    
    
    return 0;
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
    gengrid.SetBC(gmesh, 4, -2);
    gengrid.SetBC(gmesh, 5, -3);
    gengrid.SetBC(gmesh, 6, -4);
    gengrid.SetBC(gmesh, 7, -5);
    
    //    TPZVec<TPZGeoEl*> indexf;
    //    gmesh->Element(0)->Divide(indexf);
    //    gmesh->Element(1)->Divide(indexf);
    //    gmesh->Element(2)->Divide(indexf);
    //    gmesh->Element(3)->Divide(indexf);
    //    gmesh->Element(4)->Divide(indexf);
    //    gmesh->Element(12)->Divide(indexf);
    //    gmesh->Element(13)->Divide(indexf);
    //    gmesh->Element(7)->Divide(indexf);
    //    gmesh->Element(21)->Divide(indexf);
    
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
    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(8,-1,-6);
    gmesh3d->BuildConnectivity();
    
    
    int nels = gmesh3d->NElements();
    for (int iel =0; iel<nels; iel++){
        TPZGeoEl * gel = gmesh3d->Element(iel);
        if(!gel){continue;}
        if (gel->Dimension() == 2){continue;}
        TPZManVector<REAL,3> qsi(gel->Dimension()), xCenter(3,0.);
        gel->CenterPoint(gel->NSides()-1, qsi);
        gel->X(qsi,xCenter);
        if (xCenter[0]>0.5 && xCenter[1]<0.5){
            gel->SetMaterialId(2);
        }
        
        if ((xCenter[0]>0.625 && xCenter[0]<0.75 ) && (xCenter[1]>0.5 && xCenter[1]<0.625)&&(xCenter[2]>0.5 && xCenter[2]<0.75)){
            gel->SetMaterialId(2);
        }
        
        if (xCenter[0]>0.75 && (xCenter[1]>0.5 && xCenter[1]<0.75)&&(xCenter[2]>0.5)){
            gel->SetMaterialId(2);
        }
        
        if (xCenter[0]>0.875 && xCenter[1]>0.875 && xCenter[2]>0.875 ){
            gel->SetMaterialId(-7);
        }
        
        if (xCenter[0]<0.25 && xCenter[1]<0.25 && xCenter[2]<0.25 ){
            gel->SetMaterialId(-8);
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
    
    
    
    InsertFrac(gmesh3d,frac1,3);
    InsertFrac(gmesh3d,frac2,4);
    InsertFrac(gmesh3d,frac3,5);
    InsertFrac(gmesh3d,frac4,6);
    InsertFrac(gmesh3d,frac5,7);
    InsertFrac(gmesh3d,frac6,8);
    
    return gmesh3d;
}
