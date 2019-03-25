//
//  THybridzeDFN.cpp
//  HDiv
//
//  Created by Omar Durán on 3/23/19.
//

#include "THybridzeDFN.h"


/// Default constructor
THybridzeDFN::THybridzeDFN() : TPZHybridizeHDiv(){
    
}

/// Default destructor
THybridzeDFN::~THybridzeDFN(){
    
}

void THybridzeDFN::Hybridize(TPZCompMesh * mp_cmesh, int target_dim, TPZStack<TFracture> & fracture_ids){
    
    
    
}

std::tuple<int, int> THybridzeDFN::DissociateConnects(int flux_trace_id, int lagrange_mult_id, const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> & mesh_vec){
    
    TPZCompMesh * flux_cmesh = mesh_vec[0];
    
    TPZGeoElSide gleft(left.Reference());
    TPZGeoElSide gright(right.Reference());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());
    TPZInterpolatedElement *intelright = dynamic_cast<TPZInterpolatedElement *> (right.Element());
    intelleft->SetSideOrient(left.Side(), 1);
    intelright->SetSideOrient(right.Side(), 1);
    TPZStack<TPZCompElSide> equalright;
    TPZConnect & cleft = intelleft->SideConnect(0, left.Side());
    
    if (cleft.HasDependency()) {
    
        gright.EqualLevelCompElementList(equalright,1,0);
    
#ifdef PZDEBUG
        if(equalright.size() > 1)
        {
            DebugStop();
        }
        if(equalright.size()==1)
        {
            TPZGeoEl *equalgel = equalright[0].Element()->Reference();
            if(equalgel->Dimension() != flux_cmesh->Dimension()-1)
            {
                DebugStop();
            }
        }
#endif
        // reset the reference of the wrap element
        if(equalright.size()) equalright[0].Element()->Reference()->ResetReference();
        cleft.RemoveDepend();
    }
    else
    {
        int64_t index = flux_cmesh->AllocateNewConnect(cleft);
        TPZConnect &newcon = flux_cmesh->ConnectVec()[index];
        cleft.DecrementElConnected();
        newcon.ResetElConnected();
        newcon.IncrementElConnected();
        newcon.SetSequenceNumber(flux_cmesh->NConnects() - 1);
        
        int rightlocindex = intelright->SideConnectLocId(0, right.Side());
        intelright->SetConnectIndex(rightlocindex, index);
    }
    int sideorder = cleft.Order();
    flux_cmesh->SetDefaultOrder(sideorder);
    // create HDivBound on the sides of the elements
    TPZCompEl *wrap1, *wrap2;
    {
        intelright->Reference()->ResetReference();
        intelleft->LoadElementReference();
        intelleft->SetPreferredOrder(sideorder);
        TPZGeoElBC gbc(gleft, flux_trace_id);/// red flag!
        int64_t index;
        wrap1 = flux_cmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *flux_cmesh, index);
        if(cleft.Order() != sideorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap1);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelleft->Reference()->ResetReference();
        wrap1->Reference()->ResetReference();
    }
    // if the wrap of the large element was not created...
    if(equalright.size() == 0)
    {
        intelleft->Reference()->ResetReference();
        intelright->LoadElementReference();
        TPZConnect &cright = intelright->SideConnect(0,right.Side());
        int rightprevorder = cright.Order();
        intelright->SetPreferredOrder(cright.Order());
        TPZGeoElBC gbc(gright, flux_trace_id);
        int64_t index;
        wrap2 = flux_cmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *flux_cmesh, index);
        if(cright.Order() != rightprevorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap2);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelright->Reference()->ResetReference();
        wrap2->Reference()->ResetReference();
    }
    else
    {
        wrap2 = equalright[0].Element();
    }
    wrap1->LoadElementReference();
    wrap2->LoadElementReference();
    int64_t pressureindex;
    int pressureorder;
    {
        TPZGeoElBC gbc(gleft, lagrange_mult_id); /// red flag!
        pressureindex = gbc.CreatedElement()->Index();
        pressureorder = sideorder;
    }
    intelleft->LoadElementReference();
    intelright->LoadElementReference();
    return std::make_tuple(pressureindex, pressureorder);
    
}

void THybridzeDFN::CreateInterfaceElements(int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec) {

    TPZCompMesh * pressure_cmesh = mesh_vec[1];
    int dim = pressure_cmesh->Dimension();
    TPZVec<TPZCompEl *> celpressure(pressure_cmesh->NElements(), 0);
    for (int64_t el = 0; el < pressure_cmesh->NElements(); el++) {
        TPZCompEl *cel = pressure_cmesh->Element(el);
        if (cel && cel->Reference() && cel->Reference()->Dimension() == dim - 1) {
            celpressure[el] = cel;
        }
    }
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    
    TPZManVector<int64_t,3> left_mesh_indexes(1,0), right_mesh_indexes(1,0);
    for (auto cel : celpressure) {
        if (!cel) continue;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoEl *gel = cel->Reference();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        int count = 0;
        for (auto &celstackside : celstack) {
            if (celstackside.Reference().Element()->Dimension() == dim - 1) {
                TPZGeoElBC gbc(gelside, interface_id);

                TPZCompEl *right_cel = celstackside.Element();
                int n_connects = right_cel->NConnects();
                if (n_connects == 1) {
                    int64_t index;
                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstackside);
                    
                    TPZCompEl * left_cel = celside.Element();
                    TPZMultiphysicsElement * left_mp_cel = dynamic_cast<TPZMultiphysicsElement *>(left_cel);
                    if(left_mp_cel){
                        left_mesh_indexes[0]=1;
                    }
                    TPZMultiphysicsElement * right_mp_cel = dynamic_cast<TPZMultiphysicsElement *>(right_cel);
                    if(right_mp_cel){
                        right_mesh_indexes[0]=1;
                    }
                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
                    count++;
                }

            }
        }
        if (count == 1)
        {
            TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
            if(!clarge) DebugStop();
            TPZGeoElSide glarge = clarge.Reference();
            if (glarge.Element()->Dimension() == dim) {
                TPZGeoElSide neighbour = glarge.Neighbour();
                while(neighbour != glarge)
                {
                    if (neighbour.Element()->Dimension() < dim) {
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == glarge) DebugStop();
                glarge = neighbour;
            }
            clarge = glarge.Reference();
            if(!clarge) DebugStop();
            TPZGeoElBC gbc(gelside, interface_id);
            
            int64_t index;
            TPZMultiphysicsInterfaceElement *mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, clarge);
            TPZCompEl *right_cel = clarge.Element();
            
            TPZCompEl * left_cel = celside.Element();
            TPZMultiphysicsElement * left_mp_cel = dynamic_cast<TPZMultiphysicsElement *>(left_cel);
            if(left_mp_cel){
                left_mesh_indexes[0]=1;
            }
            TPZMultiphysicsElement * right_mp_cel = dynamic_cast<TPZMultiphysicsElement *>(right_cel);
            if(right_mp_cel){
                right_mesh_indexes[0]=1;
            }
            mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
            
            count++;
        }
        if (count != 2 && count != 0) {
            DebugStop();
        }
    }
    pressure_cmesh->InitializeBlock();
}
