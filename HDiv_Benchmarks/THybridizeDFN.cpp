//
//  THybridizeDFN.cpp
//  HDiv
//
//  Created by Omar Durán on 3/23/19.
//

#include "THybridizeDFN.h"


/// Default constructor
THybridizeDFN::THybridizeDFN() : TPZHybridizeHDiv(){
 
    m_fracture_ids.clear();
    m_fracture_data.resize(0);
    m_geometry = NULL;
    m_geometry_dim = 0;
}

/// Default destructor
THybridizeDFN::~THybridizeDFN(){
    
}

std::tuple<int, int> THybridizeDFN::DissociateConnects(int flux_trace_id, int lagrange_mult_id, const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> & mesh_vec){
    
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

void THybridizeDFN::CreateInterfaceElements(int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec) {

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
                    if(left_cel->Reference()->MaterialId() == 100){
                        left_mesh_indexes[0]=1;
                    }
                    TPZMultiphysicsElement * right_mp_cel = dynamic_cast<TPZMultiphysicsElement *>(right_cel);
                    if(right_cel->Reference()->MaterialId() == 100){
                        right_mesh_indexes[0]=1;
                    }
//                    mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
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
            if(left_cel->Reference()->MaterialId() == 100){
                left_mesh_indexes[0]=1;
            }
            TPZMultiphysicsElement * right_mp_cel = dynamic_cast<TPZMultiphysicsElement *>(right_cel);
            if(right_cel->Reference()->MaterialId() == 100){
                right_mesh_indexes[0]=1;
            }
//            mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
            
            count++;
        }
        if (count != 2 && count != 0) {
            DebugStop();
        }
    }
    pressure_cmesh->InitializeBlock();
}


void THybridizeDFN::SetFractureData(TPZStack<TFracture> & fracture_data){
    m_fracture_data = fracture_data;
    /// Fill m_fracture_ids
    for (auto fracture: m_fracture_data) {
        m_fracture_ids.insert(fracture.m_id);
    }
}


void THybridizeDFN::SetDimension(int dimension){
    m_geometry_dim = dimension;
}

void THybridizeDFN::ComputeAndInsertMaterials(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id){
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = mp_cmesh->MeshVector();
    
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
        lagrange_id = maxMatId + 2;
        mp_nterface_id = maxMatId + 3;
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
        if (!pressure_cmesh->FindMaterial(lagrange_id)) {
            
            if (target_dim == 2) {
                auto lagrange_mult_mat = new TPZMat1dLin(lagrange_id);
                lagrange_mult_mat->SetMaterial(xk, xb, xc, xf);
                pressure_cmesh->InsertMaterialObject(lagrange_mult_mat);
            }else if(target_dim == 3){
                auto lagrange_mult_mat = new TPZMat2dLin(lagrange_id);
                lagrange_mult_mat->SetMaterial(xk, xc, xf);
                pressure_cmesh->InsertMaterialObject(lagrange_mult_mat);
            }
        }
        
    }
}


void THybridizeDFN::ApplyHibridizationOnInternalFaces(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id){
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = mp_cmesh->MeshVector();

    TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
    TPZGeoMesh * gmesh = flux_cmesh->Reference();

    gmesh->ResetReference();
    flux_cmesh->LoadReferences();
    int64_t nel = flux_cmesh->NElements();
    std::list<std::tuple<int64_t, int> > pressures;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = flux_cmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel || (intel->Reference()->Dimension() != target_dim)) {
            continue;
        }
        
        // loop over the side of dimension dim-1
        TPZGeoEl *gel = intel->Reference();
        for (int side = gel->NCornerNodes(); side < gel->NSides() - 1; side++) {
            if (gel->SideDimension(side) != target_dim - 1) {
                continue;
            }
            
            TPZCompElSide celside(intel, side);
            TPZCompElSide neighcomp = RightElement(intel, side); // Candidate to reimplement
            if (neighcomp) {
                pressures.push_back(DissociateConnects(flux_trace_id,lagrange_id,celside, neighcomp, dfn_mixed_mesh_vec));
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

TPZManVector<int,5> & THybridizeDFN::ExtractActiveApproxSpaces(TPZCompMesh * cmesh){
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    
    TPZManVector<int,5> & active_approx_spaces = mp_cmesh->GetActiveApproximationSpaces();
    return active_approx_spaces;
}

TPZManVector<TPZCompMesh *, 3> & THybridizeDFN::GetMeshVector(TPZCompMesh * cmesh){
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    return mp_cmesh->MeshVector();
}

TPZCompMesh * THybridizeDFN::DuplicateMultiphysicsCMeshMaterials(TPZCompMesh * cmesh){
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    
    TPZGeoMesh *gmesh = cmesh->Reference();
    TPZMultiphysicsCompMesh * dfn_hybrid_cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->CopyMaterials(*dfn_hybrid_cmesh);
    return dfn_hybrid_cmesh;
}

void THybridizeDFN::CleanUpMultiphysicsCMesh(TPZCompMesh * cmesh){
    cmesh->CleanUp();
    cmesh->SetReference(NULL);
}

void THybridizeDFN::InsertMaterialsForHibridization(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id){
    
    int n_state = 1;
    TPZFNMatrix<1, STATE> xk(n_state, n_state, 0.), xb(n_state, n_state, 0.), xc(n_state, n_state, 0.), xf(n_state, 1, 0.);
    
    if (!cmesh) {
        DebugStop();
    }
    if (!cmesh->FindMaterial(flux_trace_id)) {
        
        if (target_dim == 2) {
            auto flux_trace_mat = new TPZMat1dLin(flux_trace_id);
            flux_trace_mat->SetMaterial(xk, xb, xc, xf);
            cmesh->InsertMaterialObject(flux_trace_mat);
        }else if(target_dim == 3){
            auto flux_trace_mat = new TPZMat2dLin(flux_trace_id);
            flux_trace_mat->SetMaterial(xk, xc, xf);
            cmesh->InsertMaterialObject(flux_trace_mat);
        }
    }
    
    if (!cmesh->FindMaterial(lagrange_id)) {
        
        if (target_dim == 2) {
            auto lagrange_mult_mat = new TPZMat1dLin(lagrange_id);
            lagrange_mult_mat->SetMaterial(xk, xb, xc, xf);
            cmesh->InsertMaterialObject(lagrange_mult_mat);
        }else if(target_dim == 3){
            auto lagrange_mult_mat = new TPZMat2dLin(lagrange_id);
            lagrange_mult_mat->SetMaterial(xk, xc, xf);
            cmesh->InsertMaterialObject(lagrange_mult_mat);
        }
    }
    
    if (!cmesh->FindMaterial(mp_nterface_id)) {
        auto lagrange_multiplier = new TPZLagrangeMultiplier(mp_nterface_id, target_dim - 1, n_state);
        lagrange_multiplier->SetMultiplier(1.0);
        cmesh->InsertMaterialObject(lagrange_multiplier);
    }
    
}

TPZCompMesh * THybridizeDFN::Hybridize(TPZCompMesh * cmesh, int target_dim){
    
    /// Step one
    int flux_trace_id, lagrange_id, mp_nterface_id;
    ComputeAndInsertMaterials(target_dim, cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
    
    /// Step two
    ApplyHibridizationOnInternalFaces(target_dim, cmesh, flux_trace_id, lagrange_id);
    
    /// Step three reconstruct computational mesh
    TPZManVector<int,5> & active_approx_spaces = ExtractActiveApproxSpaces(cmesh);
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = GetMeshVector(cmesh);
    TPZCompMesh * dfn_hybrid_cmesh = DuplicateMultiphysicsCMeshMaterials(cmesh);
    CleanUpMultiphysicsCMesh(cmesh);
    InsertMaterialsForHibridization(target_dim, dfn_hybrid_cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
    
    {
        TPZMultiphysicsCompMesh  * mp_dfn_mixed_mesh_vec = dynamic_cast<TPZMultiphysicsCompMesh * >(dfn_hybrid_cmesh);
        if (!mp_dfn_mixed_mesh_vec) {
            DebugStop();
        }
        mp_dfn_mixed_mesh_vec->SetDimModel(target_dim);
        mp_dfn_mixed_mesh_vec->BuildMultiphysicsSpace(active_approx_spaces, dfn_mixed_mesh_vec);
    }
    CreateInterfaceElements(mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    dfn_hybrid_cmesh->InitializeBlock();
   return dfn_hybrid_cmesh;
}
