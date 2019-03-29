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

#define NEW_V_Q

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

void THybridizeDFN::CreateInterfaceElements(int target_dim, int interface_id, TPZCompMesh *cmesh, TPZVec<TPZCompMesh *> & mesh_vec) {

    TPZCompMesh * pressure_cmesh = mesh_vec[1];
    TPZVec<TPZCompEl *> celpressure(pressure_cmesh->NElements(), 0);
    for (int64_t el = 0; el < pressure_cmesh->NElements(); el++) {
        TPZCompEl *cel = pressure_cmesh->Element(el);
        if (cel && cel->Reference() && cel->Reference()->Dimension() == target_dim - 1) {
            celpressure[el] = cel;
        }
    }
    cmesh->Reference()->ResetReference();
    for (auto &cel :cmesh->ElementVec()) {
        if(cel->Reference()->Dimension() == target_dim -1){
            cel->LoadElementReference();
        }
    }
    
    TPZManVector<int64_t,2> left_mesh_indexes(1,1), right_mesh_indexes(1,0);
    for (auto cel : celpressure) {
        if (!cel) continue;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoEl *gel = cel->Reference();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        int count = 0;
        for (auto &celstackside : celstack) {
            if (celstackside.Reference().Element()->Dimension() == target_dim - 1) {
                TPZGeoElBC gbc(gelside, interface_id);

                TPZCompEl * right_cel = celstackside.Element();
                int n_connects = right_cel->NConnects();
                if (n_connects == 1) {
                    int64_t index;
                    TPZMultiphysicsInterfaceElement * mp_interface_el = new TPZMultiphysicsInterfaceElement(*cmesh, gbc.CreatedElement(), index, celside, celstackside);
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
            if (glarge.Element()->Dimension() == target_dim) {
                TPZGeoElSide neighbour = glarge.Neighbour();
                while(neighbour != glarge)
                {
                    if (neighbour.Element()->Dimension() < target_dim) {
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
            mp_interface_el->SetLeftRightElementIndices(left_mesh_indexes,right_mesh_indexes);
            
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

void THybridizeDFN::ComputeAndInsertMaterials(int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id, int shift){
    
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
        flux_trace_id = maxMatId + 1 + shift;
        lagrange_id = maxMatId + 2 + shift;
        mp_nterface_id = maxMatId + 3 + shift;
    }
    
    /// Insert LM and flux_trace materials
    {
        TPZVec<STATE> sol(1,0);
        TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
        if (!flux_cmesh) {
            DebugStop();
        }
        if (!flux_cmesh->FindMaterial(flux_trace_id)) {
            auto flux_trace_mat = new TPZL2Projection(flux_trace_id, target_dim, n_state, sol);
            flux_trace_mat->SetScaleFactor(0.0);
            flux_cmesh->InsertMaterialObject(flux_trace_mat);
        }
        
        TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
        if (!pressure_cmesh) {
            DebugStop();
        }
        if (!pressure_cmesh->FindMaterial(lagrange_id)) {
            auto lagrange_mult_mat = new TPZL2Projection(lagrange_id, target_dim, n_state, sol);
            lagrange_mult_mat->SetScaleFactor(0.0);
            pressure_cmesh->InsertMaterialObject(lagrange_mult_mat);
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
    for (auto cel : flux_cmesh->ElementVec()) {
        if (cel->Reference()->Dimension() == target_dim and cel->NConnects() > 1) {
            cel->LoadElementReference();
        }
    }
    int64_t nel = flux_cmesh->NElements();
    std::list<std::tuple<int64_t, int> > pressures;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = flux_cmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel || (intel->Reference()->Dimension() != target_dim)) {
            continue;
        }
        
        if (intel->NConnects() == 1 ) {
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
    TPZVec<STATE> sol(1,0);
    if (!cmesh) {
        DebugStop();
    }
    if (!cmesh->FindMaterial(flux_trace_id)) {
        
        auto flux_trace_mat = new TPZL2Projection(flux_trace_id, target_dim, n_state, sol);
        flux_trace_mat->SetScaleFactor(0.0);
        cmesh->InsertMaterialObject(flux_trace_mat);
        
    }
    
    if (!cmesh->FindMaterial(lagrange_id)) {
        
        auto lagrange_mult_mat = new TPZL2Projection(lagrange_id, target_dim, n_state, sol);
        lagrange_mult_mat->SetScaleFactor(0.0);
        cmesh->InsertMaterialObject(lagrange_mult_mat);

    }
    
    if (!cmesh->FindMaterial(mp_nterface_id)) {
        auto lagrange_multiplier = new TPZLagrangeMultiplier(mp_nterface_id, target_dim - 1, n_state);
        lagrange_multiplier->SetMultiplier(1.0);
        cmesh->InsertMaterialObject(lagrange_multiplier);
    }
    
}

void THybridizeDFN::BuildMixedOperatorOnFractures(int p_order, int target_dim, TPZCompMesh * cmesh, int & flux_trace_id, int & lagrange_id, int & mp_nterface_id) {
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = mp_cmesh->MeshVector();
    
    /// insert material objects in pressure space
    {
        TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
        pressure_cmesh->Reference()->ResetReference();
        pressure_cmesh->LoadReferences();
        for(auto fracture : m_fracture_data){
            int fracture_id = fracture.m_id;
            auto fracture_material = new TPZMixedPoisson(fracture_id, target_dim-1);
            fracture_material->SetPermeability(fracture.m_kappa_normal);
            pressure_cmesh->InsertMaterialObject(fracture_material);
        }
        
    }
    
    /// Change lagrange multipliers identifiers to fractures identifiers
    {
        TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
        pressure_cmesh->Reference()->ResetReference();
        pressure_cmesh->LoadReferences();
        TPZGeoMesh  * gmesh = pressure_cmesh->Reference();
        
        {
            std::ofstream file_hybrid_mixed_q("b_mixed_cmesh_p.txt");
            pressure_cmesh->Print(file_hybrid_mixed_q);
        }
        
        
        for (auto gel : gmesh->ElementVec()) {
            
            if (!gel) {
                DebugStop();
            }
            
            if (gel->MaterialId() != lagrange_id) {
                continue;
            }
            
            TPZCompEl * cel = gel->Reference();
            if (!cel) {
                DebugStop();
            }
            
            if (gel->MaterialId() == lagrange_id) {
                int side = gel->NSides() - 1;
                TPZGeoElSide gel_side(gel,side);
                TPZStack<TPZGeoElSide> gelsides;
                gel_side.AllNeighbours(gelsides);
                
                for (auto iside : gelsides) {
                    TPZGeoEl * neigh_gel = iside.Element();
                    int gel_mat_id = neigh_gel->MaterialId();
                    std::set<int>::iterator it = m_fracture_ids.find(gel_mat_id);
                    bool fracture_detected_Q = it != m_fracture_ids.end();
                    if (fracture_detected_Q) {
                        cel->SetReference(neigh_gel->Index());
                    }
                    
                }
                
            }
            
        }
        
        {
            std::ofstream file_hybrid_mixed_q("a_mixed_cmesh_p.txt");
            pressure_cmesh->Print(file_hybrid_mixed_q);
        }
    }
    
    
    /// create flux space
    {
        TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
        flux_cmesh->Reference()->ResetReference();
//        flux_cmesh->LoadReferences();
        
//        {
//            std::ofstream file_hybrid_mixed_q("b_mixed_cmesh_q.txt");
//            flux_cmesh->Print(file_hybrid_mixed_q);
//        }
        
        std::set<int> fracture_set;
        for(auto fracture : m_fracture_data){
            int fracture_id = fracture.m_id;
            auto fracture_material = new TPZMixedPoisson(fracture_id, target_dim-1);
            fracture_material->SetPermeability(fracture.m_kappa_normal);
            flux_cmesh->InsertMaterialObject(fracture_material);
            fracture_set.insert(fracture_id);
        }
        
        
        
        
        /// Adding bc elements (This method must provide the correct material ids.)
        {
            /// Insert fractures intersections with bc
            {
                TPZGeoMesh * geometry = flux_cmesh->Reference();
                int dim = geometry->Dimension();
                std::map<std::pair<int,int>,std::pair<int,int>> surf_to_surf_side_indexes;
                std::vector<int> sides = {4,5,6,7};
                std::set<int> bc_indexes = {-1,-2,-3,-4,-5,-6};
                std::set<int> bc_indexes_1d = {-1000,-2000,-3000,-4000,-5000,-6000};
                for (auto gel : geometry->ElementVec()) {
                    
                    if (!gel) continue;
                    if (gel->Dimension() != dim - 1) continue;
                    int bc_mat_id = gel->MaterialId();
                    bool is_not_bc_memeber_Q = bc_indexes.find(bc_mat_id) ==  bc_indexes.end();
                    if (is_not_bc_memeber_Q) {
                        continue;
                    }
                    
                    
                    
                    for (auto side: sides) {
//                        std::cout << "Analysing side " << side << std::endl;
                        TPZStack<TPZGeoElSide> all_neigh;
                        TPZGeoElSide gelside(gel, side);
                        gelside.AllNeighbours(all_neigh);
                        std::set<int> surfaces;
                        surfaces.insert(gel->Index());

                        bool foundface = false;
                        bool foundbc = false;
                        int bcmatid = 0;
                        bool foundfrac = false;
                        for (auto gel_side : all_neigh) {
                            bool is_bc_member_Q = gel_side.Element()->MaterialId() == bc_mat_id;
                            if (is_bc_member_Q) {
                                foundface = true;
                                bcmatid = gel_side.Element()->MaterialId();
                            }
                            if(bc_indexes_1d.find(gel_side.Element()->MaterialId()) != bc_indexes_1d.end())
                            {
                                foundbc = true;

                            }
                            if(m_fracture_ids.find(gel_side.Element()->MaterialId()) != m_fracture_ids.end())
                            {
                                foundfrac = true;
                            }
//                            std::cout << "matid = " << gel_side.Element()->MaterialId() << std::endl;
                        }
                        if(foundface == false || foundbc == true || foundfrac == false) continue;

                        
//                        std::cout << "Creating boundary with matid " << bcmatid*1000 << std::endl;
                        for (auto gel_side : all_neigh) {
                            if(m_fracture_ids.find(gel_side.Element()->MaterialId()) != m_fracture_ids.end())
                            {
                                TPZGeoElBC gbc(gelside,bcmatid*1000);
                                break;
                            }
                        }
                    }
                }
            }
            
            fracture_set.insert(-1000);
            fracture_set.insert(-2000);
            fracture_set.insert(-3000);
            fracture_set.insert(-4000);
            fracture_set.insert(-5000);
            fracture_set.insert(-6000);
        }
        
        flux_cmesh->SetDimModel(target_dim);
        flux_cmesh->SetDefaultOrder(p_order);
        flux_cmesh->SetAllCreateFunctionsHDiv();
        flux_cmesh->AutoBuild(fracture_set);
        

        int sideorient = 1;
        for (auto gel : flux_cmesh->Reference()->ElementVec()) {
            if (!gel->Reference()) {
                continue;
            }
            
            int n_sides = gel->NSides();
            
            for (int is = 0 ; is < n_sides; is++) {
                if (gel->SideDimension(is) != target_dim-1) {
                    continue;
                }
                
                TPZInterpolatedElement * intel = dynamic_cast<TPZInterpolatedElement *>(gel->Reference());
            
                if (!intel) {
                    DebugStop();
                }
                
                intel->SetSideOrient(is, sideorient);
            }
        }
        
        flux_cmesh->InitializeBlock();
        
        
        {
            std::ofstream file_hybrid_mixed_q("a_mixed_cmesh_q.txt");
            flux_cmesh->Print(file_hybrid_mixed_q);
        }
        
        flux_cmesh->SetDimModel(target_dim); ///  for coherence
        
    }
}

void THybridizeDFN::InsertMaterialsForMixedOperatorOnFractures(int target_dim, TPZCompMesh * cmesh){
    
    
    std::vector<int> bc_indexes = {-2000,-4000,-5000,-600};
    std::vector<double> bc_vals = {0,0,2,1};
    std::vector<int> bc_type = {0,0,1,1};
    
    for(auto fracture : m_fracture_data){
        int fracture_id = fracture.m_id;
        if (!cmesh->FindMaterial(fracture_id)) {
            auto fracture_material = new TPZMatMixedPoisson3D(fracture_id, target_dim-1);
            fracture_material->SetPermeability(fracture.m_kappa_tangential);
            cmesh->InsertMaterialObject(fracture_material);
            
//            TPZFNMatrix<9,STATE> val1(1,1,0.0),val2(1,1,0.0);

//            for (int bc = 0; bc < bc_indexes.size(); bc++ ) {
//                val2(0,0)=bc_vals[bc];
//                int condType=bc_type[bc];
//                int bc_id = bc_indexes[bc];
//                TPZMaterial * face = fracture_material->CreateBC(fracture_material,bc_id,condType,val1,val2);
//                cmesh->InsertMaterialObject(face);
//            }
        }
    }
    

    
}

TPZCompMesh * THybridizeDFN::Hybridize(TPZCompMesh * cmesh, int target_dim){
    
    /// Step one
    int flux_trace_id, lagrange_id, mp_nterface_id;
    ComputeAndInsertMaterials(target_dim, cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
    
    /// Step two
    ApplyHibridizationOnInternalFaces(target_dim, cmesh, flux_trace_id, lagrange_id);
    
    
    {   /// 2D Fractures hybridization

        int target_dim = 2;

        /// Step one
        int flux_trace_2d_id, lagrange_2d_id, mp_nterface_2d_id;
        int mat_id_shift_2d = mp_nterface_id;

//        ComputeAndInsertMaterials(target_dim, cmesh, flux_trace_2d_id, lagrange_2d_id, mp_nterface_2d_id,mat_id_shift_2d);

//        /// Construction for mixed operator on fractures with provided target_dim
//        int p_order = 1;
//        BuildMixedOperatorOnFractures(p_order, target_dim, cmesh, flux_trace_id, lagrange_id, mp_nterface_id);

//        /// Step two
//        ApplyHibridizationOnInternalFaces(target_dim, cmesh, flux_trace_id, lagrange_id);

    }
    
    
    /// Step three reconstruct computational mesh
    TPZManVector<int,5> & active_approx_spaces = ExtractActiveApproxSpaces(cmesh);
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = GetMeshVector(cmesh);
    TPZCompMesh * dfn_hybrid_cmesh = DuplicateMultiphysicsCMeshMaterials(cmesh);
    CleanUpMultiphysicsCMesh(cmesh);
    InsertMaterialsForHibridization(target_dim, dfn_hybrid_cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
    
//    InsertMaterialsForMixedOperatorOnFractures(target_dim,dfn_hybrid_cmesh);
    
    {
        TPZMultiphysicsCompMesh  * mp_dfn_mixed_mesh_vec = dynamic_cast<TPZMultiphysicsCompMesh * >(dfn_hybrid_cmesh);
        if (!mp_dfn_mixed_mesh_vec) {
            DebugStop();
        }
        mp_dfn_mixed_mesh_vec->SetDimModel(target_dim);
        mp_dfn_mixed_mesh_vec->BuildMultiphysicsSpace(active_approx_spaces, dfn_mixed_mesh_vec);
    }
    dfn_mixed_mesh_vec[1]->SetDimModel(target_dim);
    CreateInterfaceElements(target_dim, mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    
//    dfn_mixed_mesh_vec[1]->SetDimModel(target_dim-1);
//    CreateInterfaceElements(mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    
    dfn_mixed_mesh_vec[1]->SetDimModel(target_dim);
    dfn_hybrid_cmesh->InitializeBlock();
   return dfn_hybrid_cmesh;
}

void THybridizeDFN::LoadReferencesByDimension(TPZCompMesh * flux_cmesh, int dim){

#ifdef PZDEBUG
    if (!flux_cmesh) {
        DebugStop();
    }
#endif
    
    TPZGeoMesh * geometry = flux_cmesh->Reference();
    geometry->ResetReference();
#ifdef PZDEBUG
    if (!geometry) {
        DebugStop();
    }
#endif
    
    for (auto cel: flux_cmesh->ElementVec()) {
        
        if (!cel) {
            continue;
        }
        
        TPZGeoEl * gel = cel->Reference();
        if (!gel) {
            continue;
        }
        
        bool check_vols = gel->Dimension() == dim && cel->NConnects() > 1; ///  Volumetric elements
        bool check_bc = gel->Dimension() == dim - 1 && cel->NConnects() == 1; ///  Boundary elements
        
        if (check_vols || check_bc) {
            cel->LoadElementReference();
        }
        
    }
    
}



TPZStack<TPZCompElSide> THybridizeDFN::AnalyzeSide(int target_dim, TPZGeoEl * gel, int side){
    
    if (!gel) {
        DebugStop();
    }

    
    TPZCompEl * cel = gel->Reference();
    if (!cel) {
        DebugStop();
    }
    
    /// Append the volumetric cel side on the stack
    TPZCompElSide cel_side(cel,side);
    
    TPZStack<TPZCompElSide> candidates;
    candidates.push_back(cel_side);
    
    /// Scattering for dim and dim-1 data
    TPZStack<TPZCompElSide> cel_side_vec_dim;
    TPZStack<TPZCompElSide> cel_side_vec_dim_m_1;
    
    TPZGeoElSide gel_side(gel,side);
    gel_side.EqualLevelCompElementList(candidates, 1, 1);
    
    
    for (auto cel_side: candidates) {
        TPZGeoEl * gel_side = cel_side.Element()->Reference();
        if (!gel_side) {
            DebugStop();
        }
        
        bool is_volumetric_Q = gel_side->Dimension() == target_dim;
        bool is_bc_Q = gel_side->Dimension() == target_dim - 1;
        
        if (is_volumetric_Q) {
            cel_side_vec_dim.push_back(cel_side);
        }else if (is_bc_Q){
            cel_side_vec_dim_m_1.push_back(cel_side);
        }
    }
    
    int n_dim_m_1_sides = cel_side_vec_dim_m_1.size();
    bool has_been_hybridize_Q = cel_side_vec_dim.size() == cel_side_vec_dim_m_1.size();
    if (has_been_hybridize_Q) {
        return candidates;
    }else if(n_dim_m_1_sides){
        DebugStop(); /// Nothing to say!
        return candidates;
    }
    
    return candidates;
    
}


std::pair<int,int> THybridizeDFN::HybridizeInSide(const TPZCompElSide &cel_side_l, const TPZCompElSide &cel_side_r, TPZCompMesh * flux_cmesh, int & flux_trace_id, int & lagrange_id){
    
    TPZGeoElSide gel_l(cel_side_l.Reference());
    TPZGeoElSide gel_r(cel_side_r.Reference());
    TPZInterpolatedElement *int_cel_l = dynamic_cast<TPZInterpolatedElement *> (cel_side_l.Element());
    TPZInterpolatedElement *int_cel_r = dynamic_cast<TPZInterpolatedElement *> (cel_side_r.Element());
    
    if (!int_cel_l || !int_cel_r) {
        DebugStop();
    }
    
    int_cel_l->SetSideOrient(cel_side_l.Side(), 1);
    int_cel_r->SetSideOrient(cel_side_r.Side(), 1);
    
    
    TPZConnect & connect_l = int_cel_l->SideConnect(0, cel_side_l.Side());
    
    TPZStack<TPZCompElSide> cel_r_neighbours;
    if (connect_l.HasDependency()) { /// Case for connect with dependency (non-conformal case) experiemntal version
        
        gel_r.EqualLevelCompElementList(cel_r_neighbours,1,0);
        
        // only one wrap element should exist
#ifdef PZDEBUG
        if(cel_r_neighbours.size() > 1)
        {
            DebugStop();
        }
        if(cel_r_neighbours.size()==1)
        {
            TPZGeoEl * gel = cel_r_neighbours[0].Element()->Reference();
            if(gel->Dimension() != flux_cmesh->Dimension()-1)
            {
                DebugStop();
            }
        }
#endif
        // reset the reference of the right flux trace element
        if(cel_r_neighbours.size()) {
            cel_r_neighbours[0].Element()->Reference()->ResetReference();
            connect_l.RemoveDepend();
        }
    }
    else
    {
        int64_t index = flux_cmesh->AllocateNewConnect(connect_l);
        TPZConnect &newcon = flux_cmesh->ConnectVec()[index];
        connect_l.DecrementElConnected();
        
        newcon.ResetElConnected();
        newcon.IncrementElConnected();
        newcon.SetSequenceNumber(flux_cmesh->NConnects() - 1);
        
        int local_connect_index = int_cel_r->SideConnectLocId(0, cel_side_r.Side());
        int_cel_r->SetConnectIndex(local_connect_index, index);
    }
    
    int connect_order_l = connect_l.Order();
    
    TPZCompEl * flux_trace_cel_l, * flux_trace_cel_r;
    
    { /// Create a computational element for the left side of flux trace
        int_cel_r->Reference()->ResetReference();
        int_cel_l->LoadElementReference();
        int_cel_l->SetPreferredOrder(connect_order_l);
        TPZGeoElBC gbc(gel_l, flux_trace_id);
        int64_t index;
        flux_trace_cel_l = flux_cmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *flux_cmesh, index);
        if(connect_l.Order() != connect_order_l)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (flux_trace_cel_l);
        int flux_trace_side = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(flux_trace_side, 1);
        int_cel_l->Reference()->ResetReference();
        flux_trace_cel_l->Reference()->ResetReference();
    }
    
    if(cel_r_neighbours.size() == 0)
    {   /// Create a computational element for the left side of flux trace
        int_cel_l->Reference()->ResetReference();
        int_cel_r->LoadElementReference();
        TPZConnect & connec_r = int_cel_r->SideConnect(0,cel_side_r.Side());
        int connect_order_r = connec_r.Order();
        int_cel_r->SetPreferredOrder(connect_order_r);
        TPZGeoElBC gbc(gel_r, flux_trace_id);
        int64_t index;
        flux_trace_cel_r = flux_cmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *flux_cmesh, index);
        if(connec_r.Order() != connect_order_r)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (flux_trace_cel_r);
        int flux_trace_side = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(flux_trace_side, 1);
        int_cel_r->Reference()->ResetReference();
        flux_trace_cel_r->Reference()->ResetReference();
    }
    else
    {
        flux_trace_cel_r = cel_r_neighbours[0].Element();
    }
    
    flux_trace_cel_l->LoadElementReference();
    flux_trace_cel_r->LoadElementReference();
    int64_t lambda_index;
    int lambda_order;
    {
        TPZGeoElBC gbc(gel_l, lagrange_id);
        lambda_index = gbc.CreatedElement()->Index();
        lambda_order = connect_order_l;
    }
    int_cel_l->LoadElementReference();
    int_cel_r->LoadElementReference();
    
    return std::make_pair(lambda_index, lambda_order);
}

int THybridizeDFN::CreateBCGeometricalElement(const TPZCompElSide & cel_side, TPZCompMesh * flux_cmesh,int & bc_impervious_id){
    int bc_gel_index;
    TPZGeoElSide gel_side(cel_side.Reference());
    TPZGeoElBC gbc(gel_side, bc_impervious_id);
    bc_gel_index = gbc.CreatedElement()->Index();
    return bc_gel_index;
}


void THybridizeDFN::ClassifyCompelSides(int target_dim, TPZCompMesh * flux_cmesh, TPZStack<std::pair<int, int>> & gel_index_and_order_lagrange_mult, TPZStack<TPZStack<TPZCompElSide>> & candidates_for_impervious_bc, int & flux_trace_id, int & lagrange_id){
    
    std::pair<int, int> gel_index_and_order;
    TPZGeoMesh * geometry = flux_cmesh->Reference();
    for (auto gel : geometry->ElementVec()) {
        
        if (!gel) {
            continue;
        }
        
        bool check = gel->Dimension() == target_dim;
        if (!check) {
            continue;
        }
        
        TPZCompEl * cel = gel->Reference();
        if (!cel) {
            DebugStop();
        }
        
        int n_sides = gel->NSides();
        for (int is = 0; is < n_sides; is++) {
            if(gel->SideDimension(is) != target_dim - 1){
                continue;
            }
            
            /// AnalyzeSide
            TPZStack<TPZCompElSide> candidates = AnalyzeSide(target_dim, gel, is);
            int n_candidates = candidates.size();
            
            bool needs_bc_Q = n_candidates == 1;
            if (needs_bc_Q) {
                candidates_for_impervious_bc.push_back(candidates);
                continue;
            }
            
            bool itself_and_bc_Q;
            if (n_candidates == 2) {
                TPZGeoEl * gel_bc = candidates[1].Element()->Reference(); ///  The first one is the element itself
                itself_and_bc_Q = gel_bc->Dimension() == target_dim - 1;
            }
            bool has_been_hybridize_Q = n_candidates == 0;
            if (has_been_hybridize_Q  || itself_and_bc_Q) {
                continue;
            }
            
            bool needs_hybridization_Q = n_candidates >= 2;
            if (needs_hybridization_Q) {
                if (n_candidates != 2) {
                    DebugStop();
                }
                gel_index_and_order = HybridizeInSide(candidates[0], candidates[1], flux_cmesh, flux_trace_id, lagrange_id);
                gel_index_and_order_lagrange_mult.push_back(gel_index_and_order);
                continue;
            }
            
        }
        
    }
    
    std::ofstream file("geometry_after_classify.txt");
    geometry->Print(file);
    
    std::ofstream file_geo("geometry_after_classify.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geometry, file_geo);
    
    flux_cmesh->InitializeBlock();
    flux_cmesh->ComputeNodElCon();
        
    
}

/// Construct a lagrange multiplier approximation space over the target dimension elements
TPZCompMesh * THybridizeDFN::Hybridize_II(TPZCompMesh * cmesh, int target_dim){
    
    TPZMultiphysicsCompMesh  * mp_cmesh = dynamic_cast<TPZMultiphysicsCompMesh * >(cmesh);
    if (!mp_cmesh) {
        DebugStop();
    }
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = mp_cmesh->MeshVector();
    
    TPZCompMesh * q_cmesh = dfn_mixed_mesh_vec[0];
    TPZCompMesh * p_cmesh = dfn_mixed_mesh_vec[1];
    
    if (!p_cmesh || !q_cmesh) {
        DebugStop();
    }
    int flux_trace_id, lagrange_id, mp_nterface_id;
    ComputeAndInsertMaterials(target_dim, mp_cmesh, flux_trace_id, lagrange_id, mp_nterface_id); /// split functionalities ...
    int bc_impervious_id = -1942;
    
    LoadReferencesByDimension(q_cmesh,target_dim);
    
    
    /// ClassifyCompelSides
    TPZStack<std::pair<int, int>> gel_index_and_order_lagrange_mult;
    TPZStack<TPZStack<TPZCompElSide>> candidates_for_impervious_bc;
    ClassifyCompelSides(target_dim, q_cmesh, gel_index_and_order_lagrange_mult, candidates_for_impervious_bc, flux_trace_id, lagrange_id);
    
    /// Create Impervious boundary elements
    TPZStack<int> bc_gel_indexes;
    for (auto candidates: candidates_for_impervious_bc) {
        int bc_index = CreateBCGeometricalElement(candidates[0], q_cmesh, bc_impervious_id);
        bc_gel_indexes.push_back(bc_index);
    }
    
    bool is_DFN_Hybridaze_Q = false;
    if (is_DFN_Hybridaze_Q) {
        DebugStop();
    }else{ /// Case for hybridization
        
        TPZGeoMesh * geometry = p_cmesh->Reference();
        geometry->ResetReference();
        
        p_cmesh->SetDimModel(target_dim-1);
        for (auto gel_index_and_order : gel_index_and_order_lagrange_mult) {
            
            int gel_index = gel_index_and_order.first;
            int cel_order = gel_index_and_order.second;
            
            TPZGeoEl * gel = geometry->Element(gel_index);
            int64_t cel_index;
            TPZCompEl * cel = p_cmesh->ApproxSpace().CreateCompEl(gel, *p_cmesh, cel_index);
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
            TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
            if (intel){
                intel->PRefine(cel_order);
            } else if (intelDisc) {
                intelDisc->SetDegree(cel_order);
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
        p_cmesh->InitializeBlock();
        p_cmesh->SetDimModel(geometry->Dimension());
        
    }
    
    /// Step three reconstruct computational mesh
    TPZManVector<int,5> & active_approx_spaces = ExtractActiveApproxSpaces(cmesh);
    TPZCompMesh * dfn_hybrid_cmesh = DuplicateMultiphysicsCMeshMaterials(cmesh);
    CleanUpMultiphysicsCMesh(cmesh);
    InsertMaterialsForHibridization(target_dim, dfn_hybrid_cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
    
    //    InsertMaterialsForMixedOperatorOnFractures(target_dim,dfn_hybrid_cmesh);
    
    {
        TPZMultiphysicsCompMesh  * mp_dfn_mixed_mesh_vec = dynamic_cast<TPZMultiphysicsCompMesh * >(dfn_hybrid_cmesh);
        if (!mp_dfn_mixed_mesh_vec) {
            DebugStop();
        }
        mp_dfn_mixed_mesh_vec->SetDimModel(target_dim);
        mp_dfn_mixed_mesh_vec->BuildMultiphysicsSpace(active_approx_spaces, dfn_mixed_mesh_vec);
    }
    CreateInterfaceElements(target_dim, mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    dfn_hybrid_cmesh->InitializeBlock();
    return dfn_hybrid_cmesh;
}
