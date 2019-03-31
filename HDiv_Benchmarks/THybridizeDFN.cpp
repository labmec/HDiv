//
//  THybridizeDFN.cpp
//  HDiv
//
//  Created by Omar Durán on 3/23/19.
//

#include "THybridizeDFN.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr dfn_logger(Logger::getLogger("DFN"));
#endif


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
//        if (count != 2 && count != 0) {
//            DebugStop();
//        }
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
    
    
    /// Change lagrange multipliers identifiers to fractures identifiers
    {
        
        /// insert material objects in pressure space
        {
            TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
            pressure_cmesh->Reference()->ResetReference();
            pressure_cmesh->LoadReferences();
            for(auto fracture : m_fracture_data){
                int fracture_id = fracture.m_id;
                auto fracture_material = new TPZMixedDarcyFlow(fracture_id, target_dim-1);
                fracture_material->SetPermeability(fracture.m_kappa_normal);
                if (!pressure_cmesh->FindMaterial(fracture_id)) {
                    pressure_cmesh->InsertMaterialObject(fracture_material);
                }
            }
            
        }
        
        TPZCompMesh * pressure_cmesh = dfn_mixed_mesh_vec[1];
        pressure_cmesh->Reference()->ResetReference();
        pressure_cmesh->LoadReferences();
        TPZGeoMesh  * gmesh = pressure_cmesh->Reference();
        
        {
            std::ofstream file_hybrid_mixed_q("b_mixed_cmesh_p.txt");
            pressure_cmesh->Print(file_hybrid_mixed_q);
        }
        
        
        // switch the material id of pressure elements if they are neighbour of fracture elements
        for (auto gel : gmesh->ElementVec()) {
            
            if (!gel) {
                DebugStop();
            }
            
            if (gel->MaterialId() != lagrange_id) {
                continue;
            }
            
            TPZCompEl * cel = gel->Reference();
            if (!cel) {
                continue;
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
                    if (fracture_detected_Q && neigh_gel->Dimension() == gel->Dimension()) {
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
    
    
    /// create boundary elements for the 2d fractures
    if(target_dim == 2)
    {
        TPZCompMesh * flux_cmesh = dfn_mixed_mesh_vec[0];
        flux_cmesh->Reference()->ResetReference();
        
        // insert the material objects if needed
        std::set<int> fracture_set;
        for(auto fracture : m_fracture_data){
            int fracture_id = fracture.m_id;
            auto fracture_material = new TPZMixedDarcyFlow(fracture_id, target_dim-1);
            fracture_material->SetPermeability(fracture.m_kappa_normal);
            if (!flux_cmesh->FindMaterial(fracture_id)) {
                flux_cmesh->InsertMaterialObject(fracture_material);
            }
        }
        
        
        /// Adding bc elements (This method must provide the correct material ids.) for the fractures that have boundary neighbours
        if(target_dim == 2)
        {
            /// Insert fractures intersections with bc
            {
                TPZGeoMesh * geometry = flux_cmesh->Reference();
                int dim = geometry->Dimension();
                std::map<std::pair<int,int>,std::pair<int,int>> surf_to_surf_side_indexes;
                std::vector<int> sides;
                if (target_dim-1==1) {
                    sides = {4,5,6,7};
                }else if (target_dim-1==0){
                    sides = {0,1};
                }
                
                std::set<int> bc_indexes = {-1,-2,-3,-4,-5,-6};
                std::set<int> bc_indexes_1d = {-1000,-2000,-3000,-4000,-5000,-6000};
                for (auto gel : geometry->ElementVec()) {
                    
                    if (!gel) continue;
                    if (gel->Dimension() != dim - 1) continue;
                    int bc_mat_id = gel->MaterialId();
                    if(bc_indexes.find(bc_mat_id) ==  bc_indexes.end())
                    {
                        continue;
                    }
                    
                    // here we have a geometric element that corresponds to a 3d boundary (2 dimensional)
                    for (auto side: sides) {
                        TPZStack<TPZGeoElSide> all_neigh;
                        TPZGeoElSide gelside(gel, side);
                        // for each side of the boundary element, look for all neighbours
                        gelside.AllNeighbours(all_neigh);
                        std::set<int> surfaces;
                        surfaces.insert(gel->Index());

                        // foundface indicates whether there is another geometric element with the same material id
                        bool foundface = false;
                        // foundbc indicates we found an element with material id of fracture bcs
                        bool foundbc = false;
                        // bcmatid will be equal the boundary material of a neighbouring face
                        int bcmatid = 0;
                        // foundfrac indicates there is a neighbouring geometric element with a fracture material id
                        bool foundfrac = false;
                        for (auto gel_side : all_neigh) {
                            bool is_bc_member_Q = gel_side.Element()->MaterialId() == bc_mat_id;
                            if (is_bc_member_Q) {
                                foundface = true;
                                // this statement is redundant???
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
                        }
                        if(foundface == false || foundbc == true || foundfrac == false) continue;
                        for (auto gel_side : all_neigh) {
                            if(m_fracture_ids.find(gel_side.Element()->MaterialId()) != m_fracture_ids.end())
                            {
                                std::cout << "Inserted a boundary element of dimension " << gelside.Dimension() <<
                                " matid " << bcmatid*1000 << std::endl;
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
        
        
            LoadReferencesByDimension(flux_cmesh, target_dim);
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
        }
        
        {
            std::ofstream file_hybrid_mixed_q("a_mixed_cmesh_q.txt");
            flux_cmesh->Print(file_hybrid_mixed_q);
        }
        
        flux_cmesh->SetDimModel(target_dim); ///  for coherence
        
    }
}

void THybridizeDFN::InsertMaterialsForMixedOperatorOnFractures(int target_dim, TPZCompMesh * cmesh){
    
    
//    std::vector<int> bc_indexes = {-1000,-2000,-3000,-4000,-5000,-600};
//    std::vector<double> bc_vals = {0,0,2,1};
//    std::vector<int> bc_type = {0,0,1,1};
    
    for(auto fracture : m_fracture_data){
        int fracture_id = fracture.m_id;
        if (!cmesh->FindMaterial(fracture_id)) {
            auto fracture_material = new TPZMixedDarcyFlow(fracture_id, target_dim-1);
            fracture_material->SetPermeability(fracture.m_kappa_tangential);
            cmesh->InsertMaterialObject(fracture_material);
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

        /// Construction for mixed operator on fractures with provided target_dim
        int p_order = 1;
        BuildMixedOperatorOnFractures(p_order, target_dim, cmesh, flux_trace_id, lagrange_id, mp_nterface_id);

//        /// Step two
//        ApplyHibridizationOnInternalFaces(target_dim, cmesh, flux_trace_id, lagrange_id);

    }
    
    
    /// Step three reconstruct computational mesh
    TPZManVector<int,5> & active_approx_spaces = ExtractActiveApproxSpaces(cmesh);
    TPZManVector<TPZCompMesh *, 3> dfn_mixed_mesh_vec = GetMeshVector(cmesh);
    TPZCompMesh * dfn_hybrid_cmesh = DuplicateMultiphysicsCMeshMaterials(cmesh);
    CleanUpMultiphysicsCMesh(cmesh);
    InsertMaterialsForHibridization(target_dim, dfn_hybrid_cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
    
    InsertMaterialsForMixedOperatorOnFractures(target_dim,dfn_hybrid_cmesh);
    
    {
        TPZMultiphysicsCompMesh  * mp_dfn_mixed_mesh_vec = dynamic_cast<TPZMultiphysicsCompMesh * >(dfn_hybrid_cmesh);
        if (!mp_dfn_mixed_mesh_vec) {
            DebugStop();
        }
        mp_dfn_mixed_mesh_vec->SetDimModel(target_dim);
        mp_dfn_mixed_mesh_vec->BuildMultiphysicsSpace(active_approx_spaces, dfn_mixed_mesh_vec);
    }
    
    CreateInterfaceElements(target_dim, mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    
    CreateInterfaceElements(target_dim-1,mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    
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
            
#ifdef LOG4CXX
            if(dfn_logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << " LoadReferencesByDimension:: ";
                sout << " gel index =  " << gel->Index();
                sout << " matid = " << gel->MaterialId();
                sout << " connect indexes =  ";
                int n_connects = cel->NConnects();
                for (int i =0; i < n_connects; i++) {
                    sout << cel->ConnectIndex(i) << " ";
                }
                LOGPZ_DEBUG(dfn_logger,sout.str())
            }
#endif
            
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
    gel_side.EqualLevelCompElementList(candidates, 1, 0);
    
#ifdef LOG4CXX
    if(dfn_logger->isDebugEnabled())
    {
        int64_t n_candidates = candidates.size();
        std::stringstream sout;
        sout << "BEFORE FILTERING ******\n";
        sout << " gel index = " << gel->Index() << " side " << side;
        sout << " target_dim - 1 = " << target_dim - 1;
        TPZInterpolatedElement * intcel = dynamic_cast<TPZInterpolatedElement *>(candidates[0].Element());
        sout << " connect index = " << intcel->SideConnectIndex(0, side);
        sout << " Candidates (index, connectindex, side) = ";
        for (int i = 0 ; i < n_candidates; i++) {
            sout <<  candidates[i].Reference().Element()->Index() << " ";
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(candidates[i].Element());
            int side = candidates[i].Side();
            sout << intel->SideConnectIndex(0, side) << " ";
            sout <<  candidates[i].Side() << " * ";
        }
        LOGPZ_DEBUG(dfn_logger,sout.str())
    }
#endif

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
        candidates.Resize(0);
        return candidates;
    }else if(n_dim_m_1_sides != 0){
        DebugStop(); /// Nothing to say!
        return candidates;
    }
    
    return candidates;
    
}


std::pair<int,int> THybridizeDFN::HybridizeInSide(const TPZVec<TPZCompElSide> &candidates, TPZCompMesh * flux_cmesh, int & flux_trace_id, int & lagrange_id){
    
    if(candidates.size() < 2) DebugStop();
    TPZGeoElSide gel_first(candidates[0].Reference());
    TPZInterpolatedElement *int_cel_first = dynamic_cast<TPZInterpolatedElement *> (candidates[0].Element());
    
    int_cel_first->SetSideOrient(candidates[0].Side(), 1);
    
    for (auto it: candidates) {
        it.Element()->Reference()->ResetReference();
    }
    TPZStack<TPZCompEl *> element_to_load;
    size_t ncand = candidates.size();
    int connect_order = -1;
    int64_t connect_index = -1;
    {
        int side = candidates[0].Side();
        int_cel_first->SetSideOrient(side, 1);
        TPZConnect & connect = int_cel_first->SideConnect(0, side);
        connect_index = int_cel_first->SideConnectIndex(0, side);
        if(connect.HasDependency()) DebugStop();
        connect_order = connect.Order();
    }

    for (int i=0; i<ncand; i++) {
        
        TPZInterpolatedElement *int_cel = dynamic_cast<TPZInterpolatedElement *> (candidates[i].Element());
        int side = candidates[i].Side();
        int_cel->SetSideOrient(side, 1);
        TPZConnect & connect = int_cel->SideConnect(0, side);
        if(i>0)
        {
            int64_t global_index = flux_cmesh->AllocateNewConnect(connect);
            TPZConnect &newcon = flux_cmesh->ConnectVec()[global_index];
            connect.DecrementElConnected();
            
            newcon.ResetElConnected();
            newcon.IncrementElConnected();
            newcon.SetSequenceNumber(flux_cmesh->NConnects() - 1);
            
            int local_connect_index = int_cel->SideConnectLocId(0, side);
            int_cel->SetConnectIndex(local_connect_index, global_index);
            
            // create the HDiv bound element
            int_cel->LoadElementReference();

            int_cel->SetPreferredOrder(connect_order);
            connect_index = global_index;
        }
        TPZGeoElBC gbc(candidates[i].Reference(), flux_trace_id);
        int64_t compel_index;
        int_cel->LoadElementReference();
        TPZCompEl *flux_trace_cel = flux_cmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *flux_cmesh, compel_index);
        if(flux_trace_cel->Connect(0).Order() != connect_order)
        {
            DebugStop();
        }
        if(flux_trace_cel->ConnectIndex(0) != connect_index)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (flux_trace_cel);
        int flux_trace_side = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(flux_trace_side, 1);
        int_cel->Reference()->ResetReference();
        flux_trace_cel->Reference()->ResetReference();

        element_to_load.Push(flux_trace_cel);
        element_to_load.Push(int_cel);
    }
    
    
    int64_t lambda_index;
    int lambda_order;
    {
        TPZGeoElBC gbc(gel_first, lagrange_id);
        lambda_index = gbc.CreatedElement()->Index();
        lambda_order = connect_order;
    }
    for (auto cel:element_to_load) {
        cel->LoadElementReference();
    }
    
    return std::make_pair(lambda_index, lambda_order);
}

int THybridizeDFN::CreateBCGeometricalElement(const TPZCompElSide & cel_side, TPZCompMesh * flux_cmesh,int & bc_impervious_id){
    int bc_gel_index;
    TPZGeoElSide gel_side(cel_side.Reference());
    TPZGeoElBC gbc(gel_side, bc_impervious_id);
    bc_gel_index = gbc.CreatedElement()->Index();
    return bc_gel_index;
}


void THybridizeDFN::ClassifyCompelSides(int target_dim, TPZCompMesh * flux_cmesh, TPZStack<std::pair<int, int>> & gel_index_and_order_lagrange_mult, int impervious_bc_id, int & flux_trace_id, int & lagrange_id){
    
    flux_cmesh->SetDimModel(target_dim);
    flux_cmesh->SetAllCreateFunctionsHDiv();
    LoadReferencesByDimension(flux_cmesh, target_dim);
    if (!flux_cmesh->FindMaterial(impervious_bc_id)) {
        int n_state = 1;
        TPZVec<STATE> sol(1,0.);
        auto impervious_mat = new TPZL2Projection(impervious_bc_id, target_dim, n_state, sol);
        impervious_mat->SetScaleFactor(0.0);
        flux_cmesh->InsertMaterialObject(impervious_mat);
    }

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
            continue;
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
                TPZGeoElBC gbc(gel,is,impervious_bc_id);
                int64_t index;
                flux_cmesh->CreateCompEl(gbc.CreatedElement(), index);
#ifdef LOG4CXX
                if(dfn_logger->isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Created a impervious boundary compel for gel index " << gel->Index()
                    << " dim " << gel->Dimension() <<
                    " side " << is << std::endl;
                    LOGPZ_DEBUG(dfn_logger, sout.str())
                }
#endif
                continue;
            }
            
            // means there is only the element/side and a boundary condition. No hybridization needed
            bool itself_and_bc_Q = false;
            if (n_candidates == 2) {
                TPZGeoEl * gel_bc = candidates[1].Element()->Reference(); ///  The first one is the element itself
                itself_and_bc_Q = gel_bc->Dimension() == target_dim - 1;
            }
            // this variable is wrong because n_candidates is never 0 (it always includes the original element/side
            bool has_been_hybridize_Q = n_candidates == 0;
            if (has_been_hybridize_Q  || itself_and_bc_Q) {
                continue;
            }
            
#ifdef LOG4CXX
            if(dfn_logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << " gel index = " << gel->Index();
                sout << " target_dim - 1 = " << target_dim - 1;
                TPZInterpolatedElement * intcel = dynamic_cast<TPZInterpolatedElement *>(candidates[0].Element());
                sout << " connect index = " << intcel->SideConnectIndex(0, is);
                sout << " Candidates = ";
                for (int i = 0 ; i < n_candidates; i++) {
                    sout <<  candidates[i].Reference().Element()->Index() << " ";
                    sout <<  candidates[i].Reference().Side() << " * ";
                }
                LOGPZ_DEBUG(dfn_logger,sout.str())
            }
#endif
            
            bool needs_hybridization_Q = n_candidates >= 2;
            if (needs_hybridization_Q) {
                gel_index_and_order = HybridizeInSide(candidates, flux_cmesh, flux_trace_id, lagrange_id);
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
    int flux_trace_id = 0, lagrange_id = 0, mp_nterface_id = 0;
    ComputeAndInsertMaterials(target_dim, mp_cmesh, flux_trace_id, lagrange_id, mp_nterface_id); /// split functionalities ...
    int bc_impervious_id = -1942;
    
    LoadReferencesByDimension(q_cmesh,target_dim);
    
    
    /// ClassifyCompelSides
    TPZStack<std::pair<int, int>> gel_index_and_order_lagrange_mult;
    // hybridize the flux elements of target_dim
    // insert elements of bc_impervious_id if the element has no neighbour
    ClassifyCompelSides(target_dim, q_cmesh, gel_index_and_order_lagrange_mult, bc_impervious_id, flux_trace_id, lagrange_id);
    
    {
        std::ofstream out("Flux_hybrid3D.txt");
        q_cmesh->Print(out);
    }
    bool is_DFN_Hybridize_Q = true;
    bool is_1d_DFN_Hybridize_Q = true; /// tototototototo
    
    // create the 2D pressure interfaces
    if (is_DFN_Hybridize_Q) { /// Case for lagrange multiplier method on dfn
        
        {
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
        
        int p_order = 1;
        // switch the material id of the pressure elements if they are neighbour of 2d fracture elements
        // create boundary elements for the 2d fracture elements
        BuildMixedOperatorOnFractures(p_order, target_dim-1, cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
        
        // now we have a flux mesh that has been hybridized and the 2d pressure elements have been set to the fracture material if needed
        q_cmesh->ComputeNodElCon();
        {
            std::ofstream frac_file("frac_flux_with_boundary_cmesh.txt");
            q_cmesh->Print(frac_file);
        }
        
        {
            // load the 2d fracture elements
            LoadReferencesByDimension(q_cmesh,target_dim-1);
            
            /// ClassifyCompelSides
            TPZStack<std::pair<int, int>> gel_index_and_order_lagrange_mult;
            // hybridize the 2d fracture elements
            ClassifyCompelSides(target_dim-1, q_cmesh, gel_index_and_order_lagrange_mult, bc_impervious_id, flux_trace_id, lagrange_id);
            
            
            {
                std::ofstream out("flux_hybrid2d.txt");
                q_cmesh->Print(out);
            }
            
            if (is_1d_DFN_Hybridize_Q) {
                // create the 1 dimensional pressure space
                {
                    TPZGeoMesh * geometry = p_cmesh->Reference();
                    geometry->ResetReference();
                    
                    p_cmesh->SetDimModel(target_dim-2);
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
                            cel->Connect(i).SetLagrangeMultiplier(3);
                        }
                        gel->ResetReference();
                    }
                    p_cmesh->InitializeBlock();
                    p_cmesh->SetDimModel(geometry->Dimension());
                }
                
                int p_order = 1;
                // switch the material id of the pressure elements if they are neighbour of 1d fracture elements
                // create boundary elements for the 1d fracture elements
                BuildMixedOperatorOnFractures(p_order, target_dim-2, cmesh, flux_trace_id, lagrange_id, mp_nterface_id);
                
                q_cmesh->ComputeNodElCon();
                {
                    std::ofstream frac_file("frac_flux_with 1d_boundary.txt");
                    q_cmesh->Print(frac_file);
                }
                
                {
                    LoadReferencesByDimension(q_cmesh,target_dim-2);
                    
                    /// ClassifyCompelSides
                    // list of pressure elements that need to be created
                    TPZStack<std::pair<int, int>> gel_index_and_order_lagrange_mult;
                    ClassifyCompelSides(target_dim-2, q_cmesh, gel_index_and_order_lagrange_mult, bc_impervious_id, flux_trace_id, lagrange_id);
                    
                
                // create the 1 dimensional pressure space
                
                    TPZGeoMesh * geometry = p_cmesh->Reference();
                    geometry->ResetReference();
                    
                    p_cmesh->SetDimModel(target_dim-2);
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
                            cel->Connect(i).SetLagrangeMultiplier(4);
                        }
                        gel->ResetReference();
                    }
                    p_cmesh->InitializeBlock();
                    p_cmesh->SetDimModel(geometry->Dimension());
                }

            }else{
                TPZGeoMesh * geometry = p_cmesh->Reference();
                geometry->ResetReference();
                
                p_cmesh->SetDimModel(target_dim-2);
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
                        cel->Connect(i).SetLagrangeMultiplier(3);
                    }
                    gel->ResetReference();
                }
                p_cmesh->InitializeBlock();
                p_cmesh->SetDimModel(geometry->Dimension());
            }
            
        }
        
        
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
    
    if(is_DFN_Hybridize_Q){
        InsertMaterialsForMixedOperatorOnFractures(target_dim,dfn_hybrid_cmesh);
    }
    
    {
        TPZMultiphysicsCompMesh  * mp_dfn_mixed_mesh_vec = dynamic_cast<TPZMultiphysicsCompMesh * >(dfn_hybrid_cmesh);
        if (!mp_dfn_mixed_mesh_vec) {
            DebugStop();
        }
        mp_dfn_mixed_mesh_vec->SetDimModel(target_dim);
        mp_dfn_mixed_mesh_vec->BuildMultiphysicsSpace(active_approx_spaces, dfn_mixed_mesh_vec);
    }
    CreateInterfaceElements(target_dim, mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    
    if(is_DFN_Hybridize_Q){
        CreateInterfaceElements(target_dim-1, mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
        CreateInterfaceElements(target_dim-2, mp_nterface_id, dfn_hybrid_cmesh, dfn_mixed_mesh_vec);
    }
    
    dfn_hybrid_cmesh->InitializeBlock();
    return dfn_hybrid_cmesh;
}
