//
//  THybridzeDFN.h
//  HDiv
//
//  Created by Omar Durán on 3/23/19.
//

#ifndef THybridzeDFN_h
#define THybridzeDFN_h

#include <stdio.h>
#include "TPZHybridizeHDiv.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzstack.h"
#include "pzintel.h"
#include "TFracture.h"

class THybridzeDFN : public TPZHybridizeHDiv {
    
private:
    
public:
    
    /// Default constructor
    THybridzeDFN();
    
    /// Default destructor
    ~THybridzeDFN();
    
    /// Construct a lagrange multiplier approximation space over the target dimension elements
    void Hybridize(TPZCompMesh * mp_cmesh, int target_dim, TPZStack<TFracture> & fracture_ids);
    
    /// Method that duplicate and dissociate the connect belonging to the right computational element side
    std::tuple<int, int> DissociateConnects(int flux_trace_id, int lagrange_mult_id, const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> & mesh_vec);
    
    /// Create interface multiphysics elements
    void CreateInterfaceElements(int interface_id, TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid);
    
};

#endif /* THybridzeDFN_h */
