
#include "TPZTracerFlow.h"


TPZTracerFlow::TPZTracerFlow(){

}

/** @brief Constructor based on a material id */
TPZTracerFlow::TPZTracerFlow(int matid, int dimension){

}

/** @brief Constructor based on a TRMMultiphase object */
TPZTracerFlow::TPZTracerFlow(const TPZTracerFlow &mat){

}

/** @brief Default destructor */
TPZTracerFlow::~TPZTracerFlow(){

}

/** @brief Set the required data at each integration point */
void TPZTracerFlow::FillDataRequirements(TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

/** @brief Set the required data at each integration point */
void TPZTracerFlow::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec){
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}



/** @brief Print out the data associated with the material */
void TPZTracerFlow::Print(std::ostream &out){
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

/** @brief Returns the variable index associated with the name */
int TPZTracerFlow::VariableIndex(const std::string &name){

    if (!strcmp("Sw", name.c_str())) return 0;
    if (!strcmp("So", name.c_str())) return 1;

    return TPZMaterial::VariableIndex(name);
}

/** @brief Returns the number of variables associated with varindex */
int TPZTracerFlow::NSolutionVariables(int var){
    switch(var) {
        case 0:
            return 1; // Scalar
        case 1:
            return 1; // Scalar

    }
    return TPZMaterial::NSolutionVariables(var);
}



// Contribute Methods being used
/** @brief Returns the solution associated with the var index */
void TPZTracerFlow::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){

    int s_w    = 0;
    REAL sw = datavec[s_w].sol[0][0];

    Solout.Resize(this->NSolutionVariables(var));

    switch(var) {
        case 0:
        {
            Solout[0] = sw;
        }
            break;
        case 1:
        {
            Solout[0] = 1.0-sw;
        }
            break;
    }


}

void TPZTracerFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 3 ) {
        std::cout << " Erro. The size of the datavec is different from 3 \n";
        DebugStop();
    }
#endif

    // Setting the phis
    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
    TPZFMatrix<REAL>  &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL>  &phiP =  datavec[1].phi;
    TPZFMatrix<REAL>  &phiS =  datavec[2].phi;
    TPZFMatrix<REAL>  &dphiS =  datavec[2].dphix;


    TPZFMatrix<REAL> &axesS = datavec[2].axes;
    int phrQ = datavec[1].fVecShapeIndex.NElements();//phiQ.Rows();
    int phrP = phiP.Rows();
    int phrS = phiS.Rows();

    int sb_a    = 0;

    TPZFNMatrix<100,STATE> phi_ss       = datavec[sb_a].phi;
    REAL sa = datavec[sb_a].sol[0][0];

    int nphis_a     = phi_ss.Rows();
    int firsts_a    = 0;

    // Time
    STATE dt = 0.5;

    // Fluid parameters
    TPZManVector<STATE, 10> rho_w,rho_o,Bw_n,Bw;


    for (int is = 0; is < nphis_a; is++)
    {

        ef(is + firsts_a) += weight * (1.0/dt) * (phi_n*(sw_n) - phi*(sw)* phiS(is,0);

        for (int js = 0; js < nphis_a; js++)
        {
            ek(is + firsts_a, js + firsts_a) += weight * (1.0/dt) * phi_n*(1.0/Bw_n[0]) * phiS(js,0) * phiS(is,0);
        }

    }

}

void TPZTracerFlow::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef){
    
    
    
    
}


/**
 * Unique identifier for serialization purposes
 */
int TPZTracerFlow::ClassId() const{

}

/**
 * Save the element data to a stream
 */
void TPZTracerFlow::Write(TPZStream &buf, int withclassid){

}

/**
 * Read the element data from a stream
 */
void TPZTracerFlow::Read(TPZStream &buf, void *context){

//}
