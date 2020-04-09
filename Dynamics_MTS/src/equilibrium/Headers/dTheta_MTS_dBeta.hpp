#ifndef dTheta_MTS_dBeta_hpp
#define dTheta_MTS_dBeta_hpp

#include "C_Matrix.h"
#include "M_Matrix_MTS.hpp"
#include "M_Matrix.h"
#include "dM_Matrix_MTS_dBeta.hpp"

class dTheta_MTS_dBeta{
    
public:
    
    dTheta_MTS_dBeta(int nuc_beads, int elec_beads, int num_states, double beta_num_beads,
                     C_Matrix &C_In, M_Matrix &M_In, M_Matrix_MTS &M_MTS_In);
    
    double  get_dTheta_MTS_dBeta();
    
    void update_dTheta_MTS_dBeta();
    
    /* Update dGamma_dBeta based on the current state of
     C, M and dM_dBeta*/
    void update_dGamma_dBeta();
    
    /* Update f_chain based on the current state of C and M.*/
    void update_f_chain();
    
    /* Update b_chain based on the current state of C and M.*/
    void update_b_chain();
    
    
private:
    
    
    /* Private data. */
    int nuc_beads; //number of enuclear beads
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    
    C_Matrix *C;
    M_Matrix_MTS *M_MTS;
    dM_Matrix_MTS_dBeta dM_MTS_dBeta;
    
    matrix<std::complex<double> > dGamma_dBeta;
    double d_Theta_MTS_dBeta;
    
    vector<matrix<std::complex<double> > > f_chain; //forward chain
    vector<matrix<std::complex<double> > >b_chain; //backwards chain
    
};



#endif
