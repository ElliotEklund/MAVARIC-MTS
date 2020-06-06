#ifndef theta_mixed_dBeta_hpp
#define theta_mixed_dBeta_hpp

#include "C_Matrix.h"
#include "M_Matrix.h"
#include "dM_Matrix_dBeta.hpp"
#include "functions.hpp"

class theta_mixed_dBeta{
    
public:
    
    theta_mixed_dBeta(int elec_beads, int num_states,
                      double beta_num_beads, C_Matrix &C_In, M_Matrix &M_In);
    
    /* Update and return the value of dtheta */
    double  get_dtheta();
    
    /* Update the value of dtheta. */
    void update_dtheta();
    
    /* Update dGamma_dBeta based on the current state of
     C, M and dM_dBeta*/
    void update_dgamma();
    
    /* Update f_chain based on the current state of C and M.*/
    void update_f_chain();
    
    /* Update b_chain based on the current state of C and M.*/
    void update_b_chain();
    
private:
    
/* Data */
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    
    matrix<std::complex<double> > dgamma; //derivative of gamma wrt beta
    double dtheta; //derivative of theta wrt beta using MTS method
    
    vector<matrix<std::complex<double> > > f_chain; //forward chain
    vector<matrix<std::complex<double> > >b_chain; //backwards chain
    
/* Objects */
    C_Matrix *C;
    M_Matrix *M;
    dM_Matrix_dBeta M_dbeta;
};

#endif
