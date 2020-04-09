#ifndef dTheta_dBeta_hpp
#define dTheta_dBeta_hpp

#include "C_Matrix.h"
#include "M_Matrix.h"
#include "dM_Matrix_dBeta.hpp"

class dTheta_dBeta{
    
public:
    dTheta_dBeta(int num_beads,int num_states, double beta_num_beads,
                 C_Matrix &C_In, M_Matrix &M_In);
    
    
    /* Return d_Theta_dBeta; this function assumes d_Theta_dBeta has been updated*/
    double  get_dTheta_dBeta();
    
    /* Update d_Theta_dBeta based on the state of C and M*/
    void update_dTheta_dBeta();
    
    /* Update dGamma_dBeta based on the current state of
     C, M and dM_dBeta*/
    void update_dGamma_dBeta();
    
    /* Update f_chain based on the current state of C and M.*/
    void update_f_chain();
    
    /* Update b_chain based on the current state of C and M.*/
    void update_b_chain();
    
private:
    
    /* Private data. */
    int num_beads; //number of ring polymer beads
    int num_states; //number of electronic states
    
    C_Matrix *C;
    M_Matrix *M;
    dM_Matrix_dBeta dM_dBeta;
    
    matrix<std::complex<double> > dGamma_dBeta;
    double d_Theta_dBeta;
    
    vector<matrix<std::complex<double> > > f_chain; //forward chain
    vector<matrix<std::complex<double> > >b_chain; //backwards chain
    
    matrix<std::complex<double> > f_temp; //dummy matrix used in update_f_chain
    matrix<std::complex<double> > b_temp; //dummy matrix used in update_b_chain
    
    matrix<std::complex<double> > dummy1; //dummy matrix used in update_dGamma_dBeta
    matrix<std::complex<double> > dummy2; //dummy matrix used in update_dGamma_dBeta

};

#endif
