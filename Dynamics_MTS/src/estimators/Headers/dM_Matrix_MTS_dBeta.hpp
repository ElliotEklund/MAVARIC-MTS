#ifndef dM_Matrix_MTS_dBeta_hpp
#define dM_Matrix_MTS_dBeta_hpp

#include "dM_Matrix_dBeta.hpp"
#include "M_Matrix.h"

class dM_Matrix_MTS_dBeta{
    
public:
    dM_Matrix_MTS_dBeta(int nuc_beads, int elec_beads, int num_states,
                        double beta_nuc_beads,M_Matrix &M_In);
    
    void update_dM_MTS_dBeta_vec();
    
    matrix<std::complex<double> >& get_dM_MTS_dBeta_alpha(int alpha);
    
private:
    
    vector<matrix<std::complex<double> > >  dM_MTS_dBeta_vec;
    dM_Matrix_dBeta dM_dBeta;
    
    int elec_beads;
    int num_states;
    int gamma;
    
};

#endif
