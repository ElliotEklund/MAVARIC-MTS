#ifndef dM_Matrix_MTS_dQ_hpp
#define dM_Matrix_MTS_dQ_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <complex>

#include "dM_Matrix_dQ.hpp"

using namespace boost::numeric::ublas;

class dM_Matrix_MTS_dQ{
    
public:
    
    dM_Matrix_MTS_dQ(int nuc_beads, int elec_beads, int num_states, dM_Matrix_dQ &dMdQ_In);
        
    matrix<std::complex<double> > get_dM_MTS_dQ_alpha(int alpha);
    
    void update_dM_MTS_dQ_vec(const vector<double> &Q);

    
private:
    
    int nuc_beads; //number of nuclear beads
    int elec_beads; //number of electronic beads
    int gamma; //nuc_beads/elec_beads
    double ONE_gamma; //1.0/gamma
    
    dM_Matrix_dQ *dMdQ;
    vector<matrix<std::complex<double> > > dM_MTS_dQ_vec;
    
};

#endif
