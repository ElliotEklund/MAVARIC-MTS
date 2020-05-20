#ifndef dM_Matrix_dQ_hpp
#define dM_Matrix_dQ_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <complex>

#include "M_Matrix.h"
#include "dStateDep_dQ.hpp"

using namespace boost::numeric::ublas;

class dM_Matrix_dQ{
    
public:
    dM_Matrix_dQ(int nuc_beads,int num_states, double betaN, M_Matrix &M_In);
    
    /* Update dM_dQ_vec to reflect the state of Q.
     It is assummed that M has already been updated*/
    void update_dM_dQ_vec(const vector<double> &Q);
    
    /* Return dMdQ_vec*/
    const vector<matrix<std::complex<double> > > & get_dMdQ_vec();
    
    /* Return dM/dQ_{\alpha} */
    const matrix<std::complex<double> > & get_dMdQ_alpha(int alpha);

    
private:
    
    /* Private Data. */
    int nuc_beads; //number of nuclear ring polymer beads
    int num_states; //number of electronic states
    
    /* beta/(number of beads) this will depend on MTS or no MTS. */
    double betaN;
    
    vector<matrix<std::complex<double> > > dM_dQ_vec; //dM_dQ_vec[i] = dM_dQ_{i}
    matrix<double> dVdQ_mat; //dVdQ_mat(i,j) = dV_{jj}/dQ_{i}
    
    /* dVcouple_dQ_vec[i](j,k) = dV_{j,k}/dQ_{i} */
    vector<matrix<double> > dVcouple_dQ_vec;
    
    matrix<std::complex<double> > dM_dQ_alpha; // dM/dQ_{alpha}
    
    M_Matrix *M;
    
    dStateDep_dQ dVdQ;
    
    
    /* Private Functions */
    
    /* Update dM_dQ_alpha to reflect the state of Q. */
    void update_dM_dQ_alpha(const double &Q, int alpha);

};

#endif
