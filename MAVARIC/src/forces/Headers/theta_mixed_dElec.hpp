#ifndef theta_mixed_dElec_hpp
#define theta_mixed_dElec_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "C_Matrix.h"
#include "M_Matrix.h"
#include "dCdelec.hpp"
#include "functions.hpp"

using namespace boost::numeric::ublas;

class theta_mixed_dElec{
    
public:
    theta_mixed_dElec(int num_states, int elec_beads,double alpha,
                      C_Matrix &C_In, M_Matrix &M_IN);
    
    void update_theta_dElec(const matrix<std::complex<double> > & x,
                                 const matrix<std::complex<double> > & p);

    void update_theta_dx_vec(const matrix<std::complex<double> > & x,
                                  const matrix<std::complex<double> > & p);
    
    void update_theta_dp_vec(const matrix<std::complex<double> > & x,
                                  const matrix<std::complex<double> > & p);

    const matrix<double> & get_theta_dx_vec();
    
    const matrix<double> & get_theta_dp_vec();

private:
    
/* Data */
    int num_states; //number of electronic states
    int elec_beads; //number of electronic beads
    
    /* dTheat_MTS_dx_vec[i,j] = dTheta/x[i,j]; derivative of theta
     w.r.t to the ith bead and jth state of x*/
    matrix<double> theta_dx_vec;
    
    /* dTheat_MTS_dp_vec[i,j] = dTheta/p[i,j]; derivative of theta
     w.r.t to the ith bead and jth state of p*/
    matrix<double> theta_dp_vec;
    
    vector<matrix<std::complex<double> > > f_chain; //forward chain
    vector<matrix<std::complex<double> > > b_chain; //backwards chain

    matrix<std::complex<double> > fchain_temp; //used in update_f_chain
    matrix<std::complex<double> > bchain_temp; //used in update_b_chain
    
    matrix<std::complex<double> > xtemp1; //used in update calls
    matrix<std::complex<double> > xtemp2; //used in update calls

    matrix<std::complex<double> > ptemp1; //used in update calls
    matrix<std::complex<double> > ptemp2; //used in update calls
    
/* Objects */
    C_Matrix *C; //C Matrix
    dCdelec dC; //derivative of C w.r.t x or p
    M_Matrix *M; //M Matrix
    
/* Functions */
    void update_f_chain(); //update f_chain
    void update_b_chain(); //update b_chain
};

#endif
