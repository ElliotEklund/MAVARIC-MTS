#ifndef dTheta_MTS_dElec_hpp
#define dTheta_MTS_dElec_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "C_Matrix.h"
#include "M_Matrix_MTS.hpp"
#include "dCdelec.hpp"
#include "functions.hpp"

using namespace boost::numeric::ublas;

class dTheta_MTS_dElec{
    
public:
    dTheta_MTS_dElec(int num_states, int elec_beads, C_Matrix &C_In,
                     M_Matrix_MTS &M_MTS_In);
    
    void update_dTheta_MTS_dElec(const matrix<std::complex<double> > & x,
                                 const matrix<std::complex<double> > & p);

    void update_dTheta_MTS_dx_vec(const matrix<std::complex<double> > & x,
                                  const matrix<std::complex<double> > & p);
    
    void update_dTheta_MTS_dp_vec(const matrix<std::complex<double> > & x,
                                  const matrix<std::complex<double> > & p);

    const matrix<double> & get_dThetaMTS_dx_vec();
    
    const matrix<double> & get_dThetaMTS_dp_vec();

private:
    
    /* Private Functions*/
    void update_f_chain(); //update f_chain
    void update_b_chain(); //update b_chain
    
    /* Private Data*/
    int num_states; //number of electronic states
    int elec_beads; //number of electronic beads
    
    C_Matrix *C;
    dCdelec dC;
    M_Matrix_MTS *M_MTS;
    
    /* dTheat_MTS_dx_vec[i,j] = dTheta/x[i,j]; derivative of theta
     w.r.t to the ith bead and jth state of x*/
    matrix<double> dTheta_MTS_dx_vec;
    
    /* dTheat_MTS_dp_vec[i,j] = dTheta/p[i,j]; derivative of theta
     w.r.t to the ith bead and jth state of p*/
    matrix<double> dTheta_MTS_dp_vec;
    
    vector<matrix<std::complex<double> > > f_chain; //forward chain
    vector<matrix<std::complex<double> > > b_chain; //backwards chain

    matrix<std::complex<double> > fchain_temp; //used in update_f_chain
    matrix<std::complex<double> > bchain_temp; //used in update_b_chain
    
    matrix<std::complex<double> > xtemp1; //used in update calls
    matrix<std::complex<double> > xtemp2; //used in update calls

    matrix<std::complex<double> > ptemp1; //used in update calls
    matrix<std::complex<double> > ptemp2; //used in update calls
};

#endif
