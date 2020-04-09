#ifndef dTheta_MTS_dQ_hpp
#define dTheta_MTS_dQ_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <cmath>

#include "C_Matrix.h"
#include "M_Matrix_MTS.hpp"
#include "dM_Matrix_MTS_dQ.hpp"

using namespace boost::numeric::ublas;

class dTheta_MTS_dQ{
    
public:
    dTheta_MTS_dQ(int num_states,int nuc_beads, int elec_beads, C_Matrix &C_In,
                  M_Matrix_MTS &M_MTS_In,dM_Matrix_MTS_dQ &dM_MTS_dQ_In);
    
    
    /* Update dTheta_MTS_dQ_vec to reflect the state of C, M_MTS, and
     dM_MTS_dQ.*/
    void update_dTheta_MTS_dQ_vec();
    
    void update_dTheta_MTS_dQ_vec(const vector<double> &Q);

    const vector<double> & get_dThetaMTS_dQ_vec();

private:
    
    /* Private Functions*/
    void update_f_chain(); //update f_chain
    void update_b_chain(); //update b_chain

    
    /* Private Data*/
    int num_states; //number of electronic states
    int elec_beads; //number of electronic beads
    int nuc_beads; //number of nuclear beads
    int gamma; //nuc_beads/elec_beads

    matrix<std::complex<double> > fchain_temp; //used in f_chain
    matrix<std::complex<double> > bchain_temp; //used in b_chain
 
    matrix<std::complex<double> > temp1; //used in update_dTheta
    matrix<std::complex<double> > temp2; //used in update_dTheta


    C_Matrix *C;
    M_Matrix_MTS *M_MTS;
    dM_Matrix_MTS_dQ * dM_MTS_dQ;
    
    /* dTheta_MTS_dQ_vec[i] = dTheta_MTS/dQ_{i}*/
    vector<double> dTheta_MTS_dQ_vec;
    
    vector<matrix<std::complex<double> > > f_chain; //forward chain
    vector<matrix<std::complex<double> > > b_chain; //backwards chain
    
};


#endif
