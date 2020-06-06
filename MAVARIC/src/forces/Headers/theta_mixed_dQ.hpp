#ifndef theta_mixed_dQ_hpp
#define theta_mixed_dQ_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <cmath>

#include "C_Matrix.h"
#include "M_Matrix.h"
#include "dM_Matrix_dQ.hpp"
#include "functions.hpp"

using namespace boost::numeric::ublas;

class theta_mixed_dQ{
    
public:
    
    theta_mixed_dQ(int num_states,int nuc_beads, int elec_beads, C_Matrix &C_In,
                  M_Matrix &M_IN,dM_Matrix_dQ &M_dQ_IN);
    
    
    /*
     Update theta_dQ_vec. This update uses the current state of C, M, and
     and updates M_dQ using Q.
     theta_dQ_vec will be updated when complete.
     */
    void update_theta_dQ(const vector<double> &Q);
    
/* Mutators */
    vector<double> & get_theta_dQ_vec();

private:
    
/* Data */
    int num_states; //number of electronic states
    int elec_beads; //number of electronic beads
    int nuc_beads; //number of nuclear beads
    int gamma; //nuc_beads/elec_beads
    vector<double> Q_trans; //Q transformed; result of W Q.
    mapped_matrix<double> W; //transformation matrix
    mapped_matrix<double> V; //transformation matrix

    /* f_chain[i] = C_{0} x M_{0} x C_{1} x M_{1} ... C_{i-1} x M_{i-1} x C_{i}
     f_chain is initialized to a zero vector when theta_mixed_dQ is created*/
    vector<matrix<std::complex<double> > > f_chain;
    
    /* b_chain[elec_beads] = I
       b_chain[i] = C_{elec_beads -(i+1)} x  M_{elec_beads -(i+1)} ...
                    x C_{elec_beads - 1} x M_{elec_beads - 1}
     b_chain is initialized to a zero vector when theta_mixed_dQ is created*/
    vector<matrix<std::complex<double> > > b_chain;
    
    matrix<std::complex<double> > fchain_temp; //used in update_f_chain()
    matrix<std::complex<double> > bchain_temp; //used in update_b_chain()
    matrix<std::complex<double> > temp1; //used in update_theta_dQ
    matrix<std::complex<double> > temp2; //used in update_theta_dQ
    
    /* theta_dQ_vec[i] = derivative of theta wrt Q_{i}
       theta_dQ_vec is initialized to a zero vector when theta_mixed_dQ is created*/
    vector<double> theta_dQ_vec;
    
/* Objects */
    C_Matrix *C; //C Matrix
    M_Matrix *M; //M Matrix
    dM_Matrix_dQ * M_dQ; //derivative of M Matrix wrt Q
    
/* Functions*/
    void update_f_chain(); //update f_chain
    void update_b_chain(); //update b_chain
};

#endif
