#ifndef rpmd_force_hpp
#define rpmd_force_hpp

#include "dSpring_dQ.hpp"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace boost::numeric::ublas;


class rpmd_force{
    
public:
    rpmd_force(int nuc_beads,double mass, double beta);
    
    /* Update dHdP to reflect the state of P1*/
    void update_dHdP(const vector<double> &P);
    
    /* Update dHdQ1 to reflect the state of Q1*/
    void update_dHdQ(const vector<double> &Q);
    
    /* Return dHdQ*/
    const vector<double> & get_dHdQ();
    
    /* Return dHdP*/
    const vector<double> & get_dHdP();
    
private:
    int nuc_beads;
    double c;
    dSpring_dQ dVspring_dQ;

    vector<double> dHdQ, dHdP; //derivative of Hamiltonian wrt Q,P
    vector<double> dVspring_dQ_vec;
    vector<double> dV0_dQ_vec;

    vector<double> dV0_dQ(const vector<double> &Q);
};

#endif
