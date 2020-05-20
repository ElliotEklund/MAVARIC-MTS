#ifndef Forces_two_particles_hpp
#define Forces_two_particles_hpp

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "dSpring_dQ.hpp"
#include "dStateIndep_dQ.hpp"

using namespace boost::numeric::ublas;

class Forces_two_particles{
    
public:
    Forces_two_particles(int nuc_beads1, int nuc_beads2, double mass, double beta);
    
    
    /* Update dHdP1 to reflect the state of P1*/
    void update_dHdP1(const vector<double> &P1);
    
    /* Update dHdP2 to reflect the state of P2*/
    void update_dHdP2(const vector<double> &P2);
    
    /* Update dHdQ1 to reflect the state of Q1*/
    void update_dHdQ1(const vector<double> &Q1,const vector<double> &Q2);
    
    /* Update dHdQ2 to reflect the state of Q2*/
    void update_dHdQ2(const vector<double> &Q1,const vector<double> &Q2);
    

    void update_Forces(const vector<double> &Q1,const vector<double> &Q2,
                       const vector<double> &P1,const vector<double> &P2);
    
    /* Return dHdQ1*/
    const vector<double> & get_dHdQ1();
    
    /* Return dHdQ2*/
    const vector<double> & get_dHdQ2();
  
    /* Return dHdP1*/
    const vector<double> & get_dHdP1();
    
    /* Return dHdP2*/
    const vector<double> & get_dHdP2();
   
    /* Functions for debugging*/
    void print_dHdQ(int sub_system);
    void print_dHdP(int sub_system);
    
private:
    
    /* Private Data. */
    double nuc_beads1, nuc_beads2;
    double c1, c2; // 1.0/(mass * N_1(2)
    double d1, d2, d12;
    double one_nuc_beads1, one_nuc_beads2;
    
    mapped_matrix<double> W; //matrix used to calculate coupling energy

    double beta_nuc_beads; //beta/nuc_beads
    double ONE_beta_nuc_beads; //1.0/beta_nuc_beads
    double TWO_beta_nuc_beads; //2.0/beta_nuc_beads
    int nuc_beads; //number of nuclear beads

    vector<double> dHdQ1, dHdP1; //derivative of Hamiltonian wrt Q,P
    vector<double> dHdQ2, dHdP2; //derivative of Hamiltonian wrt Q,P

    vector<double> dVspring_dQ1_vec;
    vector<double> dVspring_dQ2_vec;

    vector<double> dV0_dQ1_vec;
    vector<double> dV0_dQ2_vec;
    
    vector<double> dV12_dQ1_vec;
    vector<double> dV12_dQ2_vec;

    /* Objects need for force calculation*/
    dSpring_dQ dVspring_dQ1;
    dSpring_dQ dVspring_dQ2;
    
    vector<double> dV0_dQ1(const vector<double> &Q);
    
    vector<double> dV0_dQ2(const vector<double> &Q);
    
    vector<double> dV12_dQ1(const vector<double> &Q1,const vector<double> &Q2);
    
    vector<double> dV12_dQ2(const vector<double> &Q1,const vector<double> &Q2);
};


#endif
