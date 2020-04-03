#ifndef M_MATRIX_H
#define M_MATRIX_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>

#include "StateDepPots.h"

using namespace boost::numeric::ublas;

class M_Matrix{
    
public:
    M_Matrix(int num_states,int num_beads,double beta_num_beads);
    
    /* M_vec is updated to reflect the state of Q. */
    void update_M_vec(const vector<double> &Q);
    
    /* Return M_vec[alpha] */
    matrix<std::complex<double> >& get_M_alpha(int alpha);
    
    /* Return V_mat */
    const matrix<double>& get_V_mat();
    
    /* Return coupling matrix corresponding to alphath bead*/
    const matrix<double>& get_V_couple_alpha(int alpha);

private:
    
    /* Private Structs. */
    
    /* exp_mat is only used for evaluating the following function
     element-wise over a matrix: f(x) = exp(-beta_num_beads)*/
    struct exp_mat {
        double operator() (double x) const;
        double beta_num_beads; //beta/num_beads
        void set_beta_num_beads(double x); //set beta_num_beads = x
    };
    
    /* Private functions.*/
    
    /* Update M_alpha to reflect the state corresponding
     to the given bead. The potential energies defined
     in V are used. */
    void update_M_alpha(const double &Q, int bead);
    
    /* Update V_couple_vec to reflect state of Q*/
    void update_V_couple_vec(const vector<double> &Q);
    
    /* Private data. */
    int num_states; //number of states
    int num_beads; //number of beads
    double beta_num_beads; //beta/num_beads
    vector<double> exp_row; //temp vector; exp of exp_v_mat
    matrix<double> V_couple_temp; //temp matrix 
    matrix<double> M_inter; //intermediate term

    /* Intermediate M-Matrix corresponding to alphath bead.*/
    matrix<std::complex<double> > M_alpha;
    
    /* matrix with dim (nuc_beads x num_states)
     V_mat(i,j) =  V_{jj}(Q[i]) where,
     V_{jj} is the diagonal state dependent potential of state j and,
     Q[i] is the i-th bead's position*/
    matrix<double> V_mat;
    
    matrix<double> exp_V_mat; //V_mat with myExp applied to each element
    
    vector<matrix<std::complex<double> > > M_vec; //element alpha corresponds to M_{alpha}
    
    StateDepPots V; //state dependent potential energy object
    exp_mat myExp; //exp_mat object used to evaluate element-wise matrix exp
    
    vector<matrix<double> > V_couple_vec; //vector of coupling matricies
};

#endif
