#ifndef POP_ESTIMATORS_HPP
#define POP_ESTIMATORS_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>
#include <numeric>
#include <complex>

using namespace boost::numeric::ublas;

class pop_estimators{
    
public:
    pop_estimators(int elec_beads, int num_states,double alpha);
    
    /* Return the boltzmann population estimator*/
    vector<double> boltz(const matrix<std::complex<double> > gamma);
    
    /* Return the semi-classical population estimator*/
    vector<double> sc(const matrix<double> &x, const matrix<double> &p);
    
    /* Return the wigner population estimator*/
    vector<double> wigner(const matrix<double> &x, const matrix<double> &p);

private:
    
    /* Private data */
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    double sc_coef; //1.0/2.0 * elec_beads
    double wigner_coef;
    vector<double> gamma_nn; //element n correspons to gamma(n,n)
    double alpha;
    double x_alpha, p_alpha;
    
    matrix<double> x_sq; //element-wise square of x
    matrix<double> p_sq; //element-wise square of p
    matrix<double> xp_sq; // x_sq + p_sq - 0.5
    vector<double> x_sq_sum_bead; //sum over beads of x_sq
    vector<double> p_sq_sum_bead; //sum over beads of p_sq
    vector<double> x_sq_sum_state; //sum over states of x_sq
    vector<double> p_sq_sum_state; //sum over states of p_sq
    vector<double> exp_v; //vector taking exp(-x_sq_sum_state -p_sq_sum_state)
    vector<double> wig_pops; //final vector of wigner populaitons
    
    vector<double> one_elec_beads; //vector of ones one length elec_beads
    vector<double> one_num_states; //vector of ones one length num_states
    matrix<double> one_mat; //matrix of ones
};

#endif
