#ifndef GTERM_h
#define GTERM_h

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>
#include <numeric>

using namespace boost::numeric::ublas;

class GTerm{
    
public:
    
    GTerm(int num_beads, int num_states, double alpha);
    
    /* Update energy to reflect the state of x and p.
     x is a matrix electronic beads and states;
     x[i][j] is the ith bead and jth state.
     p is the same as x, but the conjugate variable*/
    void update_gTerm(const matrix<double> &x,const matrix<double> &p);
    
    /* Return energy; update_gTerm is called.
     x is a matrix electronic beads and states;
     x[i][j] is the ith bead and jth state.
     p is the same as x, but the conjugate variable*/
    double& get_gTerm(const matrix<double> &x,const matrix<double> &p);
    
    /* Return energy; assumes that energy has been updated.
     x is a matrix electronic beads and states;
     x[i][j] is the ith bead and jth state.
     p is the same as x, but the conjugate variable*/
    double& get_gTerm();

private:
    /* Private Structs. */
    struct square {
        double operator() (double Q) const;
    };
    
    /* Private data. */
    double energy; //energy of GTerm
    double alpha; //mapping variable prefactor
    double x_alpha, p_alpha; //x and p mapping variable prefactor
    matrix<double> x_squared; //each element of x squared;intermediate data
    matrix<double> p_squared; //each element of p squared;intermediate data

};

#endif
