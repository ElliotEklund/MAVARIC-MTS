#ifndef dSpring_dQ_hpp
#define dSpring_dQ_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace boost::numeric::ublas;

class dSpring_dQ{
    
public:
    
    dSpring_dQ(int nuc_beads,double mass, double beta_nuc_beads);
    
    /* Return the derivative of the spring term w.r.t to Q*/
    const vector<double> & get_dSpring_dQ(const vector<double> &Q);
    
    const vector<double> & get_dSpring_dQ2(const vector<double> &Q);
    
private:
    
    /* Private Data. */
    double beta_nuc_beads; //1.0/temp
    int nuc_beads; //number of nuclear beads
    double mass; //system mass
    double coeff; //mass/(beta_nuc_beads^2)
    
    mapped_matrix<double> W; //matrix used to calculate coupling energy
    
    vector<double> force;
    
};

#endif
