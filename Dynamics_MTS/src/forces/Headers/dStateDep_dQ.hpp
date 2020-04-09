#ifndef dStateDep_dQ_hpp
#define dStateDep_dQ_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

class dStateDep_dQ{
    
public:
    dStateDep_dQ(int nuc_beads, int num_states);
    
    /* Update dVdQ to reflect the state of Q*/
    void update_dVdQ(const vector<double> &Q);
    
    /* Return a reference to dVdQ*/
    const matrix<double> & get_dVdQ();
    
    /* Return a reference to dVcoup_dQ_vec*/
    const vector<matrix<double> > &get_dVcoup_dQ_vec(const vector<double> &Q);

    
private:
    
    /* Private Data. */
    int nuc_beads; //number of nuclear beads
    int num_states; //number of electronic states
    
    /* dVdQ(i,k) = dV_{kk}(Q_{i})/dQ_{i}  */
    matrix<double> dVdQ;
    
    /* dVcoup_dQ_vec[i](j,k) = dV_{j,k}(Q_{i})/dQ_{i}*/
    vector<matrix<double> > dVcoup_dQ_vec;
    
    
    /* Private functions. */
    
    /* Return the derivative of V11 w.r.t Q*/
    double dV11_dQ(const double &Q);
    
    /* Return the derivative of V22 w.r.t Q*/
    double dV22_dQ(const double &Q);

    /* Return the derivative of V22 w.r.t Q*/
    double dV33_dQ(const double &Q);

    /* Return the derivative of V12 w.r.t Q*/
    double dV12_dQ(const double &Q);
    
    double dV13_dQ(const double &Q);
    
    double dV23_dQ(const double &Q);
};

#endif
