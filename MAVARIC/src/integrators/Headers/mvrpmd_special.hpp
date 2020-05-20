#ifndef mvrpmd_special_hpp
#define mvrpmd_special_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <math.h>

#include "Forces_MTS.hpp"

using namespace boost::numeric::ublas;


class mvrpmd_special{
    
public:
    mvrpmd_special(Forces_MTS *F_In,int nuc_beads, int elec_beads,
                   int num_states,double dt);
    
    void step(vector<double> &Q, vector<double> &P,
              matrix<double> &x, matrix<double> &p);
    
private:
    double dt, dt_half, dt_tenth;
    vector<double> P_half;
    matrix<double> x_next, p_next;
    matrix<double> x_half, p_half;
    
    Forces_MTS *F;

    
    
};


#endif
