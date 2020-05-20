#ifndef rpmd_vv_hpp
#define rpmd_vv_hpp

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "rpmd_force.hpp"

using namespace boost::numeric::ublas;

class rpmd_vv{
    
public:
    rpmd_vv(int nuc_beads,double mass, double beta,double dt);
    
    void step(vector<double> &Q,vector<double> &P);
    
private:
    
    double dt;
    double HALF_dt;
    vector<double> P_half;
    
    rpmd_force F;
};

#endif
