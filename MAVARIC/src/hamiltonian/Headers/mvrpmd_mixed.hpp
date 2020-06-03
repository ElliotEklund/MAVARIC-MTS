#ifndef mvrpmd_mixed_hpp
#define mvrpmd_mixed_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

class mvrpmd_mixed{
    
public:
    
    virtual double get_energy(const vector<double> &Q, const matrix<double> &x,
                              const matrix<double> &p) = 0;
    
    
};


#endif
