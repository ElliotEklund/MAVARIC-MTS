/* Mapping Variable Forces Template*/

#ifndef mv_forces_temp_hpp
#define mv_forces_temp_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

class mv_forces_temp{
public:

    virtual void update_Forces(const vector<double> &Q,const vector<double> &P,
                               const matrix<double> &x,const matrix<double> &p) = 0;
    
    virtual const vector<double> & get_dHdQ() = 0;
    
    virtual const vector<double> & get_dHdP() = 0;
     
    virtual const matrix<double> & get_dHdx() = 0;

    virtual const matrix<double> & get_dHdp() = 0;
    
    virtual double get_sign(const vector<double> &Q, const matrix<double> &x,
                            const matrix<double> &p) = 0;

};

#endif
