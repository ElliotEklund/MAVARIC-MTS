#ifndef dStateIndep_dQ_hpp
#define dStateIndep_dQ_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

class dStateIndep_dQ{
    
public:
    
    dStateIndep_dQ(int nuc_beads,double mass);
    
    const vector<double> & get_dStateIndep_dQ(const vector<double> &Q);
    
private:
    
    /* Private Data. */
    
    vector<double> force;
    
    double mass; //system mass
    int nuc_beads; //number of nuclear beads
    double coef; //coefficient for potential

};

#endif
