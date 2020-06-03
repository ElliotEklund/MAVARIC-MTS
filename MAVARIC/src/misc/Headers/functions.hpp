#ifndef functions_hpp
#define functions_hpp

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

/* Returns the numerical sign of x*/
template <typename T>
double sign(T x){
    if (x >= 0){
        return 1.0;
    }
    else {
        return - 1.0;
    }
}

/* Compute the trace of a square matrix x*/
template <typename T>
T trace(const matrix<T> &x, int dim){
    
    T tr(0.0);
    for (int i=0; i<dim; i++) {
        tr += x(i,i);
    }
    return tr;
}

#endif
