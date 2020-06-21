#ifndef functions_hpp
#define functions_hpp

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>

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

/* Return true if x is nan or inf, otherwise return false*/
template<typename T>
bool is_NaN(const T &x){
    if (x != x) {
        return true;
    }
    else{
        return isinf(x);
    }
}

/* Return true if any element of v is nan or inf. Otherwise, return false*/
template <typename T>
bool contains_NaN(const vector<T> &v){
    
    int s = v.size();
    for(int i=0; i<s; i++){
        if(is_NaN(v(i))){
            return true;
        }
    }
    return false;
}
/* Return true if any element of v is nan or inf. Otherwise, return false*/
template <typename T>
bool contains_NaN(const matrix<T> &v){
    
    int num_row = v.size1();
    int num_col = v.size2();
    
    for (int row=0; row<num_row; row++) {
        for (int col=0; col<num_col; col++) {
            if (is_NaN(v(row,col))) {
                return true;
            }
        }
    }
    return false;
}

template<typename T>
bool contains_NaN(const vector<vector<T> > &v){
    
    int s = v.size();
    
    for (int i=0; i<s; i++) {
        if (contains_NaN(v(i))) {
            return true;
        }
    }
    return false;
}

#endif
