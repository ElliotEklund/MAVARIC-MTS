#ifndef AGGREGATE
#define AGGREGATE

#include <string>
#include <map>
#include <iterator>
#include <fstream>
#include "mpi.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class aggregate{
    
public:
    aggregate();
    
    void add_calc(std::string name, int num_cols, int num_rows);
    
    void collect(std::string name, int row, const vector<double> & v0,
                 const vector<double> & v,const double & sgnTheta);
    
    void print_collection(std::string name);
    
    void merge_collections(int root_process, int my_id, std::string root,
                           double dt, double ss, unsigned long long num_trajs);
    
private:
    
/* Data */
    std::map<std::string, matrix<double> * > myMap;
    
};

#endif
