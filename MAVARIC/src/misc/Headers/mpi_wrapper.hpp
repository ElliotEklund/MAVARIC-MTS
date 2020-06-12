#ifndef mpi_wrapper_hpp
#define mpi_wrapper_hpp

#include "mpi.h"
#include <string>
#include <fstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

using namespace boost::numeric::ublas;


class mpi_wrapper{
    
public:
    mpi_wrapper(int num_procs, int my_id, int root_proc);
    
    void write_vector(vector<double> &v,std::string fileName,
                                   int nuc_beadsIN,int elec_beadsIN, int num_statesIN,
                                   double betaIN,unsigned long long num_trajs_totalIN,
                                   double decorrIN);
    
private:
    
/* Data */
    int num_procs, my_id, root_proc;
};

#endif
