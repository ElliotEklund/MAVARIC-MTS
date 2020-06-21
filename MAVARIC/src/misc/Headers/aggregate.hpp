#ifndef AGGREGATE
#define AGGREGATE

#include <string>
#include <map>
#include <iterator>
#include <fstream>
#include "mpi.h"
#include <vector>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

class aggregate{
    
public:
    aggregate(int my_id,int num_procs,int root_process);

/* Functions */
    void add_calc(std::string name, int num_cols, int num_rows);
    
    void collect(std::string name, int row, const vector<double> & v0,
                 const vector<double> & v,const double & sgnTheta);
    
    void collect(std::string name, int row, const double & v0,
                            const double & v, const double & sgnTheta);
    
    void print_collection(std::string name);
    
    void merge_collections(int root_process, int my_id, std::string root,
                           double dt, double ss, unsigned long long num_trajs);
    
    
    vector<vector<double> >  sum_data(vector<double> v_error,int num_samples,
                                      int num_errors);
    
    vector<double> get_samples(std::string name,int col,int num_errors);
    
    void write_errors(std::string name,int num_samples,int num_errors,double dt,
                      int stride,std::string output_dir);
    
private:
    
/* Data */
    int my_id, num_procs, root_proc; //mpi data
    std::map<std::string, matrix<double> * > myMap;
    
/* Functions */
    vector<double> average(vector<vector<double> > v, int num_samples,
                           int num_errors,vector<vector<double> > denom);
    
    vector<double> stdev(vector<vector<double> > v, vector<double> avgs,
                         vector<vector<double> > denom,int num_samples,
                         int num_errors);
    
};

#endif
