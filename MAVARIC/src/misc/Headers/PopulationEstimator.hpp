#ifndef POPULATIONESTIMATOR_hpp
#define POPULATIONESTIMATOR_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>

using namespace boost::numeric::ublas;

class PopulationEstimator{
    
public:
    PopulationEstimator(int elec_beads, int num_states);
    
    void update_g(const matrix<double> &x, const matrix<double> &p);
    
    void update_populations(const matrix<double> &x,const matrix<double> &p);
    
    double get_pop(const int state);
    
    vector<double> get_pop();
    
    void write_populations(matrix<double> final_pops, double dt, int data_steps, int rate, std::string root);

private:
    
    /* Private data */
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    
    double coeff;
    
    matrix<double> xx; //element-wise square of x
    matrix<double> pp; //element-wise square of p
    matrix<double> g; // g = xx + pp
    vector<double> g_sum; //g summed across rows
    vector<double> ones; //vector of ones
    
    vector<double> populations; //populations[i] is the population of state i
    
};

#endif
