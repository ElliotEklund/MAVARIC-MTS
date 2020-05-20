#ifndef two_particle_Dynamics_HPP
#define two_particle_Dynamics_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <fstream>
#include <string>
#include "mpi.h"

#include "VV_two_particle.hpp"

using namespace boost::numeric::ublas;

class two_particle_Dynamics{
    
    
public:
    two_particle_Dynamics(int my_id, int root_proc, int num_procs);
    
    void run();
    
    void initialize_system(int num_beads1IN, int num_beads2IN,double massIN, double betaIN);
    
    void initialize_dynamics(double dtIN, double total_tIN,int num_trajsIN,std::string rootFolderIN);
    
private:
    int my_id;
    int root_proc;
    int num_procs;
    
    int num_beads1, num_beads2;
    double beta, mass;
    double dt, total_t;
    int num_trajs;
    std::string rootFolder;
    
    void load_var(vector<double> &X, std::string var, std::string root_path);
    
    vector<vector <double> > get_trajs(std::string file,std::string root, int num_trajs,int num_beads);
    
    void format_array(int num_beads, int num_trajs_local, vector<vector<double> > &X,
                      vector<double> &X_local);
    
    double get_centroid(const vector<double> &x, int num_beads);
    
    void write_correlation(std::string file, vector<double> &QQ);


};

#endif
