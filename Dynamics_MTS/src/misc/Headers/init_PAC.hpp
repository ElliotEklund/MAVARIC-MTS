#ifndef init_PAC_hpp
#define init_PAC_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "Theta_MTS.hpp"

#include <math.h>
#include <sstream>
#include <list>
#include <fstream>
#include <string>
#include "mpi.h"

class init_PAC{
    
public:
    init_PAC(int num_procs, int my_id, int root_proc, int nuc_beads, int elec_beads,
             int num_states, int num_trajs_total, int num_trajs_local, Theta_MTS Theta);
    
    void set_interval(int interval_IN);
    
    void set_vectors(vector<vector<double> > Q_IN,
                     vector<matrix<double> > x_IN, vector<matrix<double> > p_IN);
    
    void compute(std::string root);
    
    void compute_vecs();
    
private:
    
    int num_procs;
    int my_id;
    int root_proc;
    int nuc_beads;
    int elec_beads;
    int num_states;
    int num_trajs_total; //total number of trajectories
    int num_trajs_local; //number of trajectories per proc

    
    int interval; //number of trajectories per loop
    
    /* Q(i)(j) = jth bead position of ith trajectory*/
    vector<vector<double> > Q;
    
    /* x(i)(j,k) = kth state of jth x MV bead of ith trajectory*/
    vector<matrix<double> > x;
    
    /* p(i)(j,k) = kth state of jth p MV bead of ith trajectory*/
    vector<matrix<double> > p;
    
    vector<double> theta_vec; //theta_vec[i] corresponds to sgn_theta the ith traj
    vector<double> qqTheta_vec; //qqTheta_vec[i] = Qcent(0)Qcent(0)sgnTheta(0) for ith traj
    vector<double> ones; //vector of 1.0s

    Theta_MTS Theta;
    

    double get_centroid(const vector<double> & Q_IN);
    
};

#endif
