#ifndef init_PAC_hpp
#define init_PAC_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "theta_mixed.hpp"
#include "C_Matrix.h"
#include "M_Matrix.h"
#include "trajs_io.hpp"

#include <math.h>
#include <sstream>
#include <list>
#include <fstream>
#include <string>
#include "mpi.h"

class init_PAC{
    
public:
    init_PAC(int my_id,int num_procs, int root_proc,int num_trajs_total,
             int num_trajs_local);
    
/* Functions */

    /*
     Compute analyzes how the decorrelation length is affecting the convergence of
     zero time position auto-correlation.
     The output of compute is a file with four columns.
     Column 1: Number of samples used to generate data.
     Column 2: Sum of zero time position auto-correlation function
     Column 3: Sum of zero time sign(theta) term
     Column 4: The standard deviation of zero time position auto-correlation
     trajectory in which each point was monte carlo sampled with at least interval
     steps between every other trajectory.
     */
    void compute(std::string input_dir, std::string output_dir);
    
    /*
     Compute the zero time position auto-correlation for each trajectory.
     That is: Qcent(0)Qcent(0)sign(Theta), where Qcent is the centroid.
     theta_vec and qqTheta_vec are filled once the function as completed
     */
    void compute_vecs(std::string input_dir);
    
/* Mutators */
    void set_interval(int interval_IN);
    
    void set_vectors(vector<vector<double> > Q_IN,vector<matrix<double> > x_IN,
                     vector<matrix<double> > p_IN);
    
    void  set_system(int nuc_beadsIN,int elec_beadsIN,int num_statesIN,
                     double betaIN,double alphaIN);

private:
    
/* Data */
    int num_procs, my_id, root_proc; //mpi information
    int nuc_beads, elec_beads, num_states; //system information
    int num_trajs_total; //total number of trajectories
    int num_trajs_local; //number of trajectories per proc (num_trajs/num_proc)
    int interval; //number of trajectories per loop
    double alpha;
    double beta;
    
    /* Q(i)(j) = jth bead position of ith trajectory*/
    vector<vector<double> > Q;
    
    /* x(i)(j,k) = kth state of jth x MV bead of ith trajectory*/
    vector<matrix<double> > x;
    
    /* p(i)(j,k) = kth state of jth p MV bead of ith trajectory*/
    vector<matrix<double> > p;
    
    vector<double> theta_vec; //theta_vec[i] corresponds to sgn_theta the ith traj
    vector<double> qqTheta_vec; //qqTheta_vec[i] = Qcent(0)Qcent(0)sgnTheta(0) for ith traj
    vector<double> ones; //vector of 1.0s
    
/* Functions */
    double get_centroid(const vector<double> & Q_IN); //return centroid of Q_IN
    
};

#endif
