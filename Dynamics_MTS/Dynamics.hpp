#ifndef Dynamics_hpp
#define Dynamics_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "M_Matrix.h"
#include "C_Matrix.h"

#include "M_Matrix_MTS.hpp"
#include "dM_Matrix_MTS_dQ.hpp"
#include "dM_Matrix_dQ.hpp"

#include "Theta_MTS.hpp"
#include "dTheta_MTS_dQ.hpp"
#include "dTheta_MTS_dElec.hpp"

#include "Forces_MTS.hpp"
#include "ABM_MVRPMD.hpp"

#include <math.h>
#include <fstream>
#include "mpi.h"

class Dynamics{
    
public:
    Dynamics(int num_procs,int my_id, int root_proc,int nuc_beads, 
             int elec_beads, int num_states, double mass,
             double beta_nuc_beads, double beta_elec_beads, int num_trajs);
    
    /* Run Molecular Dynamics Simulation*/
    void runSimulation();
    
    /* Set dt to dtIN*/
    void set_dt(double dtIN);

    /* Set total_time to total_timeIN*/
    void set_total_time(double total_timeIN);

private:

    /* Private Functions*/
    
    /* Write total position auto-correlation function to file.*/
    void print_QQt(const vector<double> &QQt,double sgnTheta_total);

    /* Given a vector of bead positions, Q, return the centroid*/
    double compute_centroid(const vector<double> &Q);

    /* Given a phase spave variable, X, read in values stored in
     file corresponding to var*/
    void load_var(vector<double> &X, std::string var, std::string root_path);

    /* X is a properly formated "array" containing the elements of X_local */
    void format_array(vector<vector<double> > &X, vector<double> &X_local);
  
    /* Overloaded version to format for matricies.
     X is a properly formated "array" containing the elements of X_local */
    void format_array(vector<matrix<double> > &X, vector<double> &X_local);

    /* Private Data */
    
    int num_procs; //number of processors
    int my_id; //unique processor id
    int root_proc; //id of root processor

    int elec_beads; //number of electronic beads
    int nuc_beads; //number of nuclear beads
    int num_states; //number of electronic states
    int num_trajs; //global number of trajectories
    int num_trajs_local; //trajectories per processor
    double mass; //nuclear mass
    double beta_nuc_beads; //beta/nuc_beads
    
    double dt; //time step [a.u]
    double total_time; //total molecular dynamics simulation time [a.u]
    int num_steps; //number of steps need to complete simulation
   
    /* Q(i)(j) = jth bead position of ith trajectory*/
    vector<vector<double> > Q;
    
    /* P(i)(j) = jth bead momentum of ith trajectory*/
    vector<vector<double> > P;
    
    /* x(i)(j,k) = kth state of jth x MV bead of ith trajectory*/
    vector<matrix<double> > x;
    
    /* p(i)(j,k) = kth state of jth p MV bead of ith trajectory*/
    vector<matrix<double> > p;

    vector<double> QQt; //auto correlation function
    
    /*  Objects needed to for force calculations */
    M_Matrix M;
    C_Matrix C;
    
    M_Matrix_MTS M_MTS;
    dM_Matrix_dQ dMdQ;
    dM_Matrix_MTS_dQ dM_MTS_dQ;
    
    Theta_MTS Theta;
    dTheta_MTS_dQ dThetadQ;
    dTheta_MTS_dElec dThetadElec;
    
    Forces_MTS F;
    
    /* Functions used for debugging */
    
    /* Write all position variables of a given trajectory at time=step to file. */
    void write_Q(std::ofstream &myStream, double step, vector<double> &Q);
    
    /* Write all momentum variables of a given trajectory at time=step to file. */
    void write_P(std::ofstream &myStream, double step, vector<double> &P);
    
    /* Write all x mapping variables of a given trajectory at time=step to file. */
    void write_x(std::ofstream &myStream, double step, matrix<double> &x);
    
    /* Write all p mapping variables of a given trajectory at time=step to file. */
    void write_p(std::ofstream &myStream, double step, matrix<double> &p);
    
};

#endif
