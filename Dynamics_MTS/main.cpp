#include <iostream>
#include <fstream>
#include "mpi.h"

#include "M_Matrix.h"
#include "C_Matrix.h"

#include "M_Matrix_MTS.hpp"
#include "dM_Matrix_MTS_dQ.hpp"
#include "dM_Matrix_dQ.hpp"

#include "Theta_MTS.hpp"
#include "dTheta_MTS_dQ.hpp"
#include "dTheta_MTS_dElec.hpp"

#include "Forces_MTS.hpp"
#include "RK4_MVRPMD.hpp"
#include "ABM_MVRPMD.hpp"

#include "Dynamics.hpp"
#include "PopulationEstimator.hpp"

int main(int argc, char ** argv) {
    
   int num_procs = 1; //number of processors program is distributed over
   int my_id = 0; //unique id of each processor
   int root_process = 0; //processor 0 is default root process
 
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
 
    MPI_Comm comm = MPI_COMM_WORLD;

    /* Variable declaration*/
    int nuc_beads, elec_beads;
    int num_states;
    double temp, beta;
    double beta_elec, beta_nuc;
    double mass;
    double dt;
    int num_trajs;
    double total_time;
    double tol = 0.001;
    int energy_stride = 100;

    bool run_PAC = false;
    bool run_PopAC = false;
    bool run_energ_conserv = false;

    /* Variable initialization */
    nuc_beads = 4;
    elec_beads = 4;
    num_states = 2;
    temp = 1.0;
    beta = 1.0/temp;
    beta_elec = beta/elec_beads;
    beta_nuc = beta/nuc_beads;
    mass = 1.0;
    dt = 0.001;
    total_time = 10;
    num_trajs = 1000;
    
    //run_energ_conserv = true;
    //run_PAC = true;
    run_PopAC = true;

    /* Setup Dynamics object for simulation. */
    Dynamics myDyn(num_procs,my_id,root_process,nuc_beads,
                   elec_beads,num_states,
                   mass,beta_nuc,beta_elec,num_trajs);
    
    myDyn.set_dt(dt);
    myDyn.set_total_time(total_time);


    clock_t start = clock();
    
    if(run_energ_conserv){
        myDyn.energ_conserv(tol,energy_stride);
    }

    if(run_PAC){
        myDyn.PAC();
    }

    if(run_PopAC){
        myDyn.PopAC();
    }
    
    clock_t end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    
    if (my_id == root_process) {
        std::cout << time_taken << std::endl;
    }

    MPI_Finalize();

    return 0;
}
