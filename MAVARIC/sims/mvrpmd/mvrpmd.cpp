#include <iostream>
#include <fstream>
#include "mpi.h"

#include "MainHlpr.hpp"
#include "input_mvrpmd.hpp"
#include "equilib_mvrpmd.hpp"
#include "sampling_mvrpmd.hpp"
#include "dynamics_mvrpmd.hpp"
#include "trajs_io.hpp"
#include "aggregate.hpp"
//#include "Dynamics.hpp"

int main(int argc, char ** argv) {
    
    int num_procs = 1; //number of processors program is distributed over
    int my_id = 0; //unique id of each processor
    int root_process = 0; //processor 0 is default root process
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    
    MPI_Comm comm = MPI_COMM_WORLD;


    /* /////////////////////////////////////////////////// */
                        /* BEGIN PROCESS 1 */
    /* This process reads in all parameters stored in
     InputFiles and distributes them to their appropriate
     variables.*/

    /* Vectors used to store parameters from InputFiles */
    std::vector<double> sys_parameters;
    std::vector<double> elec_parameters;
    std::vector<double> MC_parameters;
    std::vector<double> Samp_parameters;
    std::vector<double> Dyn_parameters;

    input_mvrpmd myInput;

    std::string root = "/Users/ellioteklund/Desktop/MAVARIC_v2.0/MAVARIC/sims/mvrpmd/";
//    //std::string root = "/home/fs01/ece52/MAVARIC-MTS/MAVARIC/sims/mvrpmd/";
//
    int abort = myInput.input_file_handler(root,sys_parameters,elec_parameters,
                                           MC_parameters,Samp_parameters,Dyn_parameters);

    if (abort == -1){return -1;}

    /* Physical parameters.*/
    double temp, mass;
    int nuc_beads,num_states, elec_beads;
    double beta, beta_nuc, beta_elec;
    double alpha;

    /* Monte Carlo parameters. */
    unsigned long long num_steps, esti_rate;
    double nuc_ss, x_ss, p_ss;
    bool writePSV, readPSV;
    bool readData, writeData;
    bool runMC;

    /* Sampling parameters.*/
    bool runSamp, saveTrajs;
    unsigned long long num_trajs;
    unsigned long long decor_len;

    /* Dynamics Variables */
    bool runDyn;
    double dt;
    double total_time;
    double tol;
    int energy_stride;
    bool run_init_PAC;
    double interval;

    bool run_PAC = false;
    bool run_PopAC = false;
    bool run_energ_conserv = false;

    /* From MonteCarlo */
    runMC = MC_parameters[0];
    num_steps = MC_parameters[1];
    esti_rate = MC_parameters[2];
    writePSV = MC_parameters[3];
    writeData = MC_parameters[4];
    readPSV = MC_parameters[5];
    readData = MC_parameters[6];

    /* From SystemParameters*/
    mass = sys_parameters[0];
    nuc_beads = sys_parameters[1];
    temp = sys_parameters[2];
    nuc_ss = sys_parameters[3];

    /* From ElecParameters*/
    num_states = elec_parameters[0];
    elec_beads = elec_parameters[1];
    x_ss = elec_parameters[2];
    p_ss = elec_parameters[3];
    alpha = elec_parameters[4];

    /* From Sampling */
    runSamp = Samp_parameters[0];
    num_trajs = Samp_parameters[1];
    decor_len = Samp_parameters[2];
    saveTrajs = Samp_parameters[3];

    /* From Dynamics */
    runDyn = Dyn_parameters[0];
    total_time = Dyn_parameters[1];
    dt = Dyn_parameters[2];
    double junk = Dyn_parameters[3];
    run_PAC = Dyn_parameters[4];
    run_energ_conserv = Dyn_parameters[5];
    tol = Dyn_parameters[6];
    run_PopAC = Dyn_parameters[7];
    run_init_PAC = Dyn_parameters[8];
    interval = Dyn_parameters[9];

    /* Set derived variables */
    beta = 1.0/temp;
    beta_nuc = beta/nuc_beads;
    beta_elec = beta/elec_beads;

                        /* END PROCESS 1 */
    /* /////////////////////////////////////////////////////////*/


    /* /////////////////////////////////////////////////////////// */
                          /* BEGIN PROCESS 2 */
      /* This process runs the Monte Carlo simulation if requested.*/

    if (runMC) {
        if (my_id == root_process) {
            std::cout << std::endl << std::endl;
            std::cout << "Begin Monte Carlo Simulation" << std::endl;
            std::cout << std::endl;
        }

        equilib_mvrpmd equilibrator(my_id,root_process,num_procs,root);
        equilibrator.initialize_system(nuc_beads,elec_beads,num_states,
                                       mass,beta,alpha);

        equilibrator.initialize_files(writePSV,readPSV,writeData,readData);

        clock_t start = clock();
        
        equilibrator.run(nuc_ss,x_ss,p_ss,num_steps,esti_rate);

        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "\t Monte Carlo simulation time: " << time_taken <<
            std::endl << std::endl;
            
            std::cout << "End Monte Carlo Simulation" << std::endl;
            std::cout << std::endl;
        }
    }

                           /* END PROCESS 2 */
    /* /////////////////////////////////////////////////////////// */


    /* /////////////////////////////////////////////////////////// */
                          /* BEGIN PROCESS 3 */
      /* This process runs Sampling  if requested.*/

    if(runSamp){
        if (my_id == root_process) {
            std::cout << "Begin Sampling Simulation" << std::endl;
            std::cout << std::endl;
        }
        
        sampling_mvrpmd sampler(my_id,root_process,num_procs);
        sampler.initialize_system(nuc_beads,elec_beads,num_states,mass,beta,alpha);
        sampler.initialize_files(readPSV,saveTrajs,root);

        clock_t start = clock();

        sampler.run(nuc_ss,x_ss,p_ss,num_trajs,decor_len);

        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "\t Sampling simulation time: " << time_taken <<
            std::endl << std::endl;
            
            std::cout << "End Sampling Simulation" << std::endl;
            std::cout << std::endl << std::endl;
        }
    }

                            /* END PROCESS 3 */
  /* /////////////////////////////////////////////////////////// */


    /* /////////////////////////////////////////////////////////// */
                          /* BEGIN PROCESS 4 */
    /* This process runs the Dynamics simulation if requested.*/
    

    if(runDyn){
        if (my_id == root_process) {
            std::cout << "Begin Dynamics Simulation" << std::endl;
            std::cout << std::endl << std::endl;
        }

        /* Setup Dynamics object for simulation. */
        dynamics_mvrpmd dyn(my_id,num_procs,root_process);
        
        dyn.set_system(nuc_beads,elec_beads,num_states,mass,beta,beta_nuc,
                       beta_elec,alpha);
        dyn.set_time(dt,total_time);
        dyn.set_trajs(num_trajs,root);

        clock_t start = clock();

        if(run_energ_conserv){
            energy_stride = 100;
            dyn.energy_conserve(tol,energy_stride,root + "Output/Trajectories/",
                               root + "Output/");
        }
        if(run_PopAC){
            bool pac = true;
            bool bp = true;
            bool sp = true;
            bool wp = false;

            int pac_stride = 100;
            int bp_stride = 100;
            int sp_stride = 100;
            int wp_stride = 100;
            int num_samples = 4;
            int num_errors = 10;

            dyn.compute_ac(pac,pac_stride,bp,bp_stride,sp,sp_stride,wp,wp_stride,
                           root + "Output/Trajectories/",root + "Output/",
                           num_samples,num_errors);
        }

        if(run_init_PAC){
          dyn.iPAC(interval,root+"Output/Trajectories/",
                                root+"Output/");
        }

        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "Dynamics simulation time:" << time_taken << std::endl;
            std::cout << "End Dynamics Simulation" << std::endl;
        }
    }
    
                            /* END PROCESS 4 */
    /* /////////////////////////////////////////////////////////// */
    
    MPI_Finalize();
    
    return 0;
}
