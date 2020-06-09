#include <iostream>
#include <fstream>
#include "mpi.h"
#include <vector>

#include "csrpmd_sampling.hpp"

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
//    std::vector<double> sys_parameters;
//    std::vector<double> elec_parameters;
//    std::vector<double> MC_parameters;
//    std::vector<double> Samp_parameters;
//    std::vector<double> Dyn_parameters;

    //input_mvrpmd myInput;

    std::string root_dir = "/Users/ellioteklund/Desktop/MAVARIC_v2.0/MAVARIC/sims/csrpmd/";
//    //std::string root = "/home/fs01/ece52/MAVARIC-MTS/Dynamics_MTS/";
//
//    int abort = myInput.input_file_handler(root,sys_parameters,elec_parameters,
//                                           MC_parameters,Samp_parameters,Dyn_parameters);

//    if (abort == -1){return -1;}

    /* Physical parameters.*/
    double temp, mass;
    int nuc_beads,num_states, elec_beads;
    double beta, beta_nuc, beta_elec;
    
    int num_trajs;
    double decorr, rate, dt;
    
    bool run_samp = true;
    
    temp = 1.0;
    mass = 1.0;
    beta = 1.0/temp;
    nuc_beads = 6;
    elec_beads = 1;
    num_states = 2;
    
    num_trajs = 1;
    decorr = 200;
    rate = 10;
    dt = 0.01;
    
    
//
//    /* Monte Carlo parameters. */
//    unsigned long long num_steps, esti_rate;
//    double nuc_ss, elec_ss;
//    bool writePSV, readPSV;
//    bool readData, writeData;
//    bool runMC;
//
//    /* Sampling parameters.*/
//    bool runSamp, saveTrajs;
//    int num_trajs;
//    unsigned long long decor_len;
//
//    /* Dynamics Variables */
//    bool runDyn;
//    double dt;
//    double total_time;
//    double tol;
//    int energy_stride;
//    bool run_init_PAC;
//    double interval;
//
//    bool run_PAC = false;
//    bool run_PopAC = false;
//    bool run_energ_conserv = false;
//
//    /* From MonteCarlo */
//    runMC = MC_parameters[0];
//    num_steps = MC_parameters[1];
//    esti_rate = MC_parameters[2];
//    writePSV = MC_parameters[3];
//    writeData = MC_parameters[4];
//    readPSV = MC_parameters[5];
//    readData = MC_parameters[6];
//
//    /* From SystemParameters*/
//    mass = sys_parameters[0];
//    nuc_beads = sys_parameters[1];
//    temp = sys_parameters[2];
//    nuc_ss = sys_parameters[3];
//
//    /* From ElecParameters*/
//    num_states = elec_parameters[0];
//    elec_beads = elec_parameters[1];
//    elec_ss = elec_parameters[2];
//
//    /* From Sampling */
//    runSamp = Samp_parameters[0];
//    num_trajs = Samp_parameters[1];
//    decor_len = Samp_parameters[2];
//    saveTrajs = Samp_parameters[3];
//
//    /* From Dynamics */
//    runDyn = Dyn_parameters[0];
//    total_time = Dyn_parameters[1];
//    dt = Dyn_parameters[2];
//    double junk = Dyn_parameters[3];
//    run_PAC = Dyn_parameters[4];
//    run_energ_conserv = Dyn_parameters[5];
//    tol = Dyn_parameters[6];
//    run_PopAC = Dyn_parameters[7];
//    run_init_PAC = Dyn_parameters[8];
//    interval = Dyn_parameters[9];
//
//    /* Set derived variables */
//    beta = 1.0/temp;
//    beta_nuc = beta/nuc_beads;
//    beta_elec = beta/elec_beads;
//
//    double x_ss = 0.33;
//    double p_ss = 4.0;

    
    
    if(run_samp){
        
        if (my_id == root_process) {
            std::cout << "Begin Sampling" << std::endl;
            std::cout << std::endl;
        }
        
        csrpmd_sampling sampler(my_id,num_procs,root_process);
        sampler.set_sys_vars(nuc_beads,elec_beads,num_states,beta,mass);
        sampler.set_sample_vars(num_trajs,decorr,rate,dt);
        
        sampler.run("test");
        
    }
    
    
    
    
    
//
//
//    /* /////////////////////////////////////////////////////////// */
//                          /* BEGIN PROCESS 3 */
//      /* This process runs Sampling  if requested.*/
//
//    if(runSamp){
//        if (my_id == root_process) {
//            std::cout << "Begin Sampling Simulation" << std::endl;
//            std::cout << std::endl;
//        }
//
//        sampling_mvrpmd sampler(my_id,root_process,num_procs);
//        sampler.initialize_system(nuc_beads,elec_beads,num_states,mass,beta);
//        sampler.initialize_files(readPSV,root);
//
//        clock_t start = clock();
//
//        sampler.run(nuc_ss,x_ss,p_ss,num_trajs,decor_len);
//
//        clock_t end = clock();
//        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//
//        if (my_id == root_process) {
//            std::cout << "\t Sampling simulation time: " << time_taken <<
//            std::endl << std::endl;
//
//            std::cout << "End Sampling Simulation" << std::endl;
//            std::cout << std::endl << std::endl;
//        }
//    }
//
//                            /* END PROCESS 3 */
//  /* /////////////////////////////////////////////////////////// */
//
//
//    /* /////////////////////////////////////////////////////////// */
//                          /* BEGIN PROCESS 4 */
//    /* This process runs the Dynamics simulation if requested.*/
//
//    if(runDyn){
//        if (my_id == root_process) {
//            std::cout << "Begin Dynamics Simulation" << std::endl;
//            std::cout << std::endl << std::endl;
//        }
//
//        /* Setup Dynamics object for simulation. */
//        Dynamics myDyn(num_procs,my_id,root_process,nuc_beads,elec_beads,
//                       num_states,mass,beta_nuc,beta_elec,num_trajs,root);
//
//        myDyn.set_dt(dt);
//        myDyn.set_total_time(total_time);
//
//        clock_t start = clock();
//
//        if(run_energ_conserv){
//            energy_stride = 100;
//            myDyn.energ_conserv(tol,energy_stride);
//        }
//
//        if(run_PAC){
//            myDyn.PAC();
//        }
//
//        if(run_PopAC){
//            bool pac = true;
//            bool bp = false;
//            bool sp = false;
//            bool wp = false;
//
//            int pac_stride = 100;
//            int bp_stride = 100;
//            int sp_stride = 100;
//            int wp_stride = 100;
//
//            myDyn.PopAC(pac,pac_stride,bp,bp_stride,sp,sp_stride,wp,wp_stride);
//        }
//
//        if(run_init_PAC){
//          myDyn.compute_initPAC(interval);
//        }
//
//        clock_t end = clock();
//        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
//
//        if (my_id == root_process) {
//            std::cout << time_taken << std::endl;
//            std::cout << "End Dynamics Simulation" << std::endl;
//        }
//    }
    
                            /* END PROCESS 4 */
    /* /////////////////////////////////////////////////////////// */
    
    MPI_Finalize();
    
    return 0;
}
