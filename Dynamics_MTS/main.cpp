#include <iostream>
#include <fstream>
#include "mpi.h"

#include "Dynamics.hpp"
#include "MainHlpr.hpp"
#include "MonteCarlo_MTS.hpp"
#include "Sampling_MTS.hpp"

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
    
    MainHlpr myHlpr;
    
    std::string root = "/Users/ellioteklund/Desktop/Dynamics_MTS_git/Dynamics_MTS/";
    int abort = myHlpr.input_file_handler(root,sys_parameters,elec_parameters,MC_parameters,
                                          Samp_parameters,Dyn_parameters);
    
    if (abort == -1){
        return -1;
    }
    
    /* Physical parameters.*/
    double temp, mass;
    int nuc_beads,num_states, elec_beads;
    double beta, beta_nuc, beta_elec;
    
    /* Monte Carlo parameters. */
    unsigned long long num_steps, esti_rate;
    double nuc_ss, elec_ss;
    bool writePSV, readPSV;
    bool readData, writeData;
    bool runMC;
    
    /* Sampling parameters.*/
    bool runSamp, saveTrajs;
    int num_trajs;
    unsigned long long decor_len;
    
    /* Dynamics Variables */
    bool runDyn;
    double dt;
    double total_time;
    double tol;
    int energy_stride;
    
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
    elec_ss = elec_parameters[2];
    
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
            std::cout << "Begin Monte Carlo Simulation" << std::endl;
            std::cout << std::endl << std::endl;
        }

        MonteCarlo_MTS myMC_MTS(nuc_beads,elec_beads,mass,num_states,
                                beta_nuc,beta_elec,nuc_ss,elec_ss,root);

        myMC_MTS.set_num_steps(num_steps);
        myMC_MTS.set_esti_rate(esti_rate);
        myMC_MTS.set_write_PSV(writePSV);
        myMC_MTS.set_read_PSV(readPSV);
        myMC_MTS.set_read_Data(readData);
        myMC_MTS.set_write_Data(writeData);

        clock_t start = clock();

        myMC_MTS.runSimulation();

        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << time_taken << std::endl;
            std::cout << "End Monte Carlo Simulation" << std::endl;
            std::cout << std::endl << std::endl;
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
            std::cout << std::endl << std::endl;
        }

        Sampling_MTS mySamp(nuc_beads,elec_beads,mass,num_states,
                            beta_nuc,beta_elec,nuc_ss,elec_ss,root);

        mySamp.set_decor_len(decor_len);
        mySamp.set_num_samples(num_trajs);

        clock_t start = clock();

        mySamp.runSimulation();

        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << time_taken << std::endl;
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
        Dynamics myDyn(num_procs,my_id,root_process,nuc_beads,elec_beads,
                       num_states,mass,beta_nuc,beta_elec,num_trajs,root);

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
            std::cout << "End Dynamics Simulation" << std::endl;
        }
    }
    
                            /* END PROCESS 4 */
    /* /////////////////////////////////////////////////////////// */
    
    MPI_Finalize();
    
    return 0;
}
