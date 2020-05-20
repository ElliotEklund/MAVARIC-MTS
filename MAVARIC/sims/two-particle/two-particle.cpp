#include <iostream>
#include <fstream>
#include <vector>
#include "mpi.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include "two_particle_mc.hpp"
#include "two_particle_sampling.hpp"
#include "two_particle_Dynamics.hpp"
#include "input.hpp"

using namespace boost::numeric::ublas;

int main(int argc, char ** argv) {
    
    int num_procs = 1; //number of processors program is distributed over
    int my_id = 0; //unique id of each processor
    int root_process = 0; //processor 0 is default root process
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    
    MPI_Comm comm = MPI_COMM_WORLD;
    
    input myInput;
    
    std::string root = "/Users/ellioteklund/Desktop/Dynamics_MTS_git/Dynamics_MTS/sims/two-particle/";
    
    /* Vectors used to store parameters from InputFiles */
    std::vector<double> sys1_parameters;
    std::vector<double> sys2_parameters;
    std::vector<double> MC_parameters;
    std::vector<double> Samp_parameters;
    std::vector<double> Dyn_parameters;
    
    int abort = myInput.input_file_handler(root,sys1_parameters,sys2_parameters,MC_parameters,
                                          Samp_parameters,Dyn_parameters);

    if (abort == -1){
        MPI_Finalize();
        return -1;
    }
    
    /* Physical parameters.*/
    double temp, mass;
    int num_beads1, num_beads2;
    double beta;

    /* Monte Carlo parameters. */
    unsigned long long num_steps, esti_rate;
    double ss1, ss2;
    bool writePSV1, readPSV1;
    bool writePSV2, readPSV2;
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
    bool run_init_PAC;
    double interval;
    
    /* From MonteCarlo */
    runMC = MC_parameters[0];
    num_steps = MC_parameters[1];
    esti_rate = MC_parameters[2];
    writePSV1 = MC_parameters[3];
    writeData = MC_parameters[4];
    readPSV1 = MC_parameters[5];
    readData = MC_parameters[6];
    
    writePSV2 = writePSV1;
    readPSV2 = readPSV1;
    
    /* From SystemParameters*/
    mass = sys1_parameters[0];
    num_beads1 = sys1_parameters[1];
    temp = sys1_parameters[2];
    ss1 = sys1_parameters[3];

    /* From SystemParameters*/
    num_beads2 = sys2_parameters[1];
    ss2 = sys2_parameters[3];

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
    
    /* Set derived variables */
    beta = 1.0/temp;

    if (runMC) {
        if (my_id == root_process) {
            std::cout << std::endl << std::endl;
            std::cout << "Begin Monte Carlo Simulation" << std::endl;
            std::cout << std::endl;
        }
        
        std::string output_file = root + "Output/";
        
        two_particle_mc monte_carlo(my_id,root_process,num_procs);
        monte_carlo.initialize_files(readPSV1,readPSV2,readData,
                                     writePSV1,writePSV2,writeData,
                                     output_file);
        monte_carlo.initialize_system(num_beads1,num_beads2,mass,beta);
        
        clock_t start = clock();

        monte_carlo.runSimulation(ss1,ss2,num_steps,esti_rate);
        
        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "\t Monte Carlo simulation time: " << time_taken << std::endl << std::endl;
            std::cout << "End Monte Carlo Simulation" << std::endl;
            std::cout << std::endl;
        }
    }
    
    if (runSamp) {
        if (my_id == root_process) {
            std::cout << std::endl << std::endl;
            std::cout << "Begin Sampling" << std::endl;
            std::cout << std::endl;
        }
        
        std::string output_file = root + "Output/";
        
        two_particle_sampling sampler(my_id,root_process,num_procs);
        sampler.initialize_system(num_beads1,num_beads2,mass,beta);
        sampler.initialize_files(true,true,output_file);
        
        clock_t start = clock();
        
        sampler.runSimulation(ss1,ss2,num_trajs,decor_len);
        
        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "\t Monte Carlo simulation time: " << time_taken << std::endl << std::endl;
            std::cout << "End Monte Carlo Simulation" << std::endl;
            std::cout << std::endl;
        }
    }
    
    if (runDyn) {
            if (my_id == root_process) {
            std::cout << std::endl << std::endl;
            std::cout << "Begin Dynamics" << std::endl;
            std::cout << std::endl;
        }
        
        std::string output_file = root;
        
        two_particle_Dynamics molec_dyn(my_id,root_process,num_procs);
        molec_dyn.initialize_system(num_beads1,num_beads2,mass,beta);
        molec_dyn.initialize_dynamics(dt,total_time,num_trajs,output_file);
        
        clock_t start = clock();

        molec_dyn.run();
        
        clock_t end = clock();
        double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

        if (my_id == root_process) {
            std::cout << "\t Dynamics run time: " << time_taken << std::endl << std::endl;
            std::cout << "End Molecular Dynamics" << std::endl;
            std::cout << std::endl;
        }
    }
    
    MPI_Finalize();
    
    return 0;
}
