#include <iostream>
#include <string>
#include "mpi.h"


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "two_particle_Estimator.hpp"

using namespace boost::numeric::ublas;

int main(int argc, char ** argv){
    
    int num_procs = 1; //number of processors program is distributed over
    int my_id = 0; //unique id of each processor
    int root_process = 0; //processor 0 is default root process
    
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    
    MPI_Comm comm = MPI_COMM_WORLD;

    int nuc_beads1 = 8;
    int nuc_beads2 = 4;
    double beta = 1.0;
    double mass = 1.0;
    
    vector<double> Q1(nuc_beads1,0);
    vector<double> Q2(nuc_beads2,0);
    
    for (int i=0; i<nuc_beads1; i++) {
        Q1(i) = i +2*i;
    }
    
    for (int i=0; i<nuc_beads2; i++) {
        Q2(i) = i*i;
    }
    
    SpringEnergy E_spring1(nuc_beads1,mass,beta/nuc_beads1);
    SpringEnergy E_spring2(nuc_beads2,mass,beta/nuc_beads2);
    
    two_particle_Estimator myEstimator(nuc_beads1,nuc_beads2,beta,
                                       E_spring1,E_spring2);
    
    
    E_spring1.get_springEnergy(Q1);
    E_spring2.get_springEnergy(Q2);

    
    double esti = myEstimator.get_estimator(Q1,Q2);
    std::cout << esti << std::endl;
    
    MPI_Finalize();

    return 0;
}
