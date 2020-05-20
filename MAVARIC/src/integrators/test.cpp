#include <iostream>
#include <string>
#include "mpi.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "VV_two_particle.hpp"

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
    double dt = 0.1;
    
    vector<double> Q1(nuc_beads1,0);
    vector<double> Q2(nuc_beads2,0);
    vector<double> P1(nuc_beads1,0);
    vector<double> P2(nuc_beads2,0);
    
    for (int i=0; i<nuc_beads1; i++) {
        Q1(i) = i +2*i;
        P1(i) = i*i;
    }
    
    for (int i=0; i<nuc_beads2; i++) {
        Q2(i) = i*i;
        P2(i) = i +2*i;
    }
    
    std::cout << Q1 << std::endl;
    std::cout << Q2 << std::endl;
    std::cout << P1 << std::endl;
    std::cout << P2 << std::endl;
    
    VV_two_particle VV(nuc_beads1,nuc_beads2,mass,beta,dt);
    VV.step(Q1,Q2,P1,P2);
    
    std::cout << Q1 << std::endl;
    std::cout << Q2 << std::endl;
    std::cout << P1 << std::endl;
    std::cout << P2 << std::endl;
    
    MPI_Finalize();

    return 0;
}
