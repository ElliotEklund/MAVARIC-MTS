#include <iostream>
#include <string>
#include "mpi.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "Forces_two_particles.hpp"
#include "dSpring_dQ.hpp"

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
    
    Forces_two_particles F(nuc_beads1,nuc_beads2,mass,beta);
    
    F.update_Forces(Q1,Q2,P1,P2);
    F.print_dHdQ(1);
    F.print_dHdQ(2);
    F.print_dHdP(1);
    F.print_dHdP(2);

    
    MPI_Finalize();

    return 0;
}
