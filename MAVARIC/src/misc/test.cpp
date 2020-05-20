#include <iostream>
#include <vector>
#include <string>
#include <map>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "pop_estimators.hpp"
#include "aggregate.hpp"

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

    int num_states = 2;
    int num_beads = 5;

    matrix<double> x (num_beads,num_states,0);
    matrix<double> p (num_beads,num_states,0);

    for(int state=0; state<num_states; state++){
        for(int bead=0; bead<num_beads; bead++){
          x(bead,state) = state + bead;
          p(bead,state) = 5 + state + bead;
        }
    }

    pop_estimators myPops(num_beads,num_states);

    vector<double> sc_results (num_states,0);
    vector<double> wig_results (num_states,0);

    aggregate aggregator;
    
    std::string name1 = "wig";
    std::string name2 = "sc";

    aggregator.add_calc(name1,num_states,10);
    aggregator.add_calc(name2,num_states,10);
    
    for (int j=0; j<5; j++) {
        for (int i=0; i<10; i++) {
            wig_results = myPops.wigner(x,p);
            sc_results = myPops.sc(x,p);
            
            aggregator.collect(name1,i,wig_results,1);
            aggregator.collect(name2,i,sc_results,1);
        }
    }
    
    std::string fileName = "./";
    aggregator.merge_collections(0,my_id,fileName);
    
    MPI_Finalize();

    return 0;
}
