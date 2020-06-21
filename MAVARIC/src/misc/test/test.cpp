#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "pop_estimators.hpp"
#include "aggregate.hpp"
#include "functions.hpp"

using namespace boost::numeric::ublas;

void write(matrix<double> m, std::string fileName){
    std::ofstream myFile;
    myFile.open(fileName);
    
    int size1 = m.size1();
    int size2 = m .size2();
    
    for (int i=0; i<size1; i++) {
        for (int j=0; j<size2; j++) {
            myFile << m(i,j) << " ";
        }
        myFile << std::endl;
    }
    myFile.close();
}

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
    int num_beads = 6;

    matrix<double> x (num_beads,num_states,0);
    matrix<double> p (num_beads,num_states,0);
    
    double lower_bound = -2;
    double upper_bound = 2;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re(10);

    for(int state=0; state<num_states; state++){
        for(int bead=0; bead<num_beads; bead++){
          x(bead,state) = unif(re);
          p(bead,state) = unif(re);
        }
    }
    
    write(x,"x");
    write(p,"p");
    
    pop_estimators myPops(num_beads,num_states);
    
    vector<double> sc_results (num_states,0);
    vector<double> wig_results (num_states,0);

    sc_results = myPops.sc(x,p);
    wig_results = myPops.wigner(x,p);
    
    aggregate aggregator;

    std::string name1 = "wig";
    std::string name2 = "sc";

    int num_trials = 10;
    
    aggregator.add_calc(name1,num_states,num_trials);
    aggregator.add_calc(name2,num_states,num_trials);

    for (int j=0; j<5; j++) {
        for (int i=0; i<num_trials; i++) {
            wig_results = myPops.wigner(x,p);
            sc_results = myPops.sc(x,p);

            aggregator.collect(name1,i,wig_results,wig_results,1);
            aggregator.collect(name2,i,sc_results,sc_results,1);
        }
    }

    std::string fileName = "ag_results";
    aggregator.merge_collections(0,my_id,fileName);
    
    
    double x = 1.0/0.0;
    double y = 5.0;
    double z = x + y;
    std::cout << z << std::endl;
    
    MPI_Finalize();

    return 0;
}
