#ifndef two_particle_sampling_hpp
#define two_particle_sampling_hpp

#include <stdio.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <string>

#include "SpringEnergy.h"
#include "two_particle_Hamiltonian.hpp"
#include "two_particle_Estimator.hpp"
#include "SamplingHelper.h"

#include "nr3.h"
#include "ran.h"
#include "gamma.h"
#include "deviates.h"

using namespace boost::numeric::ublas;

class two_particle_sampling {
    
public:
    two_particle_sampling(int my_id, int root_proc, int num_procs);
    
    void initialize_system(int num_beads1IN, int num_beads2IN,
                           double massIN,double betaIN);
    
    void initialize_files(bool readPSV1, bool readPSV2,
                          std::string rootFolder);
    
    void runSimulation(double ss1, double ss2, unsigned long long num_trajs,
                       unsigned long long decorrelation);
    

    /* Generate initial conditions for Q. Q is set after calling*/
    void gen_initQ();

private:
    
    /* Private functions. */
    
    bool check_move(double energy, double energy_prop);

    inline double nuc_dist(const double rn, double step_size);
  
    inline int rand_bead(const Ullong rn, int num_beads);
    
    void gen_initQ(vector<double> &Q, int num_beads, double step_size);
    
    void save_trajs(const vector<double> &v,int size, std::string name);


    /* Private data. */
    int my_id; //unique processor id
    int root_proc; //root processor
    int num_procs; //number of processors
    std::string rootFolder;
    
    int num_beads1, num_beads2; //number of nuclear beads
    double mass, beta;
    bool readPSV1, readPSV2; //read in PSC.txt before simulation if true

    Ran myRand; //NR3 random number generator

    SamplingHelper myHelper;
};

#endif

