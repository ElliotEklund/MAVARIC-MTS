#ifndef two_particle_mc_hpp
#define two_particle_mc_hpp

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

//#include "MonteCarloHelper.h"
#include "SpringEnergy.h"
#include "two_particle_Hamiltonian.hpp"
#include "two_particle_Estimator.hpp"
#include "MonteCarloHelper.h"

#include "nr3.h"
#include "ran.h"

using namespace boost::numeric::ublas;

class two_particle_mc {
    
public:
    two_particle_mc(int my_id, int root_proc, int num_procs);
    
    void initialize_system(int num_beads1IN, int num_beads2IN,
                           double massIN,double betaIN);
    
    void initialize_files(bool readPSV1, bool readPSV2, bool readData,
                          bool writePSV1, bool writePSV2, bool writeData,
                          std::string rootFolder);
    
    void runSimulation(double ss1, double ss2, unsigned long long num_steps,
                       unsigned long long stride);
    

    /* Generate initial conditions for Q. Q is set after calling*/
    void gen_initQ();

    /* Mutators: Function is self-evident */
    void set_num_steps(unsigned long long num_steps_In);
    
    void set_esti_rate(unsigned long long esti_rate_In);
    
    void set_write_PSV(bool set_In);
    
    void set_read_PSV(bool set_In);
    
    void set_read_Data(bool set_In);

    void set_write_Data(bool set_In);
    
    void print_avg_energy(double estimator_total, double sgn_total);

private:
    
    /* Private functions. */
    
    bool check_move(double energy, double energy_prop);

    inline double nuc_dist(const double rn, double step_size);
  
    inline int rand_bead(const Ullong rn, int num_beads);
    
    void gen_initQ(vector<double> &Q, int num_beads, double step_size);



    /* Private data. */
    int my_id; //unique processor id
    int root_proc; //root processor
    int num_procs; //number of processors
    std::string rootFolder;
    
    int num_beads1, num_beads2; //number of nuclear beads
    double mass, beta;
    
    bool writePSV1, writePSV2; //write PSV to file after simulation if true
    bool readPSV1, readPSV2; //read in PSC.txt before simulation if true
    bool readData; //read in mcData.txt before simulation if true
    bool writeData; //write to mcData.txt after simulation if true

    Ran myRand; //NR3 random number generator

    MonteCarloHelper myHelper;
};

#endif

