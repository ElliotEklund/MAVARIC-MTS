#ifndef SPRINGENERGY_H
#define SPRINGENERGY_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

class SpringEnergy{
    
public:
    SpringEnergy(int num_beads, double mass, double beta_num_beads);
    
    /* Update energy to reflect the state of Q.
     Q is a vector of bead positions.*/
    void update_springEnergy(const vector<double> &Q);
    
    /* Return the spring energy; this first calls update_springEnergy
     given Q.
     Q is a vector of bead positions.*/
    double& get_springEnergy(const vector<double> &Q);
    
    /* Return energy; this assumes update_springEnergy has already
     been called. */
    double& get_springEnergy();
    
private:
        
    /* Private data. */
    const int num_beads; //number of beads
    const double mass; //mass
    const double beta_num_beads; //beta/num_beads
    const double preFactor; // mass/(2.0 * beta_num_beads^(2))
    double energy; //spring energy
    vector<double> Q_diff; //Qi - Qi+1
    mapped_matrix<double> W; //matrix used to calculate spring energy
};

#endif
