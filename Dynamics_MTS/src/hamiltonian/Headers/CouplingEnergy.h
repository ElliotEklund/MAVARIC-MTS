#ifndef COUPLINGENERGY_H
#define COUPLINGENERGY_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <iostream>

using namespace boost::numeric::ublas;

class CouplingEnergy{
    
public:
    
    CouplingEnergy(int bath_modes, int num_beads, double mass,
                   vector<double> cs, vector<double> ws);
    
    /* Update energy to reflect the state of Q and Q_bath.
     Q is a vector of bead positions corresponding the system.
     Q_bath is a matrix of bead positions corresponding to the bath, it has
     dimensions (num_beads x num_states). Q_bath[i][j] returns the position
     of the ith bead of the jth bath mode.*/
    void update_couplingEnergy(const vector<double> &Q, const matrix<double> &Q_bath);
    
    /* Return the energy associated with the system-bath coupling.
     update_couplingEnergy is called first
     Q is a vector of bead positions corresponding the system.
     Q_bath is a matrix of bead positions corresponding to the bath, it has
     dimensions (num_beads x num_states). Q_bath[i][j] returns the position
     of the ith bead of the jth bath mode.*/
    double get_couplingEnergy(const vector<double> &Q, const matrix<double> &Q_bath);
    
    /* Return the energy associated with the system-bath coupling.
     It is assummed that energy has already been updated.
     Q is a vector of bead positions corresponding the system.
     Q_bath is a matrix of bead positions corresponding to the bath, it has
     dimensions (num_beads x num_states). Q_bath[i][j] returns the position
     of the ith bead of the jth bath mode.*/
    double get_couplingEnergy();
    
private:
    
    /* Private Structs. */
    struct square {
        double operator() (double x) const;
    };
    
    /* Private data. */
    int bath_modes; //number of bath modes
    int num_beads; //number of ring polymer beads
    double mass; //bath mass
    double energy; //coupling energy
    vector<double> cs; //bath coupling strengths
    vector<double> ws; //bath mode frequencies
    
    vector<double> wwm; //ws*ws*mass; used for optimization
    vector<double> c_wwm; // cs/(mass*ws*ws); used for optimization
    
    vector<double> diff_sq; // (Q_bath[i] - Q)^2; intermediate value
};

#endif
