#ifndef COUPLINGENERGY_H
#define COUPLINGENERGY_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <iostream>

using namespace boost::numeric::ublas;

class CouplingEnergy{
    
public:
    
    CouplingEnergy(int bath_modes, int sys_beads, int bath_beads,
                   double mass, vector<double> cs, vector<double> ws);
    
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
    int sys_beads; //number of syste ring polymer beads
    int bath_beads; //number of bath ring polymer beads
    int r; //sys_beads/bath_beads
    double mass; //bath mass
    double energy; //coupling energy
    vector<double> cs; //bath coupling strengths
    vector<double> ws; //bath mode frequencies
    vector<double> wwm; //ws*ws*mass; used for optimization
    vector<double> c_wwm; // cs/(mass*ws*ws); used for optimization
    mapped_matrix<double> W; //matrix used to calculate coupling energy
    vector<double> ones_bath_beads; //unit vector of length num_beads
    vector<double> Q_sub; //result of matrix product of W and Q
    matrix<double> cQ_sub; //result of outter product of Q_sub x c_wwm
    matrix<double> dif;
    matrix<double> dif_sq; //result of squaring each element in dif
    vector<double> sum_dif_sq; //result of summing across columns of dif_sq

};

#endif
