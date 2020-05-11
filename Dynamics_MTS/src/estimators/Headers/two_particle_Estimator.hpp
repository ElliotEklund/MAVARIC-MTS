#ifndef two_particle_Estimator_hpp
#define two_particle_Estimator_hpp

#include "SpringEnergy.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

class two_particle_Estimator{
    
public:
    
    two_particle_Estimator(int num_beads1,int num_beads2, double beta,
                           SpringEnergy &V_SpringIn1,SpringEnergy &V_SpringIn2);
    
    /* Return energy estimator. This function assumes that all
     parts of the Hamiltonian have been updated. */
    double get_estimator(const vector<double> &Q1,const vector<double> &Q2);
    
private:
    
    /* Private data.*/
    
    double c1, c2, c12;
    double mom_const; //constant resulting from momentum integration
    int num_beads1, num_beads2;
    double one_num_beads1, one_num_beads2;
    
    mapped_matrix<double> W; //matrix used to calculate coupling energy
    vector<double> Q2_mapped;
    
    double inline V01(const vector<double> &Q1);
    
    double inline V02(const vector<double> &Q2);
    
    double inline Vcouple(const vector<double> &Q1,const vector<double> &Q2);

    
    SpringEnergy * V_spring1;
    SpringEnergy * V_spring2;
    
};
    
#endif
