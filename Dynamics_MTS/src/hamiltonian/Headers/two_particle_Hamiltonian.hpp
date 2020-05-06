/*
 Two particle hamiltonian defines the hamiltonian for two particles described by
 ring polymer, coupled together with a simple quadratic coupling. It is assumed that
 N1 >= N2. This is intended as a test class to sort out the 1/N problem.
 */

#ifndef two_particle_Hamiltonian_hpp
#define two_particle_Hamiltonian_hpp

#include "SpringEnergy.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

class two_particle_Hamiltonian{
    
public:
    
    two_particle_Hamiltonian(int num_beads1, int num_beads2, SpringEnergy &V_springIn1,
                             SpringEnergy &V_springIn2);
    
    /*
     Calculate the energy corresponding to the two particle hamiltonian
     */
    double get_energy(const vector<double> &Q1,const vector<double> &Q2);

//    double get_energy();
//
//    double get_energy_dyn(double mass,const vector<double> &Q, const vector<double> &P);
    
private:
    
    /* Private data*/
    int num_beads1, num_beads2;
    double one_num_beads1, one_num_beads2;
    
    SpringEnergy * V_spring1;
    SpringEnergy * V_spring2;

    double c1, c2, c12;
    
    mapped_matrix<double> W; //matrix used to calculate coupling energy
    vector<double> Q2_mapped;

    
    double inline V01(const vector<double> &Q1);
    
    double inline V02(const vector<double> &Q2);

    double inline Vcouple(const vector<double> &Q1,const vector<double> &Q2);
    
};

#endif
