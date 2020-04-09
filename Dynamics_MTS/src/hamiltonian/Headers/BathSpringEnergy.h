#ifndef BATHSPRINGENERGY_H
#define BATHSPRINGENERGY_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include "SpringEnergy.h"

using namespace boost::numeric::ublas;

class BathSpringEnergy{
    
public:
    BathSpringEnergy(int num_beads,int num_modes, double mass,
                     double beta_num_beads);
    
    /* Update enerty to reflect the state of Q_bath.
     Q_bath is a matrix of bath bead positions. The element
     Q_bath[i][j] is the ith bead of the jth bath mode.*/
    void update_bathSpringEnergy(const matrix<double> &Q_bath);

    /* Return energy; this first calls update_bathSpringEnergy.
      Q_bath is a matrix of bath bead positions. The element
      Q_bath[i][j] is the ith bead of the jth bath mode.*/
    double get_bathSpringEnergy(const matrix<double> &Q_bath);

    /* Return energy; this assumes that energy has already been
     updated.
      Q_bath is a matrix of bath bead positions. The element
      Q_bath[i][j] is the ith bead of the jth bath mode.*/
    double get_bathSpringEnergy();

private:
    
    /* Private Data. */
    int num_modes; //number of bath modes
    double energy; //bath spring energy
    
    SpringEnergy V_spring;
};

#endif
