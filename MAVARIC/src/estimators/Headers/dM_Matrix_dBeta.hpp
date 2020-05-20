#ifndef dM_Matrix_dBeta_hpp
#define dM_Matrix_dBeta_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>

#include "StateDepPots.h"
#include "M_Matrix.h"

using namespace boost::numeric::ublas;

class dM_Matrix_dBeta{
    
public:
    
    /* Constructor for standard MV-RPMD
     betaN is beta/nuc_beads*/
    dM_Matrix_dBeta(int nuc_beads, int num_states, double betaN, M_Matrix &M_In);
    
    /* Constructor for MTS MV-RPMD
     betaN is beta/elec_beads*/
    dM_Matrix_dBeta(int nuc_beads, int elec_beads, int num_states, double betaN, M_Matrix &M_In);

    
    /* Update dM_dBeta_vec based on the current state of M*/
    void update_dM_dBeta_vec();
    
    /* Return dM_dBeta_vec(alpha)*/
    const matrix<std::complex<double> > &get_dM_dBeta_alpha(int alpha);
    
private:
    
    /* Private Functions. */
    
    /* Update dM_dBeta_alpha to reflect the state of
     Q_alpha. */
    void update_dM_dBeta(int alpha);
    
    /* Private data. */
    int num_states; //number of electronic states
    int nuc_beads; //number of nuclear ring polymer beads
    int elec_beads; //number of electronic ring polymer beads
    
    double ONE_beta; // 1.0/beta
    double betaN; // beta/nuc_beads(elec_beads) depending on constructor

    double ONE_N; // 1.0/nuc_beads(elec_beads) depending on constructor
    
    M_Matrix * M;
    matrix<std::complex<double> > dM_dBeta_alpha;
    vector<matrix<std::complex<double> > > dM_dBeta_vec;
    
};

#endif
