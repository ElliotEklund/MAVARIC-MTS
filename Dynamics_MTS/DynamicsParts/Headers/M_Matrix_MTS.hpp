#ifndef M_Matrix_MTS_hpp
#define M_Matrix_MTS_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include <algorithm>

#include "M_Matrix.h"

using namespace boost::numeric::ublas;

class M_Matrix_MTS{
    
public:
    
    /* beta_num_bead is beta/elec_beads for MTS*/
    M_Matrix_MTS(int nuc_beads,int elec_beads,int num_states,
                 M_Matrix & M_In);
    
    /* Update M_MTS_vec; this implicitly calls M->update_M_vec*/
    void update_M_MTS_vec(const vector<double> &Q);
    
    /* Update M_MTS_vec; this assumes M->update_M_vec* has already been called */
    void update_M_MTS_vec();
    
    /* Return M_MTS_vec(alpha) */
    matrix<std::complex<double> >& get_M_MTS_alpha(int alpha);
    
private:
    
    /* Private data*/
    int nuc_beads; //number of nuclear beads
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    
    int gamma; //nuc_beads/elec_beads;
    
    /* M_MTS_vec(i)(j,k) = element (j,k) of bead i */
    vector<matrix<std::complex<double> > > M_MTS_vec;
    
    M_Matrix *M;
    
};

#endif
