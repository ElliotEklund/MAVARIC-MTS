#ifndef C_MATRIX_h
#define C_MATRIX_h

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <math.h>

using namespace boost::numeric::ublas;

class C_Matrix{
    
public:
    
    C_Matrix(int num_beads, int num_states, double alpha);
        
    /* C_vec is updated to reflect the state of x_mat and p_mat
     x_mat is a matrix of dim num_beads x num_states; the index
     x_mat[i][j] returns x corresponding to the ith bead and the jth state
     p_mat is a matrix of dim num_beads x num_states; the index
     p_mat[i][j] returns p corresponding to the ith bead and the jth state */
    void update_C_vec(const matrix<double> &x_mat, const matrix<double> &p_mat);
    
    /* Returns C_vec[alpha]. Assumes that C_vec has been updated.*/
    matrix<std::complex<double> >& get_C_alpha(int alpha);

private:
    
    /* Private data. */
    int num_states; //number of states
    int num_beads; //number of ring polymer beads
    double alpha; // mapping variable prefactor
    double p_alpha, x_alpha; //p and x prefactor

    matrix<std::complex<double> > x_plus_ip_mat; // x_mat + ip_mat
    matrix<std::complex<double> > x_min_ip_mat; // x_mat - ip_mat
    matrix<std::complex<double> > x_plus_ip_mat_trans; //transpose of x_min_ip_mat
    
    matrix<std::complex<double> > x_mat_c; // complex vector holding x_mac
    matrix<std::complex<double> > p_mat_c; // complex vector holding x_mac
    
    std::complex<double> unit_complex; //(0.0,1.0)
    matrix<std::complex<double> > half_identity; //0.5 * identity-matrix
    vector<matrix<std::complex<double> > > C_vec; //element alpha corresponds to C_{alpha}
    
};


#endif
