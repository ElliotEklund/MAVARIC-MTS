#ifndef THETA_H
#define THETA_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>

#include "M_Matrix.h"
#include "C_Matrix.h"

using namespace boost::numeric::ublas;

class Theta{
  
public:
    
    Theta(int num_states, int num_beads, double beta_num_beads,
          C_Matrix &C_In, M_Matrix &M_In);

    
    /* Update theta to reflect the state of Q,x, and p. Calls
     update_gamma_mat.
     Q is a vector of bead positions.
    x_mat is a matrix of dim num_beads x num_states; the index
    x_mat[i][j] returns x corresponding to the ith bead and the jth state
    p_mat is a matrix of dim num_beads x num_states; the index
    p_mat[i][j] returns p corresponding to the ith bead and the jth state */
    void update_theta(const vector<double> &Q,const matrix<double> &x,
                      const matrix<double> &p);
    
    /* Return theta; both update_gamma_mat and update_theta are called.
     Q is a vector of bead positions.
    x_mat is a matrix of dim num_beads x num_states; the index
    x_mat[i][j] returns x corresponding to the ith bead and the jth state
    p_mat is a matrix of dim num_beads x num_states; the index
    p_mat[i][j] returns p corresponding to the ith bead and the jth state */
    double& get_theta(const vector<double> &Q,const matrix<double> &x,
                     const matrix<double> &p);
    
    /* Return theta; it is assumed that gamma_mat and theta have been updated
     Q is a vector of bead positions.
    x_mat is a matrix of dim num_beads x num_states; the index
    x_mat[i][j] returns x corresponding to the ith bead and the jth state
    p_mat is a matrix of dim num_beads x num_states; the index
    p_mat[i][j] returns p corresponding to the ith bead and the jth state */
    double & get_theta();
    
    double get_signTheta();
    
private:
    
    /* Private functions. */
    
    /* Update gamma_mat to reflect the state of Q,x, and p.
     Q is a vector of bead positions.
    x_mat is a matrix of dim num_beads x num_states; the index
    x_mat[i][j] returns x corresponding to the ith bead and the jth state
    p_mat is a matrix of dim num_beads x num_states; the index
    p_mat[i][j] returns p corresponding to the ith bead and the jth state */
    void update_gamma_mat(const vector<double> &Q,const matrix<double> &x,
                          const matrix<double> &p);
    
    /* Private data. */
    int num_states; //number of electronic states
    int num_beads; //number of ring polymer beads
    double beta_num_beads; //beta/num_beads
   
    M_Matrix *M;
    C_Matrix *C;

    matrix<std::complex<double> > temp1; //temporary dummy matrix used in update_gamm_mat
    matrix<std::complex<double> > temp2; //temporary dummy matrix used in update_gamm_mat
    matrix<std::complex<double> > myI; //identity matrix
    double theta; // Re(Tr(gamma_mat))
    matrix<std::complex<double> > gamma_mat; // C1 x M1 ... CN x MN
    
};

#endif
