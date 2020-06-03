#ifndef theta_mixed_hpp
#define theta_mixed_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "M_Matrix.h"
#include "C_Matrix.h"
#include "functions.hpp"

using namespace boost::numeric::ublas;

class theta_mixed{
    
public:
    
    theta_mixed(int num_states, int nuc_beads, int elec_beads, C_Matrix &C_In,
              M_Matrix &M_In);
    
    /* Update gamma_mat given Q, x, p.*/
    void update_gamma_mat(const vector<double> &Q,const matrix<double> &x,
                                     const matrix<double> &p);
    
    /* Update theta given Q, x, p. */
    void update_theta(const vector<double> &Q,const matrix<double> &x,
                                 const matrix<double> &p);
    
    /* Update theta given Q, x, p, then return it. */
    double get_theta(const vector<double> &Q,const matrix<double> &x,
                     const matrix<double> &p);
    
    /* Return theta. */
    double get_theta();
    
    /* Return the sign of theta. */
    double get_signTheta();

    /* Update the sign of theta given Q, x, and p, then return the value. */
    double get_signTheta(const vector<double> &Q, const matrix<double> &x,
                         const matrix<double> &p);
    
    /* Return the diagonal components of gamma_mat*/
    matrix<std::complex<double> > get_gamm();

private:
    
    /* Private data. */
    int num_states; //number of electronic states
    int elec_beads; //number of electronic ring polymer beads
    int nuc_beads; //number of nuclear ring polymer beads
    double beta_num_beads; //beta/num_beads
    vector<double> Q_trans; //Q transformed; result of W Q.
    mapped_matrix<double> W; //transformation matrix
    M_Matrix *M;
    C_Matrix *C;
        
    double theta; // Re(Tr(gamma_mat))
    matrix<std::complex<double> > gamma_mat; // C1 x M1 ... CN x MN
    
};
#endif
