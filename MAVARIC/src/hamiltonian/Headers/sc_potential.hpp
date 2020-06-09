#ifndef sc_potential_hpp
#define sc_potential_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "functions.hpp"

#include <algorithm>

using namespace boost::numeric::ublas;

class sc_potential{
    
public:
    sc_potential(int nuc_beads, int elec_beads, int num_states);
    
    /* Compute non-adiabatic potential contribution to csrpmd hamiltonian*/
    void update_Vsc(const vector<double> &Q,const matrix<double> &x,
                       const matrix<double> &p);
    
    void update_Vsc_dx(const vector<double> &Q,const matrix<double> &p);
    
    void update_Vsc_dp(const vector<double> &Q,const matrix<double> &x);
    
    void update_Vsc_dQ(const vector<double> &Q,const matrix<double> &x,
                       const matrix<double> &p);


/* Mutators */
    double get_Vsc();
    
    /* Update and return pot*/
    double get_Vsc(const vector<double> &Q,const matrix<double> &x,
                   const matrix<double> &p);
    
    const matrix<double> & get_Vsc_dx(const vector<double> &Q,
                                      const matrix<double> &p);
    
    const matrix<double> & get_Vsc_dp(const vector<double> &Q,
                                      const matrix<double> &x);
    
    const vector<double> & get_Vsc_dQ(const vector<double> &Q,
                                      const matrix<double> &x,
                                      const matrix<double> &p);

private:
    
/* Data */
    int nuc_beads,elec_beads,num_states; //system information
    double pot; //final potential
    matrix<double> Vmat; //non-adiabadic potential matrix
    matrix<double> Vmat_dQ; //derivative of non-adiabadic potential matrix wrt Q
    vector<double> Vsc_dQ_temp;
    vector<double> Vsc_dQ;
    matrix<double> Vsc_dx; //derivative of sc potential wrt x
    matrix<double> Vsc_dp; //derivative of sc potential wrt x
    vector<double> Q_trans; //results of W Q
    mapped_matrix<double> W, V; //transformation matrix
    vector<double> x_temp, p_temp; //used int compute_pot

/* Functions */
    
    /* Return Vnn(Q); the nth diagonal element of the non-adiabatic potential
     matrix evaluated at Q. */
    inline double V_d(const double &Q, int n);
    
    inline double V_d_dQ(const double &Q, int n);

    
    /* Return Vnm(Q); the n,m-th non-diagonal element of the non-adiabatic
     potential matrix evaluated at Q. */
    inline double V_od(const double &Q, int n, int m);
    
    inline double V_od_dQ(const double &Q, int n, int m);
    
    /* Update the non-adiabatic coupling matrix, Vmat, given Q*/
    inline void update_Vmat(const double &Q, matrix<double> &Vmat);
    
    inline void update_Vmat_dQ(const double &Q, matrix<double> &Vmat_dQIN);
    
};

#endif
