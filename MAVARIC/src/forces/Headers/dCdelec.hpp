#ifndef dCdelec_hpp
#define dCdelec_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <complex>

using namespace boost::numeric::ublas;

class dCdelec{
  
public:
    dCdelec(int elec_beads, int num_states, double alpha);
    
    /* Update dCdx_mat to reflect the state of x and p */
    void update_dCdx_mat(const matrix<std::complex<double> > & x,const matrix<std::complex<double> > & p);
    
    /* Update dCdp_mat to reflect the state of x and p */
    void update_dCdp_mat(const matrix<std::complex<double> > & x,const matrix<std::complex<double> > & p);
    
    /* Return a matrix containing the derivative of C wrt x of state and bead passed as arguments */
    const matrix<std::complex<double> > & get_dC_dx(int bead, int state);

    /* Return a matrix containing the derivative of C wrt p of state and bead passed as arguments */
    const matrix<std::complex<double> > & get_dC_dp(int bead, int state);


private:
    
    /* Private Functions*/
    void update_dCdx(int bead, int state, const matrix<std::complex<double> > & x,
                     const matrix<std::complex<double> > & p,
                     matrix<std::complex<double> > &dC);
    
    void update_dCdp(int bead, int state, const matrix<std::complex<double> > & x,
                     const matrix<std::complex<double> > & p,
                     matrix<std::complex<double> > &dC);

    /* Private Data */
    int num_states;//number of electronic states
    int elec_beads; //number of electronic beads
    double alpha; //mapping variable prefactor
    double x_alpha, p_alpha; //x and p mapping variable prefactor
    std::complex<double> unit_complex;
    
    /* dCdx_mac(i,j) returns a matrix containing the derivative of C wrt x of state j of bead i */
    matrix<matrix<std::complex<double> > > dCdx_mat;
    
    /* dCdx_mac(i,j) returns a matrix containing the derivative of C wrt p of state j of bead i */
    matrix<matrix<std::complex<double> > > dCdp_mat;
    
    matrix<std::complex<double> > x_p_ip; // x plus unit_complex * p
    matrix<std::complex<double> > x_m_ip; // x minus unit_complex * p
    
    matrix<std::complex<double> > p_p_ix; // p plus unit_complex * m
    matrix<std::complex<double> > p_m_ix; // p minus unit_complex * m

    vector<std::complex<double> > dCdx_row;
    vector<std::complex<double> > dCdx_col;

    vector<std::complex<double> > dCdp_row;
    vector<std::complex<double> > dCdp_col;
    
    vector<std::complex<double> > x_state;
    vector<std::complex<double> > p_state;
};

#endif
