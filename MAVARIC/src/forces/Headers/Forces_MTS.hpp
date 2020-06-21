#ifndef Forces_MTS_hpp
#define Forces_MTS_hpp

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "dSpring_dQ.hpp"
#include "dStateIndep_dQ.hpp"
#include "Theta_MTS.hpp"
#include "dTheta_MTS_dQ.hpp"
#include "dTheta_MTS_dElec.hpp"
#include "mv_forces_temp.hpp"

using namespace boost::numeric::ublas;

class Forces_MTS : public mv_forces_temp{
    
public:
    Forces_MTS(int nuc_beads,int elec_beads, int num_states,
               double mass, double beta_nuc_beads,
               Theta_MTS &ThetaMTS_In, dTheta_MTS_dQ &dThetaMTS_dQ_In,
               dTheta_MTS_dElec &dThetaMTS_dElec_In);
    
    
    /* Update dHdP to reflect the state of P*/
    void update_dHdP(const vector<double> &P);
    
    /* Update dHdQ to reflect the state of Q,x,p*/
    void update_dHdQ(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);
    /* Update dHdx to reflect the state of Q,x,p*/
    void update_dHdx(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);
    /* Update dHdp to reflect the state of Q,x,p*/
    void update_dHdp(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);

    void update_Forces(const vector<double> &Q,const vector<double> &P,
                       const matrix<double> &x,const matrix<double> &p);
    
    /* Return dHdQ*/
    const vector<double> & get_dHdQ();
  
    /* Return dHdP*/
    const vector<double> & get_dHdP();
   
    /* Return dHdx*/
    const matrix<double> & get_dHdx();

    /* Return dHdp*/
    const matrix<double> & get_dHdp();
    
    double get_sign(const vector<double> &Q, const matrix<double> &x,
                        const matrix<double> &p);

    /* Functions for debugging*/
    void print_dHdQ(); 
    void print_dHdP(); 
    void print_dHdx(); 
    void print_dHdp(); 
    
private:
    
    /* Private Data. */
    double ONE_mass; //1.0/mass
    double beta_nuc_beads; //beta/nuc_beads
    double ONE_beta_nuc_beads; //1.0/beta_nuc_beads
    double TWO_beta_nuc_beads; //2.0/beta_nuc_beads
    int nuc_beads; //number of nuclear beads
    int elec_beads; //number of electronic beads
    int num_states; //number of electronic states
    double coeff_ONE_theta; //ONE_beta_nuc_beads /theta

    vector<double> dHdQ, dHdP; //derivative of Hamiltonian wrt Q,P
    matrix<double> dHdx, dHdp; //derivative of Hamiltonian wrt x,p

    vector<double> dVspring_dQ_vec;
    vector<double> dV0_dQ_vec;
    vector<double> dThetaMTS_dQ_vec;

    matrix<double> dThetaMTS_dx_vec;
    matrix<double> dThetaMTS_dp_vec;
    
    /* Objects need for force calculation*/
    dSpring_dQ dVspring_dQ;
    dStateIndep_dQ dV0_dQ;
    
    Theta_MTS * ThetaMTS;
    dTheta_MTS_dQ * dThetaMTS_dQ;
    dTheta_MTS_dElec * dThetaMTS_dElec;

};

#endif
