#ifndef ABM_MVRPMD_hpp
#define ABM_MVRPMD_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "Forces_MTS.hpp"
#include "RK4_MVRPMD.hpp"
#include "mvrpmd_mixed_forces.hpp"
#include "mv_forces_temp.hpp"
#include <iomanip>
#include <math.h>

using namespace boost::numeric::ublas;

class ABM_MVRPMD{
  
public:
    
    ABM_MVRPMD(mv_forces_temp &F_In,double dt, int num_states,
               int nuc_beads, int elec_beads);
    
    /* Initialize ABM by taking backward steps with RK4*/
    void initialize_rk4(vector<double> &Q,vector<double> &P,
                          matrix<double> &x,matrix<double> &p);
    
    /* Move phase space variables forward by dt. */
    void take_step(vector<double> &Q, vector<double> &P,
                               matrix<double> &x, matrix<double> &p);
    
private:
    
    /* Variables relevent to the physical system.*/
    int num_states; // number of electronic states
    int nuc_beads; // number of nuclear beads
    int elec_beads; //number of electronic beads
    
    /* Variables relevent to RK4 integrator. */
    double dt; // time step
    
    double h1_p; //predictor coefficient 1
    double h2_p; //predictor coefficient 2
    double h3_p; //predictor coefficient 3
    double h4_p; //predictor coefficient 4
    
    double h1_c; //corrector coefficient 1
    double h2_c; //corrector coefficient 2
    double h3_c; //corrector coefficient 3
    double h4_c; //corrector coefficient 4
    
    vector<double> f_Q_3; //derivative of Hamiltonian wrt Q, evaluated at time t-3
    vector<double> f_Q_2; //derivative of Hamiltonian wrt Q, evaluated at time t-2
    vector<double> f_Q_1; //derivative of Hamiltonian wrt Q, evaluated at time t-1
    vector<double> f_Q_0; //derivative of Hamiltonian wrt Q, evaluated at time t
    vector<double> f_Q_p1; //derivative of Hamiltonian wrt Q, evaluated at time t+1
    
    vector<double> f_P_3; //derivative of Hamiltonian wrt P, evaluated at time t-3
    vector<double> f_P_2; //derivative of Hamiltonian wrt P, evaluated at time t-2
    vector<double> f_P_1; //derivative of Hamiltonian wrt P, evaluated at time t-1
    vector<double> f_P_0; //derivative of Hamiltonian wrt P, evaluated at time t
    vector<double> f_P_p1; //derivative of Hamiltonian wrt P, evaluated at time t+1
    
    matrix<double> f_x_3; //derivative of Hamiltonian wrt x, evaluated at time t-3
    matrix<double> f_x_2; //derivative of Hamiltonian wrt x, evaluated at time t-2
    matrix<double> f_x_1; //derivative of Hamiltonian wrt x, evaluated at time t-1
    matrix<double> f_x_0; //derivative of Hamiltonian wrt x, evaluated at time t
    matrix<double> f_x_p1; //derivative of Hamiltonian wrt x, evaluated at time t+1
    
    matrix<double> f_p_3; //derivative of Hamiltonian wrt p, evaluated at time t-3
    matrix<double> f_p_2; //derivative of Hamiltonian wrt p, evaluated at time t-2
    matrix<double> f_p_1; //derivative of Hamiltonian wrt p, evaluated at time t-1
    matrix<double> f_p_0; //derivative of Hamiltonian wrt p, evaluated at time t
    matrix<double> f_p_p1; //derivative of Hamiltonian wrt p, evaluated at time t+1
    
    vector<double> pred_Q; //Q given by predictor phase
    vector<double> pred_P; //P given by predictor phase
    matrix<double> pred_x; //x given by predictor phase
    matrix<double> pred_p; //p given by predictor phase
    
    
    /* Private Functions */
    
    /* Prediction step; pred_var is updated after this call*/
    void predict(const vector<double> &Q,const vector<double> &P,
                 const matrix<double> &x,const matrix<double> &p);
    
    /* Correction step; phase space variables are updated after this call*/
    void correct(vector<double> &Q, vector<double> &P,
                             matrix<double> &x, matrix<double> &p);
    
    /* f_var_3 updated after call */
    void update_f_3(const vector<double> &Q,const vector<double> &P,
                    const matrix<double> &x,const matrix<double> &p);
    
    /* f_var_2 updated after call */
    void update_f_2(const vector<double> &Q,const vector<double> &P,
                    const matrix<double> &x,const matrix<double> &p);
    
    /* f_var_1 updated after call */
    void update_f_1(const vector<double> &Q,const vector<double> &P,
                    const matrix<double> &x,const matrix<double> &p);
    
    /* f_var_0 updated after call */
    void update_f_0(const vector<double> &Q,const vector<double> &P,
                    const matrix<double> &x,const matrix<double> &p);
    
    /* f_var_p1 updated after call */
    void update_f_p1(const vector<double> &Q,const vector<double> &P,
                    const matrix<double> &x,const matrix<double> &p);
    
    mv_forces_temp *F;
};


#endif
