#ifndef RK4_MVRPMD_hpp
#define RK4_MVRPMD_hpp

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <math.h>

#include "Forces_MTS.hpp"
#include "mvrpmd_mixed_forces.hpp"
#include "csrpmd_forces.hpp"
#include "mv_forces_temp.hpp"

using namespace boost::numeric::ublas;

class RK4_MVRPMD{
    
public:
   RK4_MVRPMD(mv_forces_temp *F_In,int nuc_beads,int elec_beads,int num_states,double dt);
    
    /* Advance the current trajectory one time step forward using RK4-MV-RPMD-SB integrator. */
    void take_step(vector<double> &Q,vector<double> &P,matrix<double> &x,
                   matrix<double> &p);
    
    
private:
    
    /* Variables relevent to the physical system.*/
    int num_states; // number of electronic states
    int nuc_beads; // number of nuclear beads
    int elec_beads;// number of electronib beads
    
    /* Variables relevent to RK4 integrator. */
    double dt; // time step
    double dt_half; // 0.5 * dt
    double coeff1, coeff2, coeff3, coeff4; // coefficients for final step of RK4
    
    /* Four k-parameters for each degree of freedom in the system. */
    vector<double> k1Q, k2Q, k3Q, k4Q;
    vector<double> k1P, k2P, k3P, k4P;
    matrix<double> k1x, k2x, k3x, k4x;
    matrix<double> k1p, k2p, k3p, k4p;
    
    mv_forces_temp *F;
    
    /* Private Functions */
    
    /* All k parameters are filled with the appropriate number of zeros after zero_ks
     is called. */
    void zero_ks();
    
    /* k1 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k1();
    
    /* k2 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k2();
    
    /* k3 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k3();
    
    /* k4 parameters for all degrees of freedom is updated based on the arguments passed. */
    void update_k4();
    
    /* All k parameters must be updated before calling update_final.
     This function ultimately advances all degrees of freedom using the k parameters. */
    void update_final(vector<double> &Q, vector<double> &P,
                      matrix<double> &x, matrix<double> &p);
    
};


#endif
