#include "dStateDep_dQ.hpp"


#define DE1 0.02
#define DE2 0.02
#define DE3 0.003

#define BETA1 0.4
#define BETA2 0.65
#define BETA3 0.65

#define RE1 4.0
#define RE2 4.5
#define RE3 6.0

#define C1 0.02
#define C2 0.0
#define C3 0.02

#define A12 0.005
#define A13 0.005
#define A23 0.0

//B is used instead of little a
#define B12 32.0
#define B13 32.0
#define B23 0.0

#define R12 3.4
#define R13 4.97
#define R23 0.0

dStateDep_dQ::dStateDep_dQ(int nuc_beads,int num_states)
    :nuc_beads(nuc_beads),num_states(num_states),
     dVdQ(nuc_beads,num_states),
     dVcoup_dQ_vec(nuc_beads)
{
    for (int bead=0; bead<nuc_beads; bead++) {
        dVcoup_dQ_vec(bead) = zero_matrix<double> (num_states,num_states);
    }

    for (int bead=0; bead<nuc_beads; bead++){
        dVdQ(bead,0) = 1.0;
        dVdQ(bead,1) = -1.0;
    }

}

double dStateDep_dQ::dV11_dQ(const double &Q){
    return 1.0;
}

double dStateDep_dQ::dV22_dQ(const double &Q){
    return -1.0;
}

//double dStateDep_dQ::dV33_dQ(const double &Q){
//    return -2.0*DE3*BETA3*exp(BETA3*(RE3-Q))*(exp(BETA3*(RE3-Q))-1.0);
//}

double dStateDep_dQ::dV12_dQ(const double &Q){
    return 0.0;
}

//ouble dStateDep_dQ::dV13_dQ(const double &Q){
//   return -2.0*B13*A13*(Q-R13)*exp(-B13*pow(Q-R13,2));
//
//
//ouble dStateDep_dQ::dV23_dQ(const double &Q){
//   return -2.0*B23*A23*(Q-R23)*exp(-B23*pow(Q-R23,2));
//

void dStateDep_dQ::update_dVdQ(const vector<double> &Q){
    
   /* Fill dV11_dQ*/
   // for (int bead=0; bead<nuc_beads; bead++) {
   //     dVdQ(bead,0) = dV11_dQ(Q(bead));
   // }
   // 
   // /* Fill dV22_dQ*/
   // for (int bead=0; bead<nuc_beads; bead++) {
   //     dVdQ(bead,1) = dV22_dQ(Q(bead));
   // }

   // /* Fill dV22_dQ*/
   // for (int bead=0; bead<nuc_beads; bead++) {
   //     dVdQ(bead,2) = dV33_dQ(Q(bead));
   // }
}

const matrix<double> & dStateDep_dQ::get_dVdQ(){return dVdQ;}

const vector<matrix<double> > & dStateDep_dQ::get_dVcoup_dQ_vec(const vector<double> &Q){
    return dVcoup_dQ_vec;
}
