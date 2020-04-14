#include "StateDepPots.h"

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

StateDepPots::StateDepPots(int num_states, int num_beads, double beta_num_beads)
    :num_states(num_states), num_beads(num_beads), beta_num_beads(beta_num_beads),
     V_mat(num_beads,num_states), V_couple_mat(num_states,num_states)
{
    std::fill(V_mat.data().begin(),V_mat.data().end(),0.0);
    std::fill(V_couple_mat.data().begin(),V_couple_mat.data().end(),1.0);

    /* Constant coupling optimization */
    V_couple_mat(1,0) = 1.0;
    V_couple_mat(0,1) = 1.0;

    V11_vec.resize(num_beads);
    V22_vec.resize(num_beads);
    //V33_vec.resize(num_beads);
}

double StateDepPots::V11::operator() (double x) const { return x;}

double StateDepPots::V22::operator() (double x) const { return -x;}

//double StateDepPots::V33::operator() (double x) const { return DE3*pow((1.0 - exp(-BETA3*(x - RE3))),2) + C3;}


//double StateDepPots::V12(const double &Q){
//    return A12*exp(-B12*(pow(Q - R12,2)));
//}
//
//double StateDepPots::V13(const double &Q){
//    return A13*exp(-B13*(pow(Q - R13,2)));
//}
//
//double StateDepPots::V23(const double &Q){
//    return A23*exp(-B23*(pow(Q - R23,2)));
//}

vector<double> StateDepPots::get_V11_vec(const vector<double> &Q){
    std::transform(Q.begin(),Q.end(),V11_vec.begin(),V11());
    return V11_vec;
}

vector<double> StateDepPots::get_V22_vec(const vector<double> &Q){
    std::transform(Q.begin(),Q.end(),V22_vec.begin(),V22());
    return V22_vec;
}

//vector<double> StateDepPots::get_V33_vec(const vector<double> &Q){
//    std::transform(Q.begin(),Q.end(),V33_vec.begin(),V22());
//    return V33_vec;
//}
matrix<double>& StateDepPots::get_V_couple_mat(const double &Q){
    
   // V_couple_mat(0,1) = V12(Q);
   // V_couple_mat(1,0) = V12(Q);

   // V_couple_mat(0,2) = V13(Q);
   // V_couple_mat(2,0) = V13(Q);

   // V_couple_mat(1,2) = V23(Q);
   // V_couple_mat(2,1) = V23(Q);

    /* Constant coupling optimization */
    return V_couple_mat;
}

matrix<double>& StateDepPots::get_V_mat(const vector<double> &Q){
    
    noalias(column(V_mat, 0)) = Q;//get_V11_vec(Q);
    noalias(column(V_mat, 1)) = -Q;//get_V22_vec(Q);
    //noalias(column(V_mat, 2)) = get_V33_vec(Q);

    return V_mat;
}
