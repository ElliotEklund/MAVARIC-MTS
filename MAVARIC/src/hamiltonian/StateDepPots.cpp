#include "StateDepPots.h"

StateDepPots::StateDepPots(int num_states, int num_beads, double beta_num_beads)
    :num_states(num_states), num_beads(num_beads), beta_num_beads(beta_num_beads),
     V_mat(num_beads,num_states,0), V_couple_mat(num_states,num_states,0)
{
    /* Constant coupling optimization */
    V_couple_mat(1,0) = 0.1;
    V_couple_mat(0,1) = 0.1;

    V11_vec.resize(num_beads);
    V22_vec.resize(num_beads);
}
double StateDepPots::V11::operator() (double x) const { return x;}

double StateDepPots::V22::operator() (double x) const { return -x;}

vector<double> StateDepPots::get_V11_vec(const vector<double> &Q){
    std::transform(Q.begin(),Q.end(),V11_vec.begin(),V11());
    return V11_vec;
}
vector<double> StateDepPots::get_V22_vec(const vector<double> &Q){
    std::transform(Q.begin(),Q.end(),V22_vec.begin(),V22());
    return V22_vec;
}

matrix<double>& StateDepPots::get_V_couple_mat(const double &Q){
    return V_couple_mat;
}
matrix<double>& StateDepPots::get_V_mat(const vector<double> &Q){
    noalias(column(V_mat, 0)) = Q;//get_V11_vec(Q);
    noalias(column(V_mat, 1)) = -Q;//get_V22_vec(Q);
    return V_mat;
}
