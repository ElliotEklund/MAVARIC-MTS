#include "GTerm.h"

GTerm::GTerm(int num_beads, int num_states)
    :x_squared(num_beads,num_states),
     p_squared(num_beads,num_states)
{}

void GTerm::update_gTerm(const matrix<double> &x,const matrix<double> &p){

    x_squared = element_prod(x,x);
    p_squared = element_prod(p,p);

    energy = 0;
    energy = std::accumulate(x_squared.data().begin(),x_squared.data().end(),energy);
    energy = std::accumulate(p_squared.data().begin(),p_squared.data().end(),energy);
}

double& GTerm::get_gTerm(const matrix<double> &x,const matrix<double> &p){
    update_gTerm(x, p);
    return energy;
}

double& GTerm::get_gTerm(){return energy;}

double GTerm::square::operator() (double Q) const { return Q*Q; }
