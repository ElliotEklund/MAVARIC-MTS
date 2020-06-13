#include "GTerm.h"

GTerm::GTerm(int num_beads, int num_states, double alpha)
    :x_squared(num_beads,num_states),
     p_squared(num_beads,num_states),
     alpha(alpha)
{
    energy = 0;
    x_alpha = alpha;
    p_alpha = 1.0/alpha;
}
void GTerm::update_gTerm(const matrix<double> &x,const matrix<double> &p){
    
    x_squared = element_prod(x,x);
    p_squared = element_prod(p,p);

    double x_sum = 0;
    double p_sum = 0;
    
    x_sum = std::accumulate(x_squared.data().begin(),x_squared.data().end(),x_sum);
    p_sum = std::accumulate(p_squared.data().begin(),p_squared.data().end(),p_sum);
    
    energy = x_alpha*x_sum + p_alpha*p_sum;
}
double& GTerm::get_gTerm(const matrix<double> &x,const matrix<double> &p){
    update_gTerm(x, p);
    return energy;
}

double& GTerm::get_gTerm(){return energy;}

double GTerm::square::operator() (double Q) const { return Q*Q; }
