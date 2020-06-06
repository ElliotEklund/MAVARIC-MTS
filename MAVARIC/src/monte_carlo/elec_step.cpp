#include "elec_step.hpp"

elec_step::elec_step(int my_id,int num_procs,int root_proc,int num_beads,
                     int num_states,double betaIN)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc),
     myRand(time(NULL) + my_id),
     num_beads(num_beads),
     num_states(num_states),
     x_prop(num_beads,num_states,0),
     p_prop(num_beads,num_states,0)
{
    set_beta(betaIN);
    x_steps_accepted = 0;
    x_steps_total = 0;
    p_steps_accepted = 0;
    p_steps_total = 0;
}
void elec_step::step_x(double energyIN, const vector<double> &Q, matrix<double> &x,
                       const matrix<double> &p){
    
    energy = energyIN;
    double mcMove=0;
    double energy_prop = 0;
    bool accept_move = false;
    
    for (int bead=0; bead<num_beads; bead++) {
        mcMove = rand_bead(myRand.int64(),num_beads);
        
        for (int state=0; state<num_states; state++) {
            x_prop(mcMove,state) = x(mcMove,state) +
                                    step_dist(myRand.doub(),x_ss);
        }
        
        energy_prop = H->get_energy(Q,x_prop,p);
        accept_move = check_move(energy,energy_prop);
        
        if (accept_move) {
            x = x_prop;
            energy = energy_prop;
            x_steps_accepted += 1;
        }
        else{
            x_prop = x;
        }
        x_steps_total += 1;
    }
}
void elec_step::step_p(double energyIN, const vector<double> &Q,
                       const matrix<double> &x, matrix<double> &p){
    
    energy = energyIN;
    double mcMove=0;
    double energy_prop = 0;
    bool accept_move = false;
    
    for (int bead=0; bead<num_beads; bead++) {
        mcMove = rand_bead(myRand.int64(),num_beads);
        
        for (int state=0; state<num_states; state++) {
            p_prop(mcMove,state) = p(mcMove,state) +
                                    step_dist(myRand.doub(),p_ss);
        }
        
        energy_prop = H->get_energy(Q,x,p_prop);
        accept_move = check_move(energy,energy_prop);
        
        if (accept_move) {
            p = p_prop;
            energy = energy_prop;
            p_steps_accepted += 1;
        }
        else{
            p_prop = p;
        }
        p_steps_total += 1;
    }
}
bool elec_step::check_move(double energy, double energy_prop){

    double delta_energy = energy_prop - energy;
    double accept_move = false;
    /* Accept new system moves if energ_prop < energy*/
    if(delta_energy < 0){
        accept_move = true;
    }
    /* Accept new system moves if inequality is met*/
    else if (myRand.doub() <= exp(-beta * delta_energy/num_beads)){
        accept_move = true;
    }
    return accept_move;
}
void elec_step::set_ss(double x_ssIN, double p_ssIN){
    x_ss = x_ssIN;
    p_ss = p_ssIN;
}
inline int elec_step::rand_bead(const Ullong rn, int num_beads){
    return rn % num_beads;
}
inline double elec_step::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
void elec_step::set_hamiltonian(mvrpmd_mixed &H_IN){H = &H_IN;}

void elec_step::set_beta(double betaIN){beta = betaIN;}

void elec_step::set_energy(double energyIN){energy = energyIN;}

double elec_step::get_energy(){return energy;}

double elec_step::get_x_steps_total(){return x_steps_total;}

double elec_step::get_x_steps_accpt(){return x_steps_accepted;}

double elec_step::get_p_steps_total(){return p_steps_total;}

double elec_step::get_p_steps_accpt(){return p_steps_accepted;}
