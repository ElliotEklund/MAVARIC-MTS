#include "MonteCarlo.hpp"

MonteCarlo::MonteCarlo(int num_beads, double mass, int num_states,
                       double beta_num_beads, double nuc_ss, double elec_ss, std::string root)
    :Q(num_beads), x(num_beads,num_states), p(num_beads,num_states),
     Q_prop(num_beads), x_prop(num_beads,num_states), p_prop(num_beads,num_states),
     mt(0),nuc_dist(-nuc_ss,nuc_ss), elec_dist(-elec_ss,elec_ss),
     unif_1(0.0,1.0), rand_bead(0,num_beads-1),
     num_beads(num_beads), num_states(num_states), beta_num_beads(beta_num_beads),
     V_spring(num_beads, mass, beta_num_beads),
     V0(num_beads, mass),
     G(num_beads, num_states),
     C(num_beads,num_states),M(num_states,num_beads,beta_num_beads),
     theta(num_states, num_beads, beta_num_beads,C,M),
     H(beta_num_beads,V_spring,V0,G,theta),
     dtheta_dBeta(num_beads,num_states,beta_num_beads,C,M),
     Esti(num_beads,beta_num_beads,V_spring,V0,theta,dtheta_dBeta),

    myHelper(root,1,1,1)
{
    gen_initQ();
    gen_initx();
    gen_initp();
}

void MonteCarlo::gen_initQ(){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = nuc_dist(mt);
    }
}

void MonteCarlo::gen_initx(){
    for (int bead=0; bead<num_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            x(bead,state) = elec_dist(mt);
        }
    }
}

void MonteCarlo::gen_initp(){
    for (int bead=0; bead<num_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            p(bead,state) = elec_dist(mt);
        }
    }
}

void MonteCarlo::runSimulation(){
    
    estimator_t.resize(num_steps/esti_rate);
    
    Q_prop = Q;
    x_prop = x;
    p_prop = p;
    
    sys_steps_accpt = 0;
    sys_steps = 0;
    
    elec_steps_accpt = 0;
    elec_steps = 0;
   
    energy = H.get_energy(Q,x,p);
    
    int esti_samples = 0;
    double estimator_total = 0;
    double sgn_total = 0;
   
    estimator = Esti.get_estimator();

    for (int step=0; step<num_steps; step++) {
        if(unif_1(mt) > 0.5){
            //Sample nuclear coordinates
            sample_nuc();
        }
        else{
            //Sample electronic coordinates
            sample_elec();
        }
        estimator_total += estimator;
        sgn_total += sgn_theta;
       
        if(step % esti_rate == 0){
            estimator_t[esti_samples] = estimator_total/(sgn_total + 1);
            esti_samples += 1;
        }
    }
    
    myHelper.print_sys_accpt(sys_steps, sys_steps_accpt);
    myHelper.print_elec_accpt(elec_steps, elec_steps_accpt);
    myHelper.print_avg_energy(estimator_total, sgn_total);
    myHelper.write_estimator(estimator_t,esti_rate);
}

void MonteCarlo::sample_nuc(){
    
    int mcMove = 0;
    
    /* Propose new nuclear moves.*/
    for (int i=0; i<num_beads; i++) {
        mcMove = rand_bead(mt);
        Q_prop(mcMove) = Q(mcMove) + nuc_dist(mt);
    }
    
    double energ_prop = H.get_energy(Q_prop,x,p);
    double esti_prop = Esti.get_estimator();
    int sgn_prop = theta.get_signTheta();
    
    /* Accept new system moves if energ_prop < energy*/
    if(energ_prop < energy){
        Q = Q_prop;
        energy = energ_prop;
        estimator = esti_prop;
        sgn_theta = sgn_prop;
        sys_steps_accpt += 1;
    }
    
    /* Accept new system moves if inequality is met*/
    else if (unif_1(mt) <= exp(-beta_num_beads * (energ_prop - energy))){
        Q = Q_prop;
        energy = energ_prop;
        estimator = esti_prop;
        sgn_theta = sgn_prop;
        sys_steps_accpt += 1;
    }
    else{Q_prop = Q;}
    
    sys_steps += 1;
}

void MonteCarlo::sample_elec(){
    
    int mcMove = 0;
    
    //Propose new electronic moves.
    for (int bead=0; bead<num_beads; bead++) {
        mcMove = rand_bead(mt);
        
        for(int state=0; state<num_states; state++){
            x_prop(mcMove,state) = x(mcMove,state) + elec_dist(mt);
            p_prop(mcMove,state) = p(mcMove,state) + elec_dist(mt);
        }
    }
    
    double energ_prop = H.get_energy(Q,x_prop,p_prop);
    double esti_prop = Esti.get_estimator();
    int sgn_prop = theta.get_signTheta();
    
    /* Accept new electronic move if energ_prop < energy */
    if(energ_prop < energy){
        x = x_prop;
        p = p_prop;
        energy = energ_prop;
        estimator = esti_prop;
        sgn_theta = sgn_prop;
        elec_steps_accpt += 1;
    }
    
    /* Accept new system moves if inequality is met*/
    else if (unif_1(mt) <= exp(-beta_num_beads * (energ_prop - energy))){
        x = x_prop;
        p = p_prop;
        energy = energ_prop;
        estimator = esti_prop;
        sgn_theta = sgn_prop;
        elec_steps_accpt += 1;
    }
    
    else{
        x_prop = x;
        p_prop = p;
    }
    
    elec_steps += 1;
}

void MonteCarlo::set_num_steps(unsigned long long num_steps_In){
    num_steps = num_steps_In;
}

void MonteCarlo::set_esti_rate(unsigned long long esti_rate_In){
    esti_rate = esti_rate_In;
}

