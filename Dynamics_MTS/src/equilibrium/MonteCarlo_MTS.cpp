#include "MonteCarlo_MTS.hpp"

MonteCarlo_MTS::MonteCarlo_MTS(int my_id, int root_proc, int num_procs, int num_beads, int elec_beads,double mass, int num_states,
                               double beta_num_beads, double beta_elec_beads, double nuc_ss, double elec_ss, std::string root)
    :/* Initialize parameters*/
     my_id(my_id), root_proc(root_proc), num_procs(num_procs),
     num_beads(num_beads), num_states(num_states),
     beta_num_beads(beta_num_beads), elec_beads(elec_beads),
    
     /* Initialize phase space variable vectors*/
     Q(num_beads), x(elec_beads,num_states), p(elec_beads,num_states),
     Q_prop(num_beads), x_prop(elec_beads,num_states), p_prop(elec_beads,num_states),

     /* Initialize random number generators */
     mt(myRand()), nuc_dist(-nuc_ss,nuc_ss), elec_dist(-elec_ss,elec_ss),
     unif_1(0.0,1.0), rand_bead(0,num_beads-1), rand_elec_bead(0,elec_beads-1),

     /* Initialize Hamiltonian parts and Estimator*/
     V_spring(num_beads, mass, beta_num_beads),
     V0(num_beads, mass),
     G(elec_beads, num_states),
     C(elec_beads,num_states),
     M(num_states,num_beads,beta_elec_beads),
     M_MTS(num_beads,elec_beads,num_states,M),
     thetaMTS(num_states, elec_beads,C,M_MTS),
     dthetaMTS_dBeta(num_beads,elec_beads,num_states,beta_elec_beads,C,M,M_MTS),
     H_MTS(beta_num_beads,V_spring,V0,G,thetaMTS),
     Esti_MTS(num_beads,beta_num_beads,V_spring,V0,thetaMTS,dthetaMTS_dBeta),

     myHelper(root,my_id,num_procs,root_proc)

{
    /* Generate random initial PSV. If readPSV is true, these will be overwritten
     with saved values during the call to runSimulation.*/
    gen_initQ();
    gen_initx();
    gen_initp();
}

void MonteCarlo_MTS::gen_initQ(){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = nuc_dist(mt);
    }
}

void MonteCarlo_MTS::gen_initx(){
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            x(bead,state) = elec_dist(mt);
        }
    }
}

void MonteCarlo_MTS::gen_initp(){
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            p(bead,state) = elec_dist(mt);
        }
    }
}

void MonteCarlo_MTS::runSimulation(){
   
    int esti_samples = 0;
    double estimator_total = 0;
    double sgn_total = 0;
    
    if(readPSV){
        myHelper.read_PSV(num_beads, elec_beads, num_states, Q, x, p);
    }
    
    if (readData) {
        myHelper.read_MC_data(sgn_total, estimator_total);
    }

    estimator_t.resize(num_steps/esti_rate);
    
    Q_prop = Q;
    x_prop = x;
    p_prop = p;

    sys_steps_accpt = 0;
    sys_steps = 0;

    elec_steps_accpt = 0;
    elec_steps = 0;

    energy = H_MTS.get_energy(Q,x,p);
    estimator = Esti_MTS.get_estimator();

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
    
    if(writePSV){
        myHelper.write_PSV(num_beads, elec_beads, num_states, Q, x, p);
    }

    if (writeData) {
        myHelper.write_MC_data(sgn_total, estimator_total);
    }
    
}

void MonteCarlo_MTS::sample_nuc(){

    int mcMove = 0;

    /* Propose new nuclear moves.*/
    for (int i=0; i<num_beads; i++) {
        mcMove = rand_bead(mt);
        Q_prop(mcMove) = Q(mcMove) + nuc_dist(mt);
    }

    double energ_prop = H_MTS.get_energy(Q_prop,x,p);
    double esti_prop = Esti_MTS.get_estimator();
    int sgn_prop = thetaMTS.get_signTheta();

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

void MonteCarlo_MTS::sample_elec(){

    int mcMove = 0;

    //Propose new electronic moves.
    for (int bead=0; bead<elec_beads; bead++) {
        mcMove = rand_elec_bead(mt);

        for(int state=0; state<num_states; state++){
            x_prop(mcMove,state) = x(mcMove,state) + elec_dist(mt);
            p_prop(mcMove,state) = p(mcMove,state) + elec_dist(mt);
        }
    }

    double energ_prop = H_MTS.get_energy(Q,x_prop,p_prop);
    double esti_prop = Esti_MTS.get_estimator();
    int sgn_prop = thetaMTS.get_signTheta();

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

void MonteCarlo_MTS::set_num_steps(unsigned long long num_steps_In){
    num_steps = num_steps_In;
}

void MonteCarlo_MTS::set_esti_rate(unsigned long long esti_rate_In){
    esti_rate = esti_rate_In;
}

void MonteCarlo_MTS::set_write_PSV(bool set_In){
    writePSV = set_In;
}

void MonteCarlo_MTS::set_read_PSV(bool set_In){
    readPSV = set_In;
}

void MonteCarlo_MTS::set_read_Data(bool set_In){
    readData = set_In;
}

void MonteCarlo_MTS::set_write_Data(bool set_In){
    writeData = set_In;
}
