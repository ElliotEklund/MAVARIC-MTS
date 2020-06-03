#include "equilib_mvrpmd.hpp"

equilib_mvrpmd::equilib_mvrpmd(int my_id, int root_proc, int num_procs,
                               std::string root_path)
    :my_id(my_id),
     root_proc(root_proc),
     num_procs(num_procs),
     myRand(time(NULL) + my_id),
     helper(root_path,my_id,num_procs,root_proc)
{}
void equilib_mvrpmd::run(double nuc_ss, double x_ss, double p_ss,
                         unsigned long long num_steps, unsigned long long stride){
    
    //Declare vectors for monte carlo moves
    vector<double> Q(nuc_beads,0);
    matrix<double> x(elec_beads,num_states,0), p(elec_beads,num_states,0);

    //Initialize vectors
    gen_initQ(Q,nuc_beads,nuc_ss);
    gen_initElec(x,elec_beads,num_states,x_ss);
    gen_initElec(p,elec_beads,num_states,p_ss);
    
    //Over-write with saved vectors if requested
    if(readPSV){helper.read_PSV(nuc_beads, elec_beads, num_states, Q, x, p);}
    
    /* Assemble Hamiltonian and Estimator*/
    SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states);
    C_Matrix C(elec_beads,num_states);
    M_Matrix M(num_states,nuc_beads,beta/elec_beads);
    M_Matrix_MTS M_MTS(nuc_beads,elec_beads,num_states,M);
    Theta_MTS thetaMTS(num_states,elec_beads,C,M_MTS);
    dTheta_MTS_dBeta dthetaMTS_dBeta(nuc_beads,elec_beads,num_states,
                                      beta/elec_beads,C,M,M_MTS);

    MVRPMD_MTS_Hamiltonian H(beta/nuc_beads,V_spring,V0,G,thetaMTS);
    MVRPMD_MTS_Estimator Esti(nuc_beads,beta/nuc_beads,V_spring,V0,
                                   thetaMTS,dthetaMTS_dBeta);
    
    /* Initialize energy and estimator variables*/
    int esti_samples = 0;
    double estimator_total(0), sgn_total(0);
    double energy(0), estimator(0), sgnTheta(0);
    vector<double> estimator_t(num_steps/stride,0);
    
    energy = H.get_energy(Q,x,p);
    estimator = Esti.get_estimator();

    /* Over-write estimator with saved value if requested*/
    if (readData) {helper.read_MC_data(sgn_total, estimator_total);}

    /* Initialize nuclear stepping procedure */
    system_step nuc_stepper(my_id,num_procs,root_proc,nuc_beads,beta);
    nuc_stepper.set_nuc_ss(nuc_ss);
    nuc_stepper.set_hamiltonian(H);

    /* Initialize electronic stepping procedure*/
    elec_step elec_stepper(my_id,num_procs,root_proc,elec_beads,num_states,beta);
    elec_stepper.set_hamiltonian(H);
    elec_stepper.set_energy(energy);
    elec_stepper.set_ss(x_ss,p_ss);

    //r is a ratio that determines how often to sample each sub-system
    double r = double(nuc_beads)/ double(nuc_beads + num_states*elec_beads);

    for (int step=0; step<num_steps; step++) {
        if (myRand.doub() < r) {
            /* Sample Nuclear Coordinates*/
            nuc_stepper.step(energy,Q,x,p);
            energy = nuc_stepper.get_energy();
        }
        else{
            /* Sample Electronic Coordinates*/
            if (myRand.doub() >= 0.5) {
                /* Sample x*/
                elec_stepper.step_x(energy,Q,x,p);
                energy = elec_stepper.get_energy();
            }
            else{
                /* Sample p*/
                elec_stepper.step_p(energy,Q,x,p);
                energy = elec_stepper.get_energy();
            }
        }
        estimator = Esti.get_estimator(Q,x,p);
        sgnTheta = thetaMTS.get_signTheta();
        
        estimator_total += estimator;
        sgn_total += sgnTheta;
        
        if(step % stride == 0){
            estimator_t[esti_samples] = estimator_total/(sgn_total + 1);
            esti_samples += 1;
        }
    }

    /* Retrive acceptance ratio information */
    unsigned long long nuc_steps_total = nuc_stepper.get_steps_total();
    unsigned long long nuc_steps_accpt = nuc_stepper.get_steps_accepted();

    unsigned long long x_steps_accpt = elec_stepper.get_x_steps_accpt();
    unsigned long long x_steps_total = elec_stepper.get_x_steps_total();

    unsigned long long p_steps_accpt = elec_stepper.get_p_steps_accpt();
    unsigned long long p_steps_total = elec_stepper.get_p_steps_total();
    
    /* Print Monte Carlo simulation Data */
    helper.print_sys_accpt(nuc_steps_total, nuc_steps_accpt);
    helper.print_elec_accpt(x_steps_total, x_steps_accpt);
    helper.print_avg_energy(estimator_total, sgn_total);
    helper.write_estimator(estimator_t,stride);
        
    /* Write Monte Carlo information to file if requested*/
    if(writePSV){helper.write_PSV(nuc_beads, elec_beads, num_states, Q, x, p);}
    if (writeData) {helper.write_MC_data(sgn_total, estimator_total);}
}
void equilib_mvrpmd::gen_initQ(vector<double> &Q, int num_beads, double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = step_dist(myRand.doub(),step_size);
    }
}
void equilib_mvrpmd::gen_initElec(matrix<double> &v, int num_beads, int num_states,
                                   double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            v(bead,state) = step_dist(myRand.doub(),step_size);
        }
    }
}
void equilib_mvrpmd::initialize_system(int nuc_beads_IN,int elec_beadsIN,
                                        int num_statesIN,double massIN,double betaIN){
    nuc_beads = nuc_beads_IN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta = betaIN;
}
void equilib_mvrpmd::initialize_files(bool writePSV_IN, bool readPSV_IN,
                                      bool writeData_IN, bool readData_IN){
    writePSV = writePSV_IN;
    readPSV = readPSV_IN;
    writeData = writeData_IN;
    readData = readData_IN;
}
inline double equilib_mvrpmd::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
void equilib_mvrpmd::set_write_PSV(bool set_In){writePSV = set_In;}

void equilib_mvrpmd::set_read_PSV(bool set_In){readPSV = set_In;}

void equilib_mvrpmd::set_read_Data(bool set_In){readData = set_In;}

void equilib_mvrpmd::set_write_Data(bool set_In){writeData = set_In;}
