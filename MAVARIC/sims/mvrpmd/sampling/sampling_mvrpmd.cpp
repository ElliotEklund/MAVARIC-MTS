#include "sampling_mvrpmd.hpp"

sampling_mvrpmd::sampling_mvrpmd(int my_id, int root_proc, int num_procs)
    :my_id(my_id),
     root_proc(root_proc),
     num_procs(num_procs),
     myRand(time(NULL) + my_id),
     helper(my_id,num_procs,root_proc)
{}
void sampling_mvrpmd::run(double nuc_ss, double x_ss, double p_ss,
                          unsigned long long num_trajs,unsigned long long decorr){
    
    /* Convert num_trajs to num_trajs per processor*/
    if (num_trajs % num_procs != 0) {
        if (my_id == root_proc) {
            std::cout << "ERROR: num_trajs is not divisible by number of procs"
            << std::endl;
        }
    }
    else{
        num_trajs = num_trajs/num_procs;
    }
    
    /* Declare vectors for monte carlo moves*/
    vector<double> Q(nuc_beads,0);
    matrix<double> x(elec_beads,num_states,0), p(elec_beads,num_states,0);
    
    gen_initQ(Q,nuc_beads,nuc_ss);
    gen_initElec(x,elec_beads,num_states,x_ss);
    gen_initElec(p,elec_beads,num_states,p_ss);
    
    if (readPSV) {
        helper.read_PSV(nuc_beads,elec_beads,num_states,Q,x,p);
    }

    /* Declare vectors to store sampled trajectories*/
    vector<double> Q_trajs(num_trajs*nuc_beads,0);
    vector<double> P_trajs(num_trajs*nuc_beads,0);
    matrix<double> x_trajs(num_trajs*elec_beads,num_states,0);
    matrix<double> p_trajs(num_trajs*elec_beads,num_states,0);

    /* Assemble Hamiltonian*/
    SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states);
    C_Matrix C(elec_beads,num_states);
    M_Matrix M(num_states,elec_beads,beta/elec_beads);
    //M_Matrix_MTS M_MTS(nuc_beads,elec_beads,num_states,M);
    //Theta_MTS thetaMTS(num_states,elec_beads,C,M_MTS);
    theta_mixed theta(num_states,nuc_beads,elec_beads,C,M);

    //MVRPMD_MTS_Hamiltonian H(beta/nuc_beads,V_spring,V0,G,thetaMTS);
    mvrpmd_mixed_ham H(beta/nuc_beads,V_spring,V0,G,theta);
    
    double energy = H.get_energy(Q,x,p);
    double energy_prop = energy;
    
    /* Initialize nuclear stepping procedure */
    system_step nuc_stepper(my_id,num_procs,root_proc,nuc_beads,beta);
    nuc_stepper.set_nuc_ss(nuc_ss);
    nuc_stepper.set_hamiltonian(H);
    nuc_stepper.set_energy(energy);
    
    /* Initialize electronic stepping procedure*/
    elec_step elec_stepper(my_id,num_procs,root_proc,elec_beads,num_states,beta);
    elec_stepper.set_hamiltonian(H);
    elec_stepper.set_energy(energy);
    elec_stepper.set_ss(x_ss,p_ss);

    //r is a ratio that determines how often to sample each sub-system
    double r = double(nuc_beads)/ double(nuc_beads + num_states*elec_beads);

    /* Main Algorithm*/
    for (int traj=0; traj<num_trajs; traj++) {
        for (int step=0; step<decorr; step++) {
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
        }
        for (int bead=0; bead<nuc_beads; bead++) {
            Q_trajs(traj*nuc_beads+bead) = Q(bead);
        }
        
        for (int bead=0; bead<elec_beads; bead++) {
            for (int state=0; state<num_states; state++) {
                x_trajs(traj*elec_beads+bead,state) = x(bead,state);
                p_trajs(traj*elec_beads+bead,state) = p(bead,state);
            }
        }
    }
    
    //double stdev = sqrt(nuc_beads*mass/beta);
    double stdev = sqrt(mass/(beta*nuc_beads));

    //system momentum distribution from Gaussian(mu, sigma, seed)
    Normaldev_BM momentum(0, stdev, myRand.int32());

    for(int i=0; i<num_trajs*nuc_beads; i++){
        P_trajs(i) = momentum.dev();
    }

    save_trajs(Q_trajs,"Output/Trajectories/Q");
    save_trajs(P_trajs,"Output/Trajectories/P");
    
    save_trajs(x_trajs,num_trajs*elec_beads*num_states,num_trajs,
               "Output/Trajectories/xElec");
    
    save_trajs(p_trajs,num_trajs*elec_beads*num_states,num_trajs,
               "Output/Trajectories/pElec");
    
    unsigned long long nuc_steps_total = nuc_stepper.get_steps_total();
    unsigned long long nuc_steps_accpt = nuc_stepper.get_steps_accepted();
    
    unsigned long long x_steps_accpt = elec_stepper.get_x_steps_accpt();
    unsigned long long x_steps_total = elec_stepper.get_x_steps_total();
    
    unsigned long long p_steps_accpt = elec_stepper.get_p_steps_accpt();
    unsigned long long p_steps_total = elec_stepper.get_p_steps_total();
    
    helper.print_sys_accpt(nuc_steps_total,nuc_steps_accpt,"Nuclear");
    helper.print_sys_accpt(x_steps_total,x_steps_accpt,"x");
    helper.print_sys_accpt(p_steps_total,p_steps_accpt,"p");
}
void sampling_mvrpmd::gen_initQ(vector<double> &Q, int num_beads, double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = step_dist(myRand.doub(),step_size);
    }
}
void sampling_mvrpmd::gen_initElec(matrix<double> &v, int num_beads, int num_states,
                                   double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            v(bead,state) = step_dist(myRand.doub(),step_size);
        }
    }
}
inline double sampling_mvrpmd::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
void sampling_mvrpmd::initialize_system(int nuc_beads_IN,int elec_beadsIN,
                                        int num_statesIN,double massIN,double betaIN){
    nuc_beads = nuc_beads_IN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta = betaIN;
}
void sampling_mvrpmd::initialize_files(bool readPSVIN, std::string rootFolderIN){
    
    readPSV = readPSVIN;
    rootFolder = rootFolderIN;
    helper.set_root(rootFolderIN);
}
void sampling_mvrpmd::save_trajs(vector<double> &v,std::string name){

    std::string fileName = rootFolder + name;
    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    
    myWrap.write_vector(v,fileName);
}
void sampling_mvrpmd::save_trajs(matrix<double> &v,int size,int num_trajs,
                                 std::string name){

    vector<double> v_transform (size,0);
    int stride = 0;
    
    for (int traj=0; traj<num_trajs; traj++) {
        for (int bead=0; bead<elec_beads; bead++) {
            stride = traj*elec_beads*num_states + bead*num_states;
            for (int state=0; state<num_states; state++) {
                v_transform(stride + state) = v(traj*elec_beads+bead,state);
            }
        }
    }
    
    std::string fileName = rootFolder + name;
    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    myWrap.write_vector(v_transform,fileName);
}
