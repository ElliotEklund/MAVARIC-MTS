#include "sampling_mvrpmd.hpp"

sampling_mvrpmd::sampling_mvrpmd(int my_id, int root_proc, int num_procs)
    :my_id(my_id),
     root_proc(root_proc),
     num_procs(num_procs),
     myRand(time(NULL) + my_id),
     helper(my_id,num_procs,root_proc)
{}

void sampling_mvrpmd::run(double nuc_ss, double elec_ss,
                          unsigned long long num_trajs,
                          unsigned long long decorr){
    
    num_trajs = num_trajs/num_procs;
    
    /* Declare vectors for monte carlo moves*/
    vector<double> Q(nuc_beads,0), Q_prop(nuc_beads,0);
    matrix<double> x(elec_beads,num_states,0), x_prop(elec_beads,num_states,0);
    matrix<double> p(elec_beads,num_states,0), p_prop(elec_beads,num_states,0);
    
    /* Declare vectors to store sampled trajectories*/
    vector<double> Q_trajs(num_trajs*nuc_beads,0);
    vector<double> P_trajs(num_trajs*nuc_beads,0);
    matrix<double> x_trajs(num_trajs*elec_beads,num_states,0);
    matrix<double> p_trajs(num_trajs*elec_beads,num_states,0);
    
    gen_initQ(Q,nuc_beads,nuc_ss);
    gen_initElec(x,elec_beads,num_states,elec_ss);
    gen_initElec(p,elec_beads,num_states,elec_ss);

    if (readPSV) {
        helper.read_PSV(nuc_beads,elec_beads,num_states,Q,x,p);
    }
    
    Q_prop = Q;
    x_prop = x;
    p_prop = p;

    /* Assemble Hamiltonian*/
    SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states);
    C_Matrix C(elec_beads,num_states);
    M_Matrix M(num_states,nuc_beads,beta/elec_beads);
    M_Matrix_MTS M_MTS(nuc_beads,elec_beads,num_states,M);
    Theta_MTS thetaMTS(num_states,elec_beads,C,M_MTS);
    MVRPMD_MTS_Hamiltonian H(beta/nuc_beads,V_spring,V0,G,thetaMTS);
    
    double energy = H.get_energy(Q,x,p);
    double energy_prop = energy;
    
    int mcMove = 0;
    bool accept_move = false;

    double r = double(nuc_beads)/double(elec_beads);
    if (r==1) {
        r = 0.5;
    }
    
    unsigned long long nuc_steps_accpt(0), nuc_steps_total(0);
    unsigned long long elec_steps_accpt(0), elec_steps_total(0);

    /* Main Algorithm*/
    for (int traj=0; traj<num_trajs; traj++) {
        for (int step=0; step<decorr; step++) {
            if (myRand.doub() > r) {
                /* Sample Nuclear Coordinates*/
                
                for (int bead=0; bead<nuc_beads; bead++) {
                    mcMove = rand_bead(myRand.int64(),nuc_beads);
                    Q_prop(mcMove) = Q(mcMove) + step_dist(myRand.doub(),nuc_ss);
                    energy_prop = H.get_energy(Q_prop,x,p);
                    accept_move = check_move(energy,energy_prop);
                    
                    if (accept_move) {
                        Q(mcMove) = Q_prop(mcMove);
                        energy = energy_prop;
                        nuc_steps_accpt += 1;
                    }
                    else{
                        Q_prop(mcMove) = Q(mcMove);
                    }
                    nuc_steps_total += 1;
                }
            }
            else{
                /* Sample Electronic Coordinates*/
                for (int bead=0; bead<elec_beads; bead++) {
                    mcMove = rand_bead(myRand.int64(),elec_beads);
                    for (int state=0; state<num_states; state++) {
                        x_prop(mcMove,state) = x(mcMove,state) +
                        step_dist(myRand.doub(),elec_ss);
                        p_prop(mcMove,state) = p(mcMove,state) +
                        step_dist(myRand.doub(),elec_ss);
                    }
                    
                    energy_prop = H.get_energy(Q,x_prop,p_prop);
                    accept_move = check_move(energy,energy_prop);
                    
                    if (accept_move) {
                        x = x_prop;
                        p = p_prop;
                        energy = energy_prop;
                        elec_steps_accpt += 1;
                    }
                    else{
                        x_prop = x;
                        p_prop = p;
                    }
                    elec_steps_total += 1;
                }
            }
        }
        for (int bead=0; bead<nuc_beads; bead++) {
            Q_trajs(traj*nuc_beads+bead) = Q(bead);
            for (int state=0; state<num_states; state++) {
                x_trajs(traj*elec_beads+bead,state) = x(bead,state);
                p_trajs(traj*elec_beads+bead,state) = p(bead,state);
            }
        }
    }
    
    //double stdev = sqrt(nuc_beads*mass/beta);
    double stdev = sqrt(mass/(beta*nuc_beads));

    //system momentum distribution from Gaussian(mu, sigma, seed)
    Normaldev_BM momentum(0, stdev, rand());

    for(int i=0; i<num_trajs*nuc_beads; i++){
        P_trajs(i) = momentum.dev();
    }

    save_trajs(Q_trajs,"Output/Trajectories/Q");
    save_trajs(P_trajs,"Output/Trajectories/P");
    save_trajs(x_trajs,num_trajs*elec_beads*num_states,num_trajs,
               "Output/Trajectories/xElec");
    save_trajs(p_trajs,num_trajs*elec_beads*num_states,num_trajs,
               "Output/Trajectories/pElec");
    
    helper.print_sys_accpt(nuc_steps_total,nuc_steps_accpt,"Nuclear");
    helper.print_sys_accpt(elec_steps_total,elec_steps_accpt,"Electronic");
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

inline int sampling_mvrpmd::rand_bead(const Ullong rn, int num_beads){
    return rn % num_beads;
}
bool sampling_mvrpmd::check_move(double energy, double energy_prop){

    double delta_energy = energy_prop - energy;
    double accept_move = false;
    /* Accept new system moves if energ_prop < energy*/
    if(delta_energy < 0){
        accept_move = true;
    }
    /* Accept new system moves if inequality is met*/
    //else if (myRand.doub() <= exp(-beta * delta_energy/nuc_beads)){
    else if (myRand.doub() <= exp(-beta * delta_energy)){

        accept_move = true;
    }
    return accept_move;
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

void sampling_mvrpmd::save_trajs(const vector<double> &v,std::string name){

    std::string fileName = rootFolder + name;
    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    
    myWrap.write_vector(v,fileName);
}

void sampling_mvrpmd::save_trajs(const matrix<double> &v,int size,int num_trajs,
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
