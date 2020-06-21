#include "sampling_mvrpmd.hpp"

sampling_mvrpmd::sampling_mvrpmd(int my_id, int root_proc, int num_procs)
    :my_id(my_id),
     root_proc(root_proc),
     num_procs(num_procs),
     myRand(time(NULL) + my_id),
     helper(my_id,num_procs,root_proc)
{
    sys_set = false;
    files_set = false;
}
int sampling_mvrpmd::run(double nuc_ss, double x_ss, double p_ss,
                          unsigned long long num_trajs,unsigned long long decorr){
    
    if (!sys_set || !files_set) {
        if (my_id==root_proc) {
            std::cout << "ERROR: sampling mvrpmd variables not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
        return -1;
    }
    
    unsigned long long num_trajs_local = 0;
    num_trajs_local = get_trajs_local(num_trajs);
    
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
    vector<double> Q_trajs(num_trajs_local*nuc_beads,0);
    vector<double> P_trajs(num_trajs_local*nuc_beads,0);
    matrix<double> x_trajs(num_trajs_local*elec_beads,num_states,0);
    matrix<double> p_trajs(num_trajs_local*elec_beads,num_states,0);

    /* Assemble Hamiltonian*/
    SpringEnergy V_spring(nuc_beads,mass,beta/nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states,alpha);
    C_Matrix C(elec_beads,num_states,alpha);
    M_Matrix M(num_states,elec_beads,beta/elec_beads);
    theta_mixed theta(num_states,nuc_beads,elec_beads,C,M);
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

    int ten_p = floor(num_trajs_local/10.0); //ten percent of trajectories
    bool get_report = true;
    
    if (ten_p == 0) {
        if (my_id==root_proc) {
            std::cout << "WARNING: Sampling progress report will not be generated "
            "because num_trajs_local < 10" << std::endl;
            
            get_report = false;
        }
    }
    
    std::ofstream progress;
    if (my_id == root_proc) {
        std::string my_file = rootFolder + "Output/samp_progress";
        progress.open(my_file.c_str());
        if (!progress.is_open()) {
            std::cout << "ERROR: could not open file " << my_file << std::endl;
        }
    }

    /* Main Algorithm*/
    for (int traj=0; traj<num_trajs_local; traj++) {
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
        if (my_id == root_proc && get_report) {
            if (traj % ten_p == 0) {
                progress <<  100 * double(traj)/num_trajs_local << "%" << std::endl;
            }
        }
    }
    
    if (my_id==root_proc) {progress.close();}

    //double stdev = sqrt(nuc_beads*mass/beta);
    double stdev = sqrt(mass/(beta*nuc_beads));
    //system momentum distribution from Gaussian(mu, sigma, seed)
    Normaldev_BM momentum(0, stdev, myRand.int32());

    for(int i=0; i<num_trajs_local*nuc_beads; i++){
        P_trajs(i) = momentum.dev();
    }
    
    if (saveTrajs) {
        if (!contains_NaN(Q_trajs)) {
            save_trajs(Q_trajs,"Output/Trajectories/Q",nuc_beads,elec_beads,num_states,
                       beta,num_trajs,decorr);
        }
        else{
            if (my_id==root_proc) {
                std::cout << "Bad value found in sampled Q trajectories." << std::endl;
            }
        }
        if (!contains_NaN(P_trajs)) {
            save_trajs(P_trajs,"Output/Trajectories/P",nuc_beads,elec_beads,num_states,
                       beta,num_trajs,decorr);
        }
        else{
            if (my_id==root_proc) {
                std::cout << "Bad value found in sampled P trajectories." << std::endl;
            }
        }
        if (!contains_NaN(x_trajs)) {
            save_trajs(x_trajs,num_trajs_local*elec_beads*num_states,num_trajs_local,
                       "Output/Trajectories/xElec",nuc_beads,elec_beads,num_states,
                       beta,num_trajs,decorr);
        }
        else{
            if (my_id==root_proc) {
                std::cout << "Bad value found in sampled x trajectories." << std::endl;
            }
        }
        if (!contains_NaN(p_trajs)) {
            save_trajs(p_trajs,num_trajs_local*elec_beads*num_states,num_trajs_local,
                       "Output/Trajectories/pElec",nuc_beads,elec_beads,num_states,
                       beta,num_trajs,decorr);
        }
        else{
            if (my_id==root_proc) {
                std::cout << "Bad value found in sampled p trajectories." << std::endl;
            }
        }
    }
    
    unsigned long long nuc_steps_total = nuc_stepper.get_steps_total();
    unsigned long long nuc_steps_accpt = nuc_stepper.get_steps_accepted();
    unsigned long long x_steps_accpt = elec_stepper.get_x_steps_accpt();
    unsigned long long x_steps_total = elec_stepper.get_x_steps_total();
    unsigned long long p_steps_accpt = elec_stepper.get_p_steps_accpt();
    unsigned long long p_steps_total = elec_stepper.get_p_steps_total();

    helper.print_sys_accpt(nuc_steps_total,nuc_steps_accpt,"Nuclear");
    helper.print_sys_accpt(x_steps_total,x_steps_accpt,"x");
    helper.print_sys_accpt(p_steps_total,p_steps_accpt,"p");
    
    return 0;
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
                                        int num_statesIN,double massIN,double betaIN,
                                        double alphaIN){
    nuc_beads = nuc_beads_IN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta = betaIN;
    alpha = alphaIN;
    sys_set = true;
}
void sampling_mvrpmd::initialize_files(bool readPSVIN,bool saveTrajsIN,
                                       std::string rootFolderIN){
    
    readPSV = readPSVIN;
    rootFolder = rootFolderIN;
    helper.set_root(rootFolderIN);
    rootFolder = rootFolderIN;
    saveTrajs = saveTrajsIN;
    files_set = true;
}
void sampling_mvrpmd::save_trajs(vector<double> &v,std::string name, int nuc_beadsIN,
                                 int elec_beadsIN, int num_statesIN,double betaIN,
                                 unsigned long long num_trajs_totalIN,double decorrIN){

    std::string fileName = rootFolder + name;
    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    
    myWrap.write_vector(v,fileName,nuc_beadsIN,elec_beadsIN,num_statesIN,betaIN,
                        num_trajs_totalIN,decorrIN);
}
void sampling_mvrpmd::save_trajs(matrix<double> &v,int size,
                                 unsigned long long num_trajsIN,std::string name,
                                 int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                                 double betaIN,unsigned long long num_trajs_totalIN,
                                 double decorrIN){

    vector<double> v_transform (size,0);
    int stride = 0;
    
    for (unsigned long long traj=0; traj<num_trajsIN; traj++) {
        for (unsigned long long bead=0; bead<elec_beads; bead++) {
            stride = traj*elec_beads*num_states + bead*num_states;
            for (unsigned long long state=0; state<num_states; state++) {
                v_transform(stride + state) = v(traj*elec_beads+bead,state);
            }
        }
    }

    std::string fileName = rootFolder + name;
    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    myWrap.write_vector(v_transform,fileName,nuc_beadsIN,elec_beadsIN,num_statesIN,
                        betaIN,num_trajs_totalIN,decorrIN);
}
unsigned long long sampling_mvrpmd::get_trajs_local(unsigned long long num_trajs_totalIN){
    
    unsigned long long num_trajs_local = 0; //number of local trajectories
    
    /* Check for more trajs than procs*/
    if (num_trajs_totalIN < num_procs) {
        if (my_id == root_proc) {
            std::cout << "ERROR: num_trajs is small than num_procs.!" << std::endl;
        }
        return 0;
    }
    else{
        unsigned long long d = num_trajs_totalIN/num_procs;
        unsigned long long r = num_trajs_totalIN%num_procs;
        
        if (my_id < num_procs - r) {
            num_trajs_local = d;
        }
        else{
            num_trajs_local = d + 1;
        }
        return num_trajs_local;
    }
}
