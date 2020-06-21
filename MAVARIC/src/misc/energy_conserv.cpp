#include "energy_conserv.hpp"

energy_conserv::energy_conserv(int my_id, int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc)
{
    is_sys_set = false;
    is_time_set = false;
}
void energy_conserv::compute(unsigned long long num_trajs_global,
                             unsigned long long num_trajs_local, double tol,
                             int energy_stride,
                             std::string input_dir, std::string output_dir){
    
    int num_steps = floor(total_time/dt);
    
    C_Matrix C(elec_beads, num_states,alpha);
    M_Matrix M(num_states, elec_beads, beta_elec_beads);
    dM_Matrix_dQ dMdQ(elec_beads, num_states, beta_elec_beads, M);
    theta_mixed theta(num_states,nuc_beads,elec_beads,C,M);
    theta_mixed_dQ theta_dQ(num_states,nuc_beads,elec_beads,C,M,dMdQ);
    theta_mixed_dElec theta_dElec(num_states,elec_beads,alpha,C,M);
    
    mvrpmd_mixed_forces F(nuc_beads, elec_beads, num_states, mass,
                          beta_nuc_beads, alpha, theta, theta_dQ, theta_dElec);
    
    ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);
    
    SpringEnergy V_spring(nuc_beads,mass,beta_nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states,alpha);
    mvrpmd_mixed_ham H(beta_nuc_beads,V_spring,V0,G,theta);
    
    vector<double> Q_traj (nuc_beads);
    vector<double> P_traj (nuc_beads);
    matrix<double> x_traj (elec_beads,num_states);
    matrix<double> p_traj (elec_beads,num_states);
    
    vector<vector<double> > Q(num_trajs_local,zero_vector<double>(nuc_beads));
    vector<vector<double> > P(num_trajs_local,zero_vector<double>(nuc_beads));
    vector<matrix<double> > x(num_trajs_local,zero_matrix<double>
                              (nuc_beads,num_states));
    vector<matrix<double> > p(num_trajs_local,zero_matrix<double>
                              (nuc_beads,num_states));
    
    std::string Q_file = input_dir + "Q";
    std::string P_file = input_dir + "P";
    std::string x_file = input_dir + "xElec";
    std::string p_file = input_dir + "pElec";
    
    Q = get_trajs_reformat(Q_file,num_trajs_global*nuc_beads,
                           num_trajs_local*nuc_beads,my_id,num_procs,
                           root_proc,num_trajs_local,nuc_beads);

    P = get_trajs_reformat(P_file,num_trajs_global*nuc_beads,
                           num_trajs_local*nuc_beads,my_id,num_procs,
                           root_proc,num_trajs_local,nuc_beads);

    x = get_trajs_reformat(x_file,num_trajs_global*elec_beads*num_states,
                           num_trajs_local*elec_beads*num_states,my_id,num_procs,
                           root_proc,num_trajs_local,elec_beads,num_states);

    p = get_trajs_reformat(p_file,num_trajs_global*elec_beads*num_states,
                           num_trajs_local*elec_beads*num_states,my_id,num_procs,
                           root_proc,num_trajs_local,elec_beads,num_states);

    bool broken = false; //true if current trajectory is broken
    int step = 0; //current step of current trajectory simulation
    double energy_init = 0; //initial energy of a given trajectory
    double energy_t = 0; //energy of a given trajectory at time t
    double badness = 0; // abs((energy_init-energy_t)/energy_init)
    /*  holds integers corresponding to which trajectories have broken*/
    std::list <int> broken_trajectories;
    
    int ten_p = floor(num_trajs_local/10.0);
    bool get_prog;
    int get_prog_int = 1;
    
    if ((my_id==root_proc) && (ten_p == 0)) {
        get_prog_int = 0;
        MPI_Bcast(&get_prog_int,1,MPI_INT,root_proc,MPI_COMM_WORLD);
    }
    
    if (get_prog_int==0) {
        get_prog = false;
    }
    
    std::ofstream progress;
    if (my_id==root_proc) {
        std::string file_name = output_dir + "conserv_prog";
        progress.open(file_name.c_str());
        if (!progress.is_open()) {
            std::cout << "ERROR: Could not open " << file_name << std::endl;
        }
        if (!get_prog) {
            std::cout << "Warning: energy conservation progress is not collected"
            "because too few trajectories are used." << std::endl;
        }
    }
    
    
    for (int traj=0; traj<num_trajs_local; traj++){

        /* Load new trajecty*/
        Q_traj = Q(traj);
        P_traj = P(traj);
        x_traj = x(traj);
        p_traj = p(traj);
        
        broken = false; //reset broken trajectory to false

        step = 0; //reset step to zero
        myABM.initialize_rk4(Q_traj, P_traj, x_traj, p_traj);
        energy_init = H.get_energy_dyn(mass,Q_traj,P_traj,x_traj,p_traj);
        
        while (step<num_steps && !broken){

            myABM.take_step(Q_traj, P_traj, x_traj, p_traj);
            step += 1;

            if (step % energy_stride == 0) {
                /* Test if trajectory has broken*/

                energy_t = H.get_energy_dyn(mass,Q_traj,P_traj,x_traj,p_traj);
                badness = abs(energy_init - energy_t)/energy_t;

                if (badness > tol || is_NaN(energy_t)) {
                    /* Trajectory is broken*/
                    broken = true;
                    broken_trajectories.push_front(traj);
                }
            }
        }
        
        if (get_prog) {
            if ((traj % ten_p == 0) && (my_id == root_proc)) {
                progress << traj/ten_p << "%" << std::endl;
            }
        }
    }
    
    if (my_id==root_proc) {
        progress.close();
    }

    write_broken(broken_trajectories,num_trajs_local,output_dir);
    write_report(broken_trajectories,num_trajs_global,tol,output_dir);
}
void energy_conserv::write_broken(std::list<int> broken,
                                  unsigned long long num_trajs_local,
                                  std::string output_dir){
    
    int num_broken = broken.size();
    vector<int> broken_array (num_broken,0);
    std::list<int>::iterator it;
    int ii = 0;
    int traj_off_set = my_id * num_trajs_local; //correct trajectory for offset
    
    /* Convert list to array. Cannot get list to work with MPI :(*/
    for (it = broken.begin(); it != broken.end(); it++) {
        broken_array(ii) = *it + traj_off_set;
        ii++;
    }
    
    int num_broke_glo = 0;
    MPI_Reduce(&num_broken,&num_broke_glo,1,MPI_INT,MPI_SUM,
               root_proc,MPI_COMM_WORLD);
    
    int all_sizes [num_procs]; //collection of v_size across all procs
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Gather(&num_broken,1,MPI_INT,&all_sizes,1,
               MPI_INT,root_proc,MPI_COMM_WORLD);

    vector<int> v_final(num_broke_glo,0);
    int displs[num_procs]; //buffer space for a given proc
    displs[0] = 0;

    for (int i=1; i<num_procs; i++) {
        displs[i] = displs[i-1] + all_sizes[i-1];
    }

    MPI_Gatherv(&broken_array(0),num_broken, MPI_INT,
                &v_final[0], all_sizes, displs, MPI_INT,
                root_proc, MPI_COMM_WORLD);

    if (my_id==root_proc) {
        std::ofstream myFile;
        std::string fileName = output_dir + "broken";
        myFile.open(fileName.c_str());
        if (!myFile.is_open()) {
            std::cout << "ERROR: Could not open " << fileName << std::endl;
        }
        
        std::vector<int> broke_sort(num_broke_glo);
        for (int i=0; i<num_broke_glo; i++) {
            broke_sort[i] = v_final[i];
        }
        std::sort(broke_sort.begin(),broke_sort.end());
        for (int i=0; i<num_broke_glo; i++) {
            myFile << broke_sort[i] << std::endl;
        }
    }
}
void energy_conserv::write_report(std::list <int> broken_trajectories,
                                  unsigned long long num_trajs_global,double tol,
                                  std::string output_dir){
    
    int num_broke_loc = broken_trajectories.size();
    int num_broke_glo = 0;
    
    MPI_Reduce(&num_broke_loc,&num_broke_glo,1,MPI_INT,MPI_SUM,root_proc,MPI_COMM_WORLD);
    if (my_id == root_proc){
        std::cout << "Percent broken: " << 100 * num_broke_glo/double(num_trajs_global)
                  << std::endl;
        
        std::ofstream my_stream;
        std::string file_name = output_dir + "conserv_report";
        my_stream.open(file_name.c_str());
        
        if(!my_stream.is_open()) {
            std::cout << "ERROR: Could not open file " << file_name << std::endl;
        }
        
        my_stream << "dt:" << dt << std::endl;
        my_stream << "run time:" << total_time << std::endl;
        my_stream << "tolerance:" << tol << std::endl;
        my_stream << "percent broken:" << 100 * num_broke_glo/double(num_trajs_global)
                    << std::endl;
        
        my_stream.close();
    }
}
void energy_conserv::set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                                 double massIN,double betaIN, double beta_nuc_beadsIN,
                                 double beta_elec_beadsIN, double alphaIN){
    nuc_beads = nuc_beadsIN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta = betaIN;
    beta_nuc_beads = beta_nuc_beadsIN;
    beta_elec_beads = beta_nuc_beadsIN;
    alpha = alphaIN;
    is_sys_set = true;
}
void energy_conserv::set_time(double dtIN, double total_timeIN){
    dt = dtIN;
    total_time = total_timeIN;
    is_time_set = true;
}
