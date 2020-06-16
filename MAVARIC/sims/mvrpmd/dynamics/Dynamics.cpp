#include "Dynamics.hpp"

Dynamics::Dynamics(int num_procs, int my_id, int root_proc)

    :num_procs(num_procs),
     my_id(my_id),
     root_proc(root_proc),
        
     dt_is_set(false),
     total_time_is_set(false),
    
//     Q(num_trajs_local,zero_vector<double>(nuc_beads)),
//     P(num_trajs_local,zero_vector<double>(nuc_beads)),
//     x(num_trajs_local,zero_matrix<double>(elec_beads,num_states)),
//     p(num_trajs_local,zero_matrix<double>(elec_beads,num_states)),
//
//     C(elec_beads, num_states,alpha),
//     M(num_states, elec_beads, beta_elec_beads),
//
//     M_MTS(nuc_beads, elec_beads, num_states, M),
//     dMdQ(elec_beads, num_states, beta_elec_beads, M),
//
//    theta(num_states,nuc_beads,elec_beads,C,M),
//    theta_dQ(num_states,nuc_beads,elec_beads,C,M,dMdQ),
//    theta_dElec(num_states,elec_beads,alpha,C,M),
//
//    F(nuc_beads, elec_beads, num_states, mass,
//      beta_nuc_beads, alpha, theta, theta_dQ, theta_dElec)
{

    /* Global vectors hold all trajectories on the root processor.
     * These trajectories are distributed to the local vector.*/
    vector<double> Q_global = zero_vector<double> (num_trajs*nuc_beads);
    vector<double> P_global = zero_vector<double> (num_trajs*nuc_beads);
    vector<double> x_global = zero_vector<double> (num_trajs*elec_beads*num_states);
    vector<double> p_global = zero_vector<double> (num_trajs*elec_beads*num_states);
    
    /* Local vectors hold trajectories that each local processor will
     * work with.*/
    vector<double> Q_local = zero_vector<double> (num_trajs_local*nuc_beads);
    vector<double> P_local = zero_vector<double> (num_trajs_local*nuc_beads);
    vector<double> x_local = zero_vector<double> (num_trajs_local*elec_beads*num_states);
    vector<double> p_local = zero_vector<double> (num_trajs_local*elec_beads*num_states);
    
    std::string fileName = root + "Output/Trajectories/";
    
    /* Read in trajectories to global PSV vectors  */
    if(my_id == root_proc){
        load_var(Q_global,"Q",fileName);
        load_var(P_global,"P",fileName);
        load_var(x_global,"xElec",fileName);
        load_var(p_global,"pElec",fileName);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* Distribute global trajectories across all processors. */
    MPI_Scatter(&Q_global[0],num_trajs_local*nuc_beads,MPI_DOUBLE,&Q_local[0],
                num_trajs_local*nuc_beads,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);
    
    MPI_Scatter(&P_global[0],num_trajs_local*nuc_beads,MPI_DOUBLE,&P_local[0],
                num_trajs_local*nuc_beads,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);
    
    MPI_Scatter(&x_global[0],num_trajs_local*elec_beads*num_states,MPI_DOUBLE,&x_local[0],
                num_trajs_local*elec_beads*num_states,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);
    
    MPI_Scatter(&p_global[0],num_trajs_local*elec_beads*num_states,MPI_DOUBLE,&p_local[0],
                num_trajs_local*elec_beads*num_states,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);
    
    /* Trajectories are now distributed across all processors, but still need to be
     * correctly formated. */
    format_array(Q,Q_local);
    format_array(P,P_local);
    format_array(x,x_local);
    format_array(p,p_local);
}
void Dynamics::compute_initPAC(int interval){
    
    init_PAC my_init_PAC(num_procs,my_id,root_proc,nuc_beads,elec_beads,
                         num_states,num_trajs,num_trajs_local,theta);
    
    my_init_PAC.set_interval(interval);
    my_init_PAC.set_vectors(Q,x,p);
    my_init_PAC.compute(root);
}
void Dynamics::energ_conserv(double tol, int energy_stride){
    
    if (!total_time) {
        std::cout << "ERROR: total_time must be set before a simulation can be run." << std::endl;
    }

    ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);

    SpringEnergy V_spring(nuc_beads,mass,beta_nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states,alpha);
    mvrpmd_mixed_ham H(beta_nuc_beads,V_spring,V0,G,theta);

    vector<double> Q_traj (nuc_beads);
    vector<double> P_traj (nuc_beads);
    matrix<double> x_traj (elec_beads,num_states);
    matrix<double> p_traj (elec_beads,num_states);

    bool broken = false; //true if current trajectory is broken
    int step = 0; //current step of current trajectory simulation
    double energy_init = 0; //initial energy of a given trajectory
    double energy_t = 0; //energy of a given trajectory at time t
    double badness = 0; // abs((energy_init-energy_t)/energy_init)
    std::list <int> broken_trajectories; //holds integers corresponding to which trajectories have broken
    
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

                if (badness > tol) {
                    /* Trajectory is broken*/
                    broken = true;
                    broken_trajectories.push_front(traj);
                }
            }
        }
    }
    
    write_broken(broken_trajectories,root);

    int num_broke_loc = broken_trajectories.size();
    int num_broke_glo = 0;

    MPI_Reduce(&num_broke_loc,&num_broke_glo,1,MPI_INT,MPI_SUM,root_proc,MPI_COMM_WORLD);
    if (my_id == root_proc){
        std::cout << "Percent broken: " << 100 * num_broke_glo/double(num_trajs) << std::endl;
        
        std::ofstream my_stream;
        std::string file_name = root + "Output/conserv_report";
        my_stream.open(file_name.c_str());
        
        if(!my_stream.is_open()) {
            std::cout << "ERROR: Could not open file " << file_name << std::endl;
        }
        
        my_stream << "dt:" << dt << std::endl;
        my_stream << "run time:" << total_time << std::endl;
        my_stream << "tolerance:" << tol << std::endl;
        my_stream << "percent broken:" << 100 * num_broke_glo/double(num_trajs) << std::endl;
        
        my_stream.close();
    }
}
void Dynamics::PopAC(bool pac, int pac_stride, bool bp, int bp_stride,
                     bool sp, int sp_stride, bool wp, int wp_stride){
    
    if (!total_time) {
        std::cout << "ERROR: total_time must be set before a simulation can be run." << std::endl;
    }
    
    C_Matrix C(elec_beads, num_states,alpha),
    M_Matrix M(num_states, elec_beads, beta_elec_beads),
    dM_Matrix_dQ dMdQ(elec_beads, num_states, beta_elec_beads, M),
    theta_mixed theta(num_states,nuc_beads,elec_beads,C,M),
    theta_mixed_dQ theta_dQ(num_states,nuc_beads,elec_beads,C,M,dMdQ),
    theta_mixed_dElec theta_dElec(num_states,elec_beads,alpha,C,M),
    
    F(nuc_beads, elec_beads, num_states, mass,
      beta_nuc_beads, alpha, theta, theta_dQ, theta_dElec)
    
    ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);

    vector<double> Q_traj (nuc_beads), P_traj (nuc_beads);
    matrix<double> x_traj (elec_beads,num_states), p_traj (elec_beads,num_states);

    aggregate myAggregator;
    pop_estimators myPops(elec_beads,num_states,alpha);

    vector<double> pac_v0, bp_v0, sp_v0, wp_v0;
    vector<double> pac_v, bp_v, sp_v, wp_v;

    if (pac){
        myAggregator.add_calc("position",1,num_steps/pac_stride);
        pac_v.resize(1);
        pac_v0.resize(1);
    }
    if (bp){
        myAggregator.add_calc("boltzman",num_states,num_steps/bp_stride);
        bp_v.resize(num_states);
        bp_v0.resize(num_states);
    }
    if (sp){
        myAggregator.add_calc("semi_classic",num_states,num_steps/sp_stride);
        sp_v.resize(num_states);
        sp_v0.resize(num_states);
    }
    if (wp){
        myAggregator.add_calc("wigner",num_states,num_steps/wp_stride);
        wp_v.resize(num_states);
        wp_v0.resize(num_states);
    }

    double sgnTheta = 0; //sign of Theta for a trajectory
    
    int ten_p = floor(double(num_trajs_local)/10.0); //ten percent of trajs
    std::ofstream progress;
    if (my_id == root_proc) {
        std::string file_name = root + "Output/dyn_progress";
        progress.open(file_name.c_str());
    }
    
    for (int traj=0; traj<num_trajs_local; traj++) {
        /* Load new trajecty*/
        Q_traj = Q(traj);
        P_traj = P(traj);
        x_traj = x(traj);
        p_traj = p(traj);

        myABM.initialize_rk4(Q_traj, P_traj, x_traj, p_traj);
        sgnTheta = F.get_sign(Q_traj,x_traj,p_traj);

        if (pac){
            pac_v0(0) = compute_centroid(Q_traj);
            myAggregator.collect("position",0,pac_v0,pac_v0,sgnTheta);
        }
        if (bp){
            bp_v0 = myPops.boltz(theta.get_gamm());
            myAggregator.collect("boltzman",0,bp_v0,bp_v0,sgnTheta);
        }
        if (sp){
            //sp_v0 = myPops.sc(x_traj,p_traj);
            /* Use boltzman for time zero */
            sp_v0 = myPops.boltz(theta.get_gamm());
            sp_v = myPops.sc(x_traj,p_traj);
            myAggregator.collect("semi_classic",0,sp_v0,sp_v,sgnTheta);
        }
        if (wp){
            wp_v0 = myPops.wigner(x_traj,p_traj);
            myAggregator.collect("wigner",0,wp_v0,wp_v0,sgnTheta);
        }

        for (int step=1; step<num_steps; step++) {

            myABM.take_step(Q_traj, P_traj, x_traj, p_traj);

            if (pac){
                if (step % pac_stride == 0) {
                    pac_v(0) = compute_centroid(Q_traj);
                    myAggregator.collect("position",step/sp_stride,pac_v0,pac_v,sgnTheta);
                }
            }
            if (bp){
                if (step % bp_stride == 0) {
                    bp_v = myPops.boltz(theta.get_gamm());
                    myAggregator.collect("boltzman",step/bp_stride,bp_v0,bp_v,sgnTheta);
                }
            }
            if (sp){
                if (step % sp_stride == 0) {
                    sp_v = myPops.sc(x_traj,p_traj);
                    myAggregator.collect("semi_classic",step/sp_stride,sp_v0,sp_v,sgnTheta);
                }
            }
            if (wp){
                if (step % wp_stride == 0) {
                    wp_v = myPops.wigner(x_traj,p_traj);
                    myAggregator.collect("wigner",step/wp_stride,wp_v0,wp_v,sgnTheta);
                }
            }
        }
        
        if (traj % ten_p == 0) {
            if (my_id == root_proc) {
                progress << 100 * (double) traj /num_trajs_local << "%" << std::endl;
            }
            /* Save progress */
            std::string fileName = root + "Output/";
            myAggregator.merge_collections(root_proc,my_id,fileName,dt,sp_stride,
                                           traj*num_procs);
        }
    }
    
    if (my_id==root_proc) {
        progress.close();
    }

    std::string fileName = root + "Output/";
    myAggregator.merge_collections(root_proc,my_id,fileName,dt,sp_stride,num_trajs);
}
void Dynamics::write_broken(std::list<int> broken,std::string file_root){
    
    std::ostringstream quick_convert;
    quick_convert << my_id;
    
    std::string file_name = file_root + "Output/broken" + quick_convert.str();
    
    std::ofstream myFile;
    myFile.open(file_name.c_str());
    
    if (!myFile.is_open()) {
        std::cout << "Could not open " << file_name << std::endl;
    }
    
    int traj_off_set = my_id * num_trajs_local; //correct trajectory for offset
    int broken_array [broken.size()];
    
    std::list<int>::iterator it;
    int ii = 0;

    /* Convert list to array. Cannot get list to work with MPI :(*/
    for (it = broken.begin(); it != broken.end(); it++) {
        broken_array[ii] = *it + traj_off_set;
        ii++;
    }
    
    /* The first line of each file is the number of broken trajectories */
    if (broken.size() > 0) {
        myFile << broken.size() << std::endl;
        for (int i=0; i<broken.size(); i++) {
            myFile << broken_array[i] << std::endl;
        }
    }
    else {
        myFile << 0 << std::endl;
    }
    
    myFile.close();
}

double Dynamics::compute_centroid(const vector<double> &Q){
    double centroid = sum(Q);
    return centroid/nuc_beads;
}
void Dynamics::load_var(vector<double> &X, std::string var, std::string root_path){
      
    int num = X.size(); //size of vector

    /* Add root to specific file indicator*/
    std::string file_name =  root_path + var;
    std::ifstream myFile;
    myFile.open(file_name.c_str());
        
    if(!myFile.is_open()) {
        std::cout << "Could not open file " << var << std::endl;
    }
 
    /* Skip past file comment lines*/
    std::string line; //current line
    int num_lines = 0; //number of comment lines
    bool keep_going = true;
    
    while (keep_going) {
        std::getline(myFile,line);
        std::size_t found = line.find("#",0);
        if (found != std::string::npos) {
            //line is a comment line; it can be skipped
            num_lines +=1;
        }
        else{keep_going = false;}
    }
    
    myFile.close();
    myFile.open(file_name.c_str());
    std::string junk;
    
    /* Skip over comment lines */
    for (int i=0; i<num_lines; i++) {
        myFile >> junk;
    }
    
    /* Read in data*/
    for (int i=0; i<num; i++) {
        myFile >> X(i);
    }
 
    myFile.close();
}
void Dynamics::set_dt(double dtIN){
    dt = dtIN;
    dt_is_set = true;
}
void Dynamics::set_total_time(double total_timeIN){
    
    if (!dt_is_set) {
        std::cout << "ERROR: dt must be set before setting total_time." << std::endl;
    }
    
    else{
        total_time = total_timeIN;
        num_steps = total_time/dt;
        total_time_is_set = true;
    }
}
void Dynamics::format_array(vector<vector<double> > &X, vector<double> &X_local){

    int s = 0; //stride

    for(int traj=0; traj<num_trajs_local; traj++){
        X(traj).resize(nuc_beads);
        s = traj*nuc_beads;

        for(int bead=0; bead<nuc_beads; bead++){
            X(traj)(bead) = X_local(s + bead);
        }
    }
}
void Dynamics::format_array(vector<matrix<double> > &X, vector<double> &X_local){

    int s1 = 0; //stride
    int s2 = 0; //stride

    for(int traj=0; traj<num_trajs_local; traj++){
        X(traj).resize(elec_beads,num_states);
        s1 = traj*elec_beads*num_states;

        for(int bead=0; bead<elec_beads; bead++){
            s2 = bead*num_states;
            
            for(int state=0; state<num_states; state++){
                X(traj)(bead,state) = X_local(s1 + s2 + state);
            }
        }
    }
}
void Dynamics::set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                          double massIN, double beta_nuc_beadsIN,
                          double beta_elec_beadsIN, double alphaIN){
    
    nuc_beads = nuc_beadsIN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN:
    beta_nuc_beads = beta_nuc_beadsIN;
    beta_elec_beads = beta_nuc_beadsIN;
    alpha = alphaIN;
}
