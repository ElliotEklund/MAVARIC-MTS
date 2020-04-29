#include "Dynamics.hpp"

Dynamics::Dynamics(int num_procs, int my_id, int root_proc,int nuc_beads, 
                   int elec_beads, int num_states, double mass,
                   double beta_nuc_beads, double beta_elec_beads, int num_trajs,std::string root)
    :num_procs(num_procs),my_id(my_id),root_proc(root_proc),
     root(root),
        
     nuc_beads(nuc_beads),elec_beads(elec_beads),num_states(num_states),
     num_trajs(num_trajs),num_trajs_local(num_trajs/num_procs), mass(mass),
     beta_nuc_beads(beta_nuc_beads),
     dt_is_set(false), total_time_is_set(false),
    
     Q(num_trajs_local,zero_vector<double>(nuc_beads)),P(num_trajs_local,zero_vector<double>(nuc_beads)),
     x(num_trajs_local,zero_matrix<double>(elec_beads,num_states)),p(num_trajs_local,zero_matrix<double>(elec_beads,num_states)),

     C(elec_beads, num_states), M(num_states, nuc_beads, beta_nuc_beads),

     M_MTS(nuc_beads, elec_beads, num_states, M),
     dMdQ(nuc_beads, num_states, beta_elec_beads, M),
     dM_MTS_dQ(nuc_beads, elec_beads, num_states, dMdQ),

     Theta(num_states, elec_beads, C, M_MTS),
     dThetadQ(num_states, nuc_beads, elec_beads, C, M_MTS, dM_MTS_dQ),
     dThetadElec(num_states, elec_beads, C, M_MTS),

     F(nuc_beads, elec_beads, num_states, mass,
       beta_nuc_beads, Theta, dThetadQ, dThetadElec)
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
    
    std::string fileName = root + "/Results/Trajectories/";

    /* Read in trajectories to global PSV vectors  */
    if(my_id == root_proc){
        load_var(Q_global,"Q0",fileName);
        load_var(P_global,"P0",fileName);
        load_var(x_global,"xelec0",fileName);
        load_var(p_global,"pelec0",fileName);
    }
    
   MPI_Barrier(MPI_COMM_WORLD);

   /* Distribute global trajectories across all processors. */
   MPI_Scatter(&Q_global[0],num_trajs_local*nuc_beads,MPI_DOUBLE,&Q_local[0],num_trajs_local*nuc_beads,
     MPI_DOUBLE,root_proc,MPI_COMM_WORLD);

   MPI_Scatter(&P_global[0],num_trajs_local*nuc_beads,MPI_DOUBLE,&P_local[0],num_trajs_local*nuc_beads,
     MPI_DOUBLE,root_proc,MPI_COMM_WORLD);

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
                         num_states,num_trajs,num_trajs_local,Theta);
    
    my_init_PAC.set_interval(interval);
    my_init_PAC.set_vectors(Q,x,p);
    my_init_PAC.compute(root);
}

void Dynamics::PAC(){
    
    if (!total_time) {
        std::cout << "ERROR: total_time must be set before a simulation can be run." << std::endl;
    }
    
    int data_count = int(total_time * 10); //resolution of data being sampled
    int rate = num_steps/data_count;

    ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);
    
    vector<double> Q_traj (nuc_beads);
    vector<double> P_traj (nuc_beads);
    matrix<double> x_traj (elec_beads,num_states);
    matrix<double> p_traj (elec_beads,num_states);
    
    double Qcent_0 = 0; //position centroid evaluated at t=0
    double Qcent_t = 0; //position centroid evaluated at t
    double sgnTheta = 0; //sign of Theta for a trajectory
    double sgnTheta_total = 0; //sum of sgnTheta across all trajectories
    double Qcent0_X_sgnTheta = 0; //Qcent_0 * sgnTheta
    
    QQt = zero_vector<double>(data_count);
    int i_data = 0;

    for (int traj=0; traj<num_trajs_local; traj++){

        /* Load new trajecty*/
        Q_traj = Q(traj);
        P_traj = P(traj);
        x_traj = x(traj);
        p_traj = p(traj);

        Qcent_0 = compute_centroid(Q_traj);

        myABM.initialize_rk4(Q_traj, P_traj, x_traj, p_traj);

        sgnTheta = F.get_sgnTheta(Q_traj,x_traj,p_traj);

        QQt(0) += Qcent_0 * Qcent_0 * sgnTheta;

        sgnTheta_total += sgnTheta;
        Qcent0_X_sgnTheta = Qcent_0 * sgnTheta;

        i_data = 1;


        for (int step=1; step<num_steps; step++) {

            myABM.take_step(Q_traj, P_traj, x_traj, p_traj);

            if(step%rate == 0){
                Qcent_t = compute_centroid(Q_traj);
                QQt(i_data) += Qcent_t * Qcent0_X_sgnTheta;
                i_data += 1;
            }
        }
    }

    print_QQt(QQt,sgnTheta_total);

}

void Dynamics::energ_conserv(double tol, int energy_stride){
    
    if (!total_time) {
        std::cout << "ERROR: total_time must be set before a simulation can be run." << std::endl;
    }

    ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);

    SpringEnergy V_spring(nuc_beads,mass,beta_nuc_beads);
    StateIndepPot V0(nuc_beads,mass);
    GTerm G(elec_beads,num_states);
    MVRPMD_MTS_Hamiltonian myHam(beta_nuc_beads,V_spring,V0,G,Theta);

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
        energy_init = myHam.get_energy_dyn(mass,Q_traj,P_traj,x_traj,p_traj);
        
        while (step<num_steps && !broken){

            myABM.take_step(Q_traj, P_traj, x_traj, p_traj);
            step += 1;

            if (step % energy_stride == 0) {
                /* Test if trajectory has broken*/

                energy_t = myHam.get_energy_dyn(mass,Q_traj,P_traj,x_traj,p_traj);
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
    }
    
}

void Dynamics::PopAC(){

    if (!total_time) {
        std::cout << "ERROR: total_time must be set before a simulation can be run." << std::endl;
    }
    
    PopulationEstimator popEsti(elec_beads,num_states);
    
    int data_count = int(total_time * 10); //resolution of data being sampled
    int rate = num_steps/data_count;
    double sgnTheta = 0;
    double sgnTheta_total = 0;
    int i_data = 0;
    
    matrix<double> PPt(data_count,num_states,0);
    vector<double> pop_0(num_states,0); //PPt at t=0
    vector<double> pop_t(num_states,0); //PPt at t=t

    
    ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);
    
    vector<double> Q_traj (nuc_beads);
    vector<double> P_traj (nuc_beads);
    matrix<double> x_traj (elec_beads,num_states);
    matrix<double> p_traj (elec_beads,num_states);
    
    for (int traj=0; traj<num_trajs_local; traj++){
        
        /* Load new trajecty*/
        Q_traj = Q(traj);
        P_traj = P(traj);
        x_traj = x(traj);
        p_traj = p(traj);
        
        myABM.initialize_rk4(Q_traj, P_traj, x_traj, p_traj);
        
        sgnTheta = F.get_sgnTheta(Q_traj,x_traj,p_traj);
        sgnTheta_total += sgnTheta;

        popEsti.update_populations(x_traj,p_traj);
        pop_0 = popEsti.get_pop();
        row(PPt,0) += element_prod(pop_0,pop_0)*sgnTheta;
        pop_0 = pop_0*sgnTheta;
        
        i_data = 1;
        
        for (int step=1; step<num_steps; step++) {
            
            myABM.take_step(Q_traj, P_traj, x_traj, p_traj);
            
            if(step%rate == 0){
                
                popEsti.update_populations(x_traj,p_traj);
                pop_t = popEsti.get_pop();
                row(PPt,i_data) += element_prod(pop_0,pop_t);
                
                i_data += 1;
            }
        }
    }
    
    PPt = PPt/sgnTheta_total;
    popEsti.write_populations(PPt,dt,data_count,rate,root);
}

void Dynamics::write_broken(std::list<int> broken,std::string file_root){
    
    std::ostringstream quick_convert;
    quick_convert << my_id;
    
    std::string file_name = file_root + "/Results/broken" + quick_convert.str();

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

void Dynamics::print_QQt(vector<double> &QQt,double sgnTheta_total){

    int data_count = total_time*10;
    int rate = 1/(dt*10);

    double global_sgn_theta_total = 0;
    vector<double> global_QQt (data_count);
    
    MPI_Reduce(&sgnTheta_total,&global_sgn_theta_total,1,
    MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    MPI_Reduce(&QQt[0],&global_QQt[0],data_count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if(my_id == root_proc){
        std::string fileName = root + "/Results/PAC";
        
        std::ofstream myFile;
        myFile.open(fileName.c_str());
        
        if (!myFile.is_open()) {
            std::cout << "ERROR: Could not open " << fileName << std::endl;
        }
        
    
        for(int i=0; i<QQt.size(); i++){
            myFile << i*rate*dt << " " << global_QQt(i)/global_sgn_theta_total << std::endl;
        }
    
        myFile.close();
        std::cout << "Successfully wrote pos_auto_corr to Results." << std::endl;
    }
}

double Dynamics::compute_centroid(const vector<double> &Q){
    double centroid = sum(Q);
    return centroid/nuc_beads;
}

void Dynamics::load_var(vector<double> &X, std::string var, std::string root_path){
      
    /* Add root to specific file indicator*/
    std::string file_name =  root_path + var;
    
    std::ifstream myFile;
    myFile.open(file_name.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file " << var << std::endl;
    }
 
    int num = X.size(); //size of vector
    
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

/* Functions used of debugging. Will be removed in final version*/

void Dynamics::write_Q(std::ofstream &myStream, double step, vector<double> &Q){

  myStream << step << " ";

  for(int bead=0; bead<nuc_beads; bead++){
      myStream << Q[bead] << " ";
  }

  myStream << std::endl;

}

void Dynamics::write_P(std::ofstream &myStream, double step, vector<double> &P){

  myStream << step << " ";

  for(int bead=0; bead<nuc_beads; bead++){
      myStream << P[bead] << " ";
  }

  myStream << std::endl;

}

void Dynamics::write_x(std::ofstream &myStream, double step, matrix<double> &x){

  myStream << step << " ";


  for(int bead=0; bead<nuc_beads; bead++){
      for(int state=0; state<num_states; state++){
          myStream << x(bead,state) << " ";
      }
  }

  myStream << std::endl;

}

void Dynamics::write_p(std::ofstream &myStream, double step, matrix<double> &p){

  myStream << step << " ";

  for(int bead=0; bead<nuc_beads; bead++){
      for(int state=0; state<num_states; state++){
          myStream << p(bead,state) << " ";
      }
  }

  myStream << std::endl;

}
