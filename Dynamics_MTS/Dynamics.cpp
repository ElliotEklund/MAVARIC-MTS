#include "Dynamics.hpp"

Dynamics::Dynamics(int num_procs, int my_id, int root_proc,int nuc_beads, 
                   int elec_beads, int num_states, double mass,
                   double beta_nuc_beads, double beta_elec_beads, int num_trajs)
    :num_procs(num_procs),my_id(my_id),root_proc(root_proc),
        
     nuc_beads(nuc_beads),elec_beads(elec_beads),num_states(num_states),
     num_trajs(num_trajs),num_trajs_local(num_trajs/num_procs),
    
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
    vector<double> x_local = zero_vector<double> (num_trajs_local*elec_beads*nuc_beads);
    vector<double> p_local = zero_vector<double> (num_trajs_local*elec_beads*nuc_beads);
    
    std::string root_path = "/Users/ellioteklund/Desktop/Dynamics_MTS_AlphaDEBUG/Dynamics_MTS/Results/Trajectories/";
    //std::string root_path = "/home/elliot/Desktop/Dynamics_MTS_AlphaDEBUG/Dynamics_MTS/Results/Trajectories/";
    
    /* Read in trajectories to global PSV vectors  */
    if(my_id == root_proc){
        load_var(Q_global,"Q0",root_path);
        load_var(P_global,"P0",root_path);
        load_var(x_global,"xelec0",root_path);
        load_var(p_global,"pelec0",root_path);
    }
    
   MPI_Barrier(MPI_COMM_WORLD);

   /* Distribute global trajectories across all processors. */
   MPI_Scatter(&Q_global[0],num_trajs_local*nuc_beads,MPI_DOUBLE,&Q_local[0],num_trajs*nuc_beads,
     MPI_DOUBLE,root_proc,MPI_COMM_WORLD);

   MPI_Scatter(&P_global[0],num_trajs_local*nuc_beads,MPI_DOUBLE,&P_local[0],num_trajs*nuc_beads,
     MPI_DOUBLE,root_proc,MPI_COMM_WORLD);

   MPI_Scatter(&x_global[0],num_trajs_local*elec_beads*num_states,MPI_DOUBLE,&x_local[0],
     num_trajs*elec_beads*num_states,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);

   MPI_Scatter(&p_global[0],num_trajs_local*elec_beads*num_states,MPI_DOUBLE,&p_local[0],
     num_trajs*elec_beads*num_states,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);

    /* Trajectories are now distributed across all processors, but still need to be
     * correctly formated. */
   format_array(Q,Q_local);
   format_array(P,P_local);
   format_array(x,x_local);
   format_array(p,p_local);
}

void Dynamics::runSimulation(){
    
    int num_steps = int(total_time/dt);
//    int data_count = int(total_time * 10); //resolution of data being sampled
//    int rate = num_steps/data_count;

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
    
    //QQt = zero_vector<double>(data_count);
    int i_data = 0;
    QQt = zero_vector<double>(num_steps);

    std::ofstream traj_file;
    traj_file.open("./traj_file_MTS");

    for (int traj=0; traj<num_trajs_local; traj++){
        
        //std::cout << std::endl << "------------------ Traj " << traj << "----------------" << std::endl << std::endl;
       
        /* Load new trajecty*/
        Q_traj = Q(traj);
        P_traj = P(traj);
        x_traj = x(traj);
        p_traj = p(traj);


//        std::cout << "Q" << std::endl;
//        std::cout << Q_traj << std::endl;
//
//        std::cout << "P" << std::endl;
//        std::cout << P_traj << std::endl;
//
//        std::cout << "x" << std::endl;
//        std::cout << x_traj << std::endl;
//
//        std::cout << "p" << std::endl;
//        std::cout << p_traj << std::endl;


        Qcent_0 = compute_centroid(Q_traj);

        myABM.initialize_rk4(Q_traj, P_traj, x_traj, p_traj);
    
        sgnTheta = F.get_sgnTheta(Q_traj,x_traj,p_traj);
        
        QQt(0) += Qcent_0 * Qcent_0 * sgnTheta;
 
        sgnTheta_total += sgnTheta;
        Qcent0_X_sgnTheta = Qcent_0 * sgnTheta;

        i_data = 1;


        traj_file << Qcent_0 * Qcent_0 *sgnTheta << " ";


        for (int step=1; step<num_steps; step++) {
           
            //std::cout << "step: " << step << std::endl;

            myABM.take_step(Q_traj, P_traj, x_traj, p_traj);
            Qcent_t = compute_centroid(Q_traj);
            QQt(step) += Qcent_t * Qcent0_X_sgnTheta;
            traj_file << Qcent_t * Qcent0_X_sgnTheta << " " ;


//            std::cout << "Q" << std::endl;
//            std::cout << Q_traj << std::endl;
//
//            std::cout << "P" << std::endl;
//            std::cout << P_traj << std::endl;
//
//            std::cout << "x" << std::endl;
//            std::cout << x_traj << std::endl;
//
//            std::cout << "p" << std::endl;
//            std::cout << p_traj << std::endl;


           // if(step%rate == 0){
           //     Qcent_t = compute_centroid(Q_traj);
           //     QQt(i_data) += Qcent_t * Qcent0_X_sgnTheta;
           //     i_data += 1;
           // }
        }
        traj_file << std::endl;
    }

    print_QQt(QQt,sgnTheta_total);
    traj_file.close();

}

void Dynamics::print_QQt(const vector<double> &QQt,double sgnTheta_total){


    int data_count = total_time/dt;

    //int data_count = total_time*10;
    //int rate = 1/(dt*10);

    double global_sgn_theta_total = 0;
    vector<double> global_QQt (data_count);
    
    MPI_Reduce(&sgnTheta_total,&global_sgn_theta_total,1,
    MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    MPI_Reduce(&QQt[0],&global_QQt[0],data_count,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if(my_id == root_proc){
        std::ofstream myFile;
        myFile.open("./pos_auto_corr");
    
        for(int i=0; i<QQt.size(); i++){
            //myFile << i*rate*dt << " " << global_QQt(i)/global_sgn_theta_total << std::endl;
            myFile << i*dt << " " << global_QQt(i)/global_sgn_theta_total << std::endl;
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
    myFile.open(file_name);
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file " << var << std::endl;
    }
 
    int num = X.size(); //size of vector
    
    for (int i=0; i<num; i++) {
        myFile >> X(i);
    }
    myFile.close();
}

void Dynamics::set_dt(double dtIN){dt = dtIN;}

void Dynamics::set_total_time(double total_timeIN){total_time = total_timeIN;}

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
