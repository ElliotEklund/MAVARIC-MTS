#include "two_particle_Dynamics.hpp"

two_particle_Dynamics::two_particle_Dynamics(int my_id, int root_proc, int num_procs)
    :my_id(my_id), root_proc(root_proc), num_procs(num_procs)
{
}
void two_particle_Dynamics::initialize_system(int num_beads1IN, int num_beads2IN,
                                              double massIN, double betaIN){
    num_beads1 = num_beads1IN;
    num_beads2 = num_beads2IN;
    mass = massIN;
    beta = betaIN;
}
void two_particle_Dynamics::initialize_dynamics(double dtIN, double total_tIN,
                                                int num_trajsIN,std::string rootFolderIN){
    dt = dtIN;
    total_t = total_tIN;
    num_trajs = num_trajsIN;
    rootFolder = rootFolderIN;
}
vector<vector <double> > two_particle_Dynamics::get_trajs(std::string file,std::string rootPath,
                                                int num_trajs,int num_beads){

    int num_trajs_local = num_trajs/num_procs;
    vector<double> x_global = zero_vector<double> (num_trajs*num_beads);
    vector<double> x_local  = zero_vector<double> (num_trajs_local*num_beads);
    vector<vector<double> > x (num_trajs_local);
    std::string root = rootPath + "OutPut/";
    
    for (int traj=0; traj<num_trajs_local; traj++) {
        x(traj).resize(num_beads,0);
    }
    
    if(my_id == root_proc){
        load_var(x_global,file,root);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* Distribute global trajectories across all processors. */
    MPI_Scatter(&x_global[0],num_trajs_local*num_beads,MPI_DOUBLE,&x_local[0],
                num_trajs_local*num_beads,MPI_DOUBLE,root_proc,MPI_COMM_WORLD);
    
    format_array(num_beads,num_trajs_local,x,x_local);
    
    return x;
}
void two_particle_Dynamics::load_var(vector<double> &X, std::string var, std::string root_path){
      
    /* Add root to specific file indicator*/
    std::string file_name =  root_path + var;
    std::ifstream myFile;
    myFile.open(file_name.c_str());
    
    if(!myFile.is_open()) {
        std::cout << "ERROR: Could not open file " << var << std::endl;
    }
 
    int num = X.size(); //size of vector
    
    for (int i=0; i<num; i++) {
        myFile >> X(i);
    }
    myFile.close();
}
void two_particle_Dynamics::format_array(int num_beads, int num_trajs_local,
                                         vector<vector<double> > &X,
                                         vector<double> &X_local){

    int s = 0; //stride

    for(int traj=0; traj<num_trajs_local; traj++){
        X(traj).resize(num_beads);
        s = traj*num_beads;

        for(int bead=0; bead<num_beads; bead++){
            X(traj)(bead) = X_local(s + bead);
        }
    }
}
void two_particle_Dynamics::run(){

    VV_two_particle VV(num_beads1,num_beads2,mass,beta,dt);
    int num_steps = int(total_t/dt);
    int num_trajs_local = num_trajs/num_procs;
    
    vector<double> cqq1(num_steps,0);
    vector<double> cqq2(num_steps,0);

    vector<vector<double> > Q1(num_trajs_local), Q2(num_trajs_local);
    vector<vector<double> > P1(num_trajs_local), P2(num_trajs_local);
    
    for (int traj=0; traj<num_trajs_local; traj++) {
        Q1(traj).resize(num_beads1,0);
        Q2(traj).resize(num_beads2,0);
        P1(traj).resize(num_beads1,0);
        P2(traj).resize(num_beads2,0);
    }
    
    vector<double> Q1_traj(num_beads1,0), Q2_traj(num_beads2);
    vector<double> P1_traj(num_beads1,0), P2_traj(num_beads2);

    Q1 = get_trajs("Q1",rootFolder,num_trajs,num_beads1);
    Q2 = get_trajs("Q2",rootFolder,num_trajs,num_beads2);
    P1 = get_trajs("P1",rootFolder,num_trajs,num_beads1);
    P2 = get_trajs("P2",rootFolder,num_trajs,num_beads2);
    
    double cent1_0 = 0;
    double cent2_0 = 0;
    double cent1_t = 0;
    double cent2_t = 0;

    for (int traj=0; traj<num_trajs_local; traj++) {
        
        Q1_traj = Q1(traj);
        Q2_traj = Q2(traj);
        P1_traj = P1(traj);
        P2_traj = P2(traj);
        
        cent1_0 = get_centroid(Q1_traj,num_beads1);
        cent2_0 = get_centroid(Q2_traj,num_beads2);
        cent1_t = cent1_0;
        cent2_t = cent2_0;

        for (int step=0; step<num_steps; step++) {
            
            cqq1(step) += cent1_0*cent1_t;
            cqq2(step) += cent2_0*cent2_t;
            
            VV.step(Q1_traj,Q2_traj,P1_traj,P2_traj);
            cent1_t = get_centroid(Q1_traj,num_beads1);
            cent2_t = get_centroid(Q2_traj,num_beads2);
        }
    }
    
    write_correlation("PAC1",cqq1);
    write_correlation("PAC2",cqq2);
}
double two_particle_Dynamics::get_centroid(const vector<double> &x, int num_beads){
    return sum(x)/num_beads;
}
void two_particle_Dynamics::write_correlation(std::string file, vector<double> &QQ){
    
    int num_steps = int(total_t/dt);
    vector<double> QQ_global(num_steps,0);
    MPI_Reduce(&QQ[0],&QQ_global[0],num_steps,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
    if (my_id==root_proc) {
        std::string fileName = rootFolder + "Output/" + file;
        
        std::ofstream myFile;
        myFile.open(fileName.c_str());
        
        if (!myFile.is_open()) {
            std::cout << "ERROR: Could not open " << fileName << std::endl;
        }
        
        for(int i=0; i<QQ.size(); i++){
            myFile << i*dt << " " << QQ_global(i)/num_trajs << std::endl;
        }
        
        myFile.close();
        std::cout << "Successfully wrote pos_auto_corr to Results." << std::endl;
    }
}
