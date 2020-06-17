#include "init_PAC.hpp"

init_PAC::init_PAC(int my_id,int num_procs,int root_proc,int num_trajs_total,
                   int num_trajs_local)
    :num_procs(num_procs),
     my_id(my_id),
     root_proc(root_proc),
   
     num_trajs_total(num_trajs_total),
     num_trajs_local(num_trajs_local),

     theta_vec(num_trajs_local,0),
     qqTheta_vec(num_trajs_local,0)
{}
void init_PAC::compute_vecs(std::string input_dir){
    
    C_Matrix C(elec_beads, num_states,alpha);
    M_Matrix M(num_states, elec_beads, beta/elec_beads);
    theta_mixed Theta(num_states,nuc_beads,elec_beads,C,M);
    
    vector<vector<double> > Q(num_trajs_local,zero_vector<double>(nuc_beads));
    vector<matrix<double> > x(num_trajs_local,zero_matrix<double>(elec_beads,num_states));
    vector<matrix<double> > p(num_trajs_local,zero_matrix<double>(elec_beads,num_states));
    
    Q = get_trajs_reformat(input_dir + "Q",num_trajs_total*nuc_beads,
                           num_trajs_local*nuc_beads,my_id,num_procs,root_proc,
                           num_trajs_local,nuc_beads);
    
    x = get_trajs_reformat(input_dir + "xElec",num_trajs_total*elec_beads*num_states,
                           num_trajs_local*elec_beads*num_states,my_id,num_procs,
                           root_proc,num_trajs_local,elec_beads,num_states);
    
    p = get_trajs_reformat(input_dir + "pElec",num_trajs_total*elec_beads*num_states,
                           num_trajs_local*elec_beads*num_states,my_id,num_procs,
                           root_proc,num_trajs_local,elec_beads,num_states);
    
    vector<double> Q_traj (nuc_beads);
    matrix<double> x_traj (elec_beads,num_states);
    matrix<double> p_traj (elec_beads,num_states);
    
    double sgntheta = 0;
    double centroid = 0;
    
    for (int traj=0; traj<num_trajs_local; traj++) {
        
        Q_traj = Q(traj);
        x_traj = x(traj);
        p_traj = p(traj);
        
        sgntheta = Theta.get_signTheta(Q_traj,x_traj,p_traj);
        theta_vec[traj] = sgntheta;
        centroid = get_centroid(Q_traj);
        qqTheta_vec[traj] = centroid*centroid*sgntheta;
    }
}
void init_PAC::compute(std::string input_dir,std::string output_dir){
    
    compute_vecs(input_dir);
    
    int num_samples = num_trajs_total/interval; //number of samples that will be generated
    int batch_size = interval/num_procs; //number of samples per processor
    double sgntheta;
    double qqTheta;
    int stride;
    
    vector<double> theta_samples (num_samples,0.0);
    vector<double> qqTheta_samples (num_samples,0.0);
    
    /* Slide through all possible samples at given interval */
    for (int sample=0; sample<num_samples; sample++) {
        sgntheta = 0;
        qqTheta = 0;

        /* sgntheta and qqTheta use samples separated by length interval*/
        for (int traj=0; traj<batch_size; traj++) {
            sgntheta += theta_vec[sample  + traj*num_samples];
            qqTheta += qqTheta_vec[sample + traj*num_samples];

        }
        theta_samples[sample] = sgntheta;
        qqTheta_samples[sample] = qqTheta;
    }
        
    vector<double> theta_samples_final (num_samples,0);
    vector<double> qqTheta_samples_final (num_samples,0);
    
    MPI_Reduce(&theta_samples[0],&theta_samples_final[0],num_samples,
               MPI_DOUBLE,MPI_SUM,root_proc,MPI_COMM_WORLD);
    
    MPI_Reduce(&qqTheta_samples[0],&qqTheta_samples_final[0],num_samples,
               MPI_DOUBLE,MPI_SUM,root_proc,MPI_COMM_WORLD);
    
    if (my_id == root_proc) {

        std::string fileName = output_dir + "iPAC";
        std::ofstream myFile;
        myFile.open(fileName.c_str());
        
        if (!myFile.is_open()) {
            std::cout << "ERROR: Could not save " << fileName << std::endl;
        }
        
        myFile << "#interval:" << interval << std::endl;
        myFile << "#samples per data point:" << num_samples << std::endl;
        
        double theta_sum = 0;
        double qq_sum = 0;
        
        for (int i=0; i<num_samples; i++) {
            theta_sum += theta_samples_final[i];
            qq_sum += qqTheta_samples_final[i];
            double error = 0;
            
            for (int j=0; j<i; j++) {
                error = pow((qq_sum/theta_sum) - (qqTheta_samples_final[i]/theta_samples_final[i]),2);
            }
            
            myFile << (i+1)*interval << " " << qq_sum/theta_sum << " " << theta_sum
            << " " << sqrt(error/(i+1)) << std::endl;
        }
        myFile.close();
    }
}
double init_PAC::get_centroid(const vector<double> & Q_IN){
    double centroid = inner_prod(Q_IN,ones)/nuc_beads;
    return centroid;
}

void init_PAC::set_interval(int interval_IN){interval = interval_IN;}

void init_PAC::set_vectors(vector<vector<double> > Q_IN, vector<matrix<double> > x_IN,
                           vector<matrix<double> > p_IN){
    Q = Q_IN;
    x = x_IN;
    p = p_IN;
}
void init_PAC::set_system(int nuc_beadsIN,int elec_beadsIN,int num_statesIN,
                          double betaIN,double alphaIN){
    nuc_beads = nuc_beadsIN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    beta = betaIN;
    alpha = alphaIN;
    ones.resize(nuc_beads);
    
    for (int i=0; i<nuc_beads; i++) {
        ones(i) = 1.0;
    }
}
