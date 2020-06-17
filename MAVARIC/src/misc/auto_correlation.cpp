#include "auto_correlation.hpp"

auto_correlation::auto_correlation(int my_id,int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc)
{
    is_sys_set = false;
    is_req_set = false;
    is_time_set = false;
}
int auto_correlation::compute(unsigned long long num_trajs_global,
                               unsigned long long num_trajs_local,
                               std::string input_dir, std::string output_dir,
                               int num_samples, int num_errors){
    
    if (num_samples> num_procs) {
        if(my_id == root_proc){
            std::cout << "ERROR: num_samples must be >= num_procs. Aborting "
            "calculations." << std::endl;
        }
        return -1;
    }
    
    /* Initialize Forces and Integrator*/
    C_Matrix C(elec_beads, num_states,alpha);
    M_Matrix M(num_states, elec_beads, beta_elec_beads);
    dM_Matrix_dQ dMdQ(elec_beads, num_states, beta_elec_beads, M);
    theta_mixed theta(num_states,nuc_beads,elec_beads,C,M);
    theta_mixed_dQ theta_dQ(num_states,nuc_beads,elec_beads,C,M,dMdQ);
    theta_mixed_dElec theta_dElec(num_states,elec_beads,alpha,C,M);

    mvrpmd_mixed_forces F(nuc_beads, elec_beads, num_states, mass,
                          beta_nuc_beads, alpha, theta, theta_dQ, theta_dElec);
    
    ABM_MVRPMD myABM(F,dt,num_states,nuc_beads,elec_beads);

    /* Initialize trajectories*/
    vector<double> Q_traj (nuc_beads), P_traj (nuc_beads);
    matrix<double> x_traj (elec_beads,num_states), p_traj (elec_beads,num_states);
    
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
    
    /* Initialize auto-correlation function machinery*/
    aggregate myAggregator(my_id,num_procs,root_proc);
    pop_estimators myPops(elec_beads,num_states,alpha);

    vector<double> pac_v0, bp_v0, sp_v0, wp_v0;
    vector<double> pac_v, bp_v, sp_v, wp_v;
    
    int num_steps = floor(total_time/dt);
    int ten_p = floor(num_trajs_local/10.0);
    std::ofstream progress;
    if (my_id==root_proc) {
        std::string fileName = output_dir + "dyn_progress";
        progress.open(fileName.c_str());
        if (!progress.is_open()) {
            std::cout << "ERROR: Could now open " << fileName << std::endl;
        }
    }

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
            /* Use semi classical for time zero */
            bp_v0 = myPops.sc(x_traj,p_traj);
            bp_v = myPops.boltz(theta.get_gamm());
            myAggregator.collect("boltzman",0,bp_v0,bp_v,sgnTheta);
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
            myAggregator.merge_collections(root_proc,my_id,output_dir,dt,sp_stride,
                                           traj*num_procs);
            if (pac) {
                myAggregator.write_errors("position",num_samples,num_errors,dt,pac_stride,output_dir);
            }
            if (bp) {
                myAggregator.write_errors("boltzman",num_samples,num_errors,dt,pac_stride,output_dir);
            }
            if (sp) {
                myAggregator.write_errors("semi_classic",num_samples,num_errors,dt,pac_stride,output_dir);
            }
            if (wp) {
                myAggregator.write_errors("wigner",num_samples,num_errors,dt,pac_stride,output_dir);
            }
        }
    }

    if (my_id==root_proc) {
        progress.close();
    }
    
    /* NOTE sp_strides currently sets the stride for all other ac functions*/
    myAggregator.merge_collections(root_proc,my_id,output_dir,dt,sp_stride,
                                   num_trajs_global);
    if (pac) {
        myAggregator.write_errors("position",num_samples,num_errors,dt,pac_stride,output_dir);}
    if (bp) {
        myAggregator.write_errors("boltzman",num_samples,num_errors,dt,pac_stride,output_dir);}
    if (sp) {
        myAggregator.write_errors("semi_classic",num_samples,num_errors,dt,pac_stride,output_dir);}
    if (wp) {
        myAggregator.write_errors("wigner",num_samples,num_errors,dt,pac_stride,output_dir);}
    
    return 0;
}
double auto_correlation::compute_centroid(const vector<double> &Q){
    double centroid = sum(Q);
    return centroid/nuc_beads;
}
void auto_correlation::set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
                                 double massIN, double beta_nuc_beadsIN,
                                 double beta_elec_beadsIN, double alphaIN){
    nuc_beads = nuc_beadsIN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    mass = massIN;
    beta_nuc_beads = beta_nuc_beadsIN;
    beta_elec_beads = beta_nuc_beadsIN;
    alpha = alphaIN;
    is_sys_set = true;
}
void auto_correlation::request_calcs(bool pacIN, int pac_strideIN, bool bpIN,
                                     int bp_strideIN,bool spIN, int sp_strideIN,
                                     bool wpIN,int wp_strideIN){

    pac = pacIN;
    bp = bpIN;
    wp= wpIN;
    sp = spIN;
    pac_stride = pac_strideIN;
    bp_stride = bp_strideIN;
    wp_stride = wp_strideIN;
    sp_stride = sp_strideIN;
    is_req_set = true;
}
void auto_correlation::set_time(double dtIN, double total_timeIN){
    dt = dtIN;
    total_time = total_timeIN;
    is_time_set = true;
}
