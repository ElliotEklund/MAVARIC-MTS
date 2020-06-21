#include "dynamics_mvrpmd.hpp"

dynamics_mvrpmd::dynamics_mvrpmd(int my_id, int num_procs, int root_proc)
    :my_id(my_id),
     num_procs(num_procs),
     root_proc(root_proc)
{
    is_sys_set = false;
    is_time_set = false;
    is_trajs_set = false;
}
int dynamics_mvrpmd::pre_comp(){
    int abort = 0;
    
    if (!is_sys_set) {
        abort = -1;
        if (my_id==root_proc) {
            std::cout << "ERROR: System variables are not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
    }
    if (!is_time_set) {
        abort = -1;
        if (my_id==root_proc) {
            std::cout << "ERROR: Time variables are not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
    }
    if (!is_trajs_set) {
        abort = -1;
        if (my_id==root_proc) {
            std::cout << "ERROR: Trajectory variables are not set!" << std::endl;
            std::cout << "Aborting calculation." << std::endl;
        }
    }
    return abort;
}
int dynamics_mvrpmd::compute_ac(bool pac, int pac_stride, bool bp,
                                 int bp_stride,bool sp, int sp_stride,
                                 bool wp,int wp_stride,std::string input_dir,
                                 std::string output_dir,int num_samples,
                                 int num_errors){
    
    if(pre_comp()==-1){
        return -1;
    }
    
    auto_correlation my_ac(my_id,num_procs,root_proc);
    my_ac.set_system(nuc_beads,elec_beads,num_states,mass,
                     beta_nuc_beads,beta_elec_beads,alpha);
    
    my_ac.request_calcs(pac,pac_stride,bp,bp_stride,sp,sp_stride,wp,wp_stride);
    my_ac.set_time(dt,total_time);
    
    my_ac.compute(num_trajs_global,num_trajs_local,input_dir,output_dir,
                  num_samples,num_errors);
    
    return 0;
}

int dynamics_mvrpmd::energy_conserve(double tol, int energy_stride,
                                      std::string input_dir,
                                      std::string output_dir){
    
    if(pre_comp()==-1){
        return -1;
    }
    
    energy_conserv my_conserv(my_id,num_procs,root_proc);
    
    my_conserv.set_system(nuc_beads,elec_beads,num_states,mass,beta,
                          beta_nuc_beads,beta_elec_beads,alpha);
    
    my_conserv.set_time(dt,total_time);
    
    my_conserv.compute(num_trajs_global,num_trajs_local,tol,energy_stride,
                       input_dir,output_dir);
    
    return 0;
}

int dynamics_mvrpmd::iPAC(int interval,std::string input_dir,std::string output_dir){
    
    if (my_id == root_proc) {
        std::cout << "HACK ALERT: iPAC calculation only supports even ratios of "
        "num_trajs, interval, and num_procs for the time being." << std::endl;
    }
    
    if(pre_comp()==-1){
        return -1;
    }
    
    init_PAC my_init_PAC(my_id,num_procs,root_proc,num_trajs_global,num_trajs_local);
    my_init_PAC.set_system(nuc_beads,elec_beads,num_states,beta,alpha);
    my_init_PAC.set_interval(interval);
    my_init_PAC.compute(input_dir,output_dir);
    
    return 0;
}
void dynamics_mvrpmd::set_system(int nuc_beadsIN, int elec_beadsIN, int num_statesIN,
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
void dynamics_mvrpmd::set_time(double dtIN, double total_timeIN){
    dt = dtIN;
    total_time = total_timeIN;
    is_time_set = true;
}
void dynamics_mvrpmd::set_trajs(unsigned long long num_trajs_globalIN,
                                std::string root_path){
    
    num_trajs_global = num_trajs_globalIN;
    num_trajs_local = 0; //number of local trajectories
    
    /* Check for more trajs than procs*/
    if (num_trajs_global < num_procs) {
        if (my_id == root_proc) {
            std::cout << "ERROR: num_trajs is small than num_procs.!" << std::endl;
        }
    }
    else{
        unsigned long long d = num_trajs_global/num_procs;
        unsigned long long r = num_trajs_global%num_procs;
        
        if (my_id < num_procs - r) {
            num_trajs_local = d;
        }
        else{
            num_trajs_local = d + 1;
        }
        is_trajs_set = true;
    }
}
