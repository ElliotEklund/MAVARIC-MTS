#include "csrpmd_sampling.hpp"

csrpmd_sampling::csrpmd_sampling(int my_id, int num_procs, int root_proc)
    :my_id(my_id),num_procs(num_procs),root_proc(root_proc),
     myRand(time(NULL) + my_id)
{
    sys_set = false;
    sample_set = false;
}
void csrpmd_sampling::run(std::string file_name){
    
    if (my_id==root_proc) {
        if (!sys_set) {
            std::cout << "ERROR: system variables have not been set!" <<
            std::endl;}
        if (!sample_set) {
            std::cout << "ERROR: sampling variables have not been set!" <<
            std::endl;}
        if ((num_trajs % num_procs) != 0) {
            std::cout << "ERROR: num_trajs must be divisible by num_trajs!" <<
            std::endl;
        }
    }

    csrpmd_forces F(nuc_beads,elec_beads,num_states,mass,beta/nuc_beads);
    mv_forces_temp *Fp = &F;
    RK4_MVRPMD rk4(Fp,nuc_beads,elec_beads,num_states,dt);
    
    int num_trajs_pn = num_trajs/num_procs; //num_trajs per node
    int rate_steps = int(rate/dt); //number of steps of size dt between kicks
    int num_kicks = int(decorr/rate); //number of kicks for given decorrelation
    
    double stdev = sqrt(mass/(beta*nuc_beads));
    Normaldev_BM momentum(0, stdev, myRand.int32());
    
    vector<double> Q_traj(nuc_beads,0.0),P_traj(nuc_beads,0.0);
    matrix<double> x_traj(elec_beads,num_states,0.0);
    matrix<double> p_traj(elec_beads,num_states,0.0);
    
    vector<double> Q(num_trajs_pn*nuc_beads,0);
    vector<double> P(num_trajs_pn*nuc_beads,0);
    matrix<double> x(num_trajs_pn*elec_beads,num_states,0);
    matrix<double> p(num_trajs_pn*elec_beads,num_states,0);
    
    init_Q(Q_traj,2.0);
    init_P(P_traj);
    init_elec(x_traj,p_traj);
    //TODO: Initialize x an p some how
    
    for (int traj=0; traj<num_trajs_pn; traj++) {
        for (int kick=0; kick<num_kicks; kick++) {
            //kick momentum
            for (int bead=0; bead<nuc_beads; bead++) {
                P_traj(bead) = momentum.dev();
            }
            for (int step=0; step<rate_steps; step++) {
                rk4.take_step(Q_traj,P_traj,x_traj,p_traj);
            }
        }
        //save trajectory
        for (int bead=0; bead<nuc_beads; bead++) {
            Q(traj*nuc_beads+bead) = Q_traj(bead);
             P(traj*nuc_beads+bead) = P_traj(bead);
         }
        for (int bead=0; bead<elec_beads; bead++) {
            for (int state=0; state<num_states; state++) {
                x(traj*elec_beads+bead,state) = x_traj(bead,state);
                p(traj*elec_beads+bead,state) = p_traj(bead,state);
            }
        }
    }
    
    std::cout << Q << std::endl;
    std::cout << P << std::endl;
    std::cout << x << std::endl;
    std::cout << p << std::endl;

    save_trajs(Q,file_name + "Q");
    save_trajs(P,file_name + "P");
    save_trajs(x,elec_beads*num_states,num_trajs_pn,file_name + "xElec");
    save_trajs(p,elec_beads*num_states,num_trajs_pn,file_name + "pElec");

}

void csrpmd_sampling::init_Q(vector<double> &Q, double ss){
    for (int bead=0; bead<nuc_beads; bead++) {
        Q(bead) = step_dist(myRand.doub(),ss);
    }
}
void csrpmd_sampling::init_P(vector<double> &P){
    
    //system momentum distribution from Gaussian(mu, sigma, seed)
    double stdev = sqrt(mass/(beta*nuc_beads));
    Normaldev_BM momentum(0, stdev, myRand.int32());
    
    for (int bead=0; bead<nuc_beads; bead++) {
        P(bead) = momentum.dev();
    }
}
void csrpmd_sampling::init_elec(matrix<double> &x, matrix<double> &p){
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            x(bead) = step_dist(myRand.doub(),2.0);
            p(bead) = step_dist(myRand.doub(),2.0);
        }
    }
}
inline double csrpmd_sampling::step_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
void csrpmd_sampling::set_sys_vars(int nuc_beadsIN, int elec_beadsIN,
                                   int num_statesIN,double betaIN,double massIN){
    
    nuc_beads = nuc_beadsIN;
    elec_beads = elec_beadsIN;
    num_states = num_statesIN;
    beta = betaIN;
    mass = massIN;
    sys_set = true;
}
void csrpmd_sampling::set_sample_vars(int num_trajsIN, double decorrIN, double rateIN,
                                      double dtIN){
    num_trajs = num_trajsIN;
    decorr = decorrIN;
    rate = rateIN;
    dt = dtIN;
    sample_set = true;
}
void csrpmd_sampling::save_trajs(vector<double> &v,std::string name){

    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    std::cout << "SAVE ME!!" << std::endl;
    myWrap.write_vector(v,name);
}
void csrpmd_sampling::save_trajs(matrix<double> &v,int size,int trajs,
                                 std::string name){

    vector<double> v_transform (size,0);
    int stride = 0;

    for (int traj=0; traj<trajs; traj++) {
        for (int bead=0; bead<elec_beads; bead++) {
            stride = traj*elec_beads*num_states + bead*num_states;
            for (int state=0; state<num_states; state++) {
                v_transform(stride + state) = v(traj*elec_beads+bead,state);
            }
        }
    }

    std::string fileName = name;
    mpi_wrapper myWrap(num_procs,my_id,root_proc);
    myWrap.write_vector(v_transform,fileName);
}
