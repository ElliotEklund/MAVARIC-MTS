#include "two_particle_sampling.hpp"

two_particle_sampling::two_particle_sampling(int my_id, int root_proc, int num_procs)
    :/* Initialize parameters*/
     my_id(my_id), root_proc(root_proc), num_procs(num_procs),

     /* Initialize random number generators */
     myRand(time(NULL) + my_id),
     myHelper(my_id,num_procs,root_proc)
{}

void two_particle_sampling::gen_initQ(vector<double> &Q, int num_beads, double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = nuc_dist(myRand.doub(),step_size);
    }
}
void two_particle_sampling::initialize_system(int num_beads1IN, int num_beads2IN,
                                              double massIN, double betaIN){
    num_beads1 = num_beads1IN;
    num_beads2 = num_beads2IN;
    mass = massIN;
    beta = betaIN;
}
void two_particle_sampling::initialize_files(bool readPSV1IN, bool readPSV2IN,
                                             std::string rootFolderIN){
    
    readPSV1 = readPSV1IN;
    readPSV2 = readPSV2IN;
    rootFolder = rootFolderIN;
    myHelper.set_root(rootFolderIN);
}
void two_particle_sampling::runSimulation(double ss1, double ss2, unsigned long long num_trajs,
                                          unsigned long long decorrelation){
    
    SpringEnergy V_spring1(num_beads1,mass,beta/num_beads1);
    SpringEnergy V_spring2(num_beads2,mass,beta/num_beads2);
    two_particle_Hamiltonian H(num_beads1,num_beads2,V_spring1,V_spring2);

    vector<double> Q1(num_beads1,0);
    vector<double> Q_prop1(num_beads1,0);
    vector<double> Q2(num_beads2,0);
    vector<double> Q_prop2(num_beads2,0);
    vector<double> Q1_trajs(num_trajs*num_beads1,0);
    vector<double> Q2_trajs(num_trajs*num_beads2,0);
    vector<double> P1_trajs(num_trajs*num_beads1,0);
    vector<double> P2_trajs(num_trajs*num_beads2,0);

    gen_initQ(Q1,num_beads1,ss1);
    gen_initQ(Q2,num_beads2,ss2);
    
    if (readPSV1) {myHelper.read_PSV(num_beads1,Q1,"1");}
    if (readPSV2) {myHelper.read_PSV(num_beads2,Q2,"2");}

    Q_prop1 = Q1;
    Q_prop2 = Q2;
 
    double energy = H.get_energy(Q1,Q2);
    double energy_prop = energy;
    
    int mcMove = 0;
    bool accept_move = false;

    double r = double(num_beads2)/double(num_beads1);
    if (r==1) {
        r = 0.5;
    }

    unsigned long long sys1_steps_accpt = 0;
    unsigned long long sys1_steps_total = 0;
    unsigned long long sys2_steps_accpt = 0;
    unsigned long long sys2_steps_total = 0;
    
    for (int traj=0; traj<num_trajs; traj++) {
        for (int step=0; step<decorrelation; step++) {
            if(myRand.doub() > r){
                
                //Sample system1 coordinates
                for (int i=0; i<num_beads1; i++) {
                    mcMove = rand_bead(myRand.int64(),num_beads1);
                    Q_prop1(mcMove) = Q1(mcMove) + nuc_dist(myRand.doub(),ss1);
                    energy_prop = H.get_energy(Q_prop1,Q2);
                    accept_move = check_move(energy,energy_prop);
                    
                    if (accept_move) {
                        Q1(mcMove) = Q_prop1(mcMove);
                        energy = energy_prop;
                        sys1_steps_accpt += 1;
                    }
                    else{
                        Q_prop1(mcMove) = Q1(mcMove);
                    }
                    sys1_steps_total += 1;
                }
            }
            else{
                
                //Sample system2 coordinates
                for (int i=0; i<num_beads1; i++) {
                    mcMove = rand_bead(myRand.int64(),num_beads2);
                    Q_prop2(mcMove) = Q2(mcMove) + nuc_dist(myRand.doub(),ss2);
                    energy_prop = H.get_energy(Q1,Q_prop2);
                    accept_move = check_move(energy,energy_prop);
                    
                    if (accept_move) {
                        Q2(mcMove) = Q_prop2(mcMove);
                        energy = energy_prop;
                        sys2_steps_accpt += 1;
                    }
                    else{
                        Q_prop2(mcMove) = Q2(mcMove);
                    }
                    sys2_steps_total += 1;
                }
            }
        }
        
        for (int bead=0; bead<num_beads1; bead++) {
            Q1_trajs(traj*num_beads1+bead) = Q1(bead);
        }
        for (int bead=0; bead<num_beads2; bead++) {
            Q2_trajs(traj*num_beads2+bead) = Q2(bead);
        }
    }
    
    double stdev1 = sqrt(mass/(beta*num_beads1));
    double stdev2 = sqrt(mass/(beta*num_beads2));

    Normaldev_BM mom1(0, stdev1, rand()); //system momentum distribution from Gaussian(mu, sigma, seed)
    Normaldev_BM mom2(0, stdev2, rand()); //system momentum distribution from Gaussian(mu, sigma, seed)

    for(int i=0; i<num_trajs*num_beads1; i++){
        P1_trajs(i) = mom1.dev();//mom1.dev();
    }
    for(int i=0; i<num_trajs*num_beads2; i++){
        P2_trajs(i) = mom2.dev();
    }
    
    myHelper.print_sys_accpt(sys1_steps_total,sys1_steps_accpt,"1");
    myHelper.print_sys_accpt(sys2_steps_total,sys2_steps_accpt,"2");
    
    save_trajs(Q1_trajs,num_beads1*num_trajs,"Q1");
    save_trajs(Q2_trajs,num_beads2*num_trajs,"Q2");
    save_trajs(P1_trajs,num_beads1*num_trajs,"P1");
    save_trajs(P2_trajs,num_beads2*num_trajs,"P2");
}
bool two_particle_sampling::check_move(double energy, double energy_prop){

    double delta_energy = energy_prop - energy;
    double accept_move = false;
    /* Accept new system moves if energ_prop < energy*/
    if(delta_energy < 0){
        accept_move = true;
    }
    /* Accept new system moves if inequality is met*/
    else if (myRand.doub() <= exp(-beta * delta_energy)){
        accept_move = true;
    }
    return accept_move;
}
inline double two_particle_sampling::nuc_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
inline int two_particle_sampling::rand_bead(const Ullong rn, int num_beads){
    return rn % num_beads;
}
void two_particle_sampling::save_trajs(const vector<double> &v,int size, std::string name){

    std::string fileName = rootFolder + name;
    std::ofstream myFile;
    myFile.open(fileName);
    
    if(!myFile.is_open()){
        std::cout << "ERROR: Could not open " << fileName << std::endl;
    }

    for (int i=0; i<size; i++) {
        myFile << v(i) << std::endl;
    }
    
    myFile.close();
}
