#include "two_particle_mc.hpp"

two_particle_mc::two_particle_mc(int my_id, int root_proc, int num_procs)
    :/* Initialize parameters*/
     my_id(my_id), root_proc(root_proc), num_procs(num_procs),

     /* Initialize random number generators */
     myRand(time(NULL) + my_id),
     myHelper(my_id,num_procs,root_proc)
{
    /* Generate random initial PSV. If readPSV is true, these will be overwritten
     with saved values during the call to runSimulation.*/
}

void two_particle_mc::gen_initQ(vector<double> &Q, int num_beads, double step_size){
    for (int bead=0; bead<num_beads; bead++) {
        Q(bead) = nuc_dist(myRand.doub(),step_size);
    }
}
void two_particle_mc::initialize_system(int num_beads1IN, int num_beads2IN,
                                        double massIN, double betaIN){
    num_beads1 = num_beads1IN;
    num_beads2 = num_beads2IN;
    mass = massIN;
    beta = betaIN;
}

void two_particle_mc::initialize_files(bool readPSV1IN, bool readPSV2IN, bool readDataIN,
                                       bool writePSV1IN, bool writePSV2IN, bool writeDataIN,
                                       std::string rootFolderIN){
    
    readPSV1 = readPSV1IN;
    readPSV2 = readPSV2IN;
    readData = readDataIN;
    writePSV1 = writePSV1IN;
    writePSV2 = writePSV2IN;
    writeData = writeDataIN;
    rootFolder = rootFolderIN;
    myHelper.set_root(rootFolderIN);
}


void two_particle_mc::runSimulation(double ss1, double ss2, unsigned long long num_steps,
                                    unsigned long long stride){
    
    unsigned long long num_steps_total = 0;
    
    SpringEnergy V_spring1(num_beads1,mass,beta/num_beads1);
    SpringEnergy V_spring2(num_beads2,mass,beta/num_beads2);
    two_particle_Hamiltonian H(num_beads1,num_beads2,V_spring1,V_spring2);
    two_particle_Estimator Esti(num_beads1,num_beads2,beta,V_spring1,V_spring2);

    vector<double> Q1(num_beads1,0);
    vector<double> Q_prop1(num_beads1,0);
    vector<double> Q2(num_beads2,0);
    vector<double> Q_prop2(num_beads2,0);

    gen_initQ(Q1,num_beads1,ss1);
    gen_initQ(Q2,num_beads2,ss2);
    
    if (readPSV1) {myHelper.read_PSV(num_beads1,Q1,"1");}
    if (readPSV2) {myHelper.read_PSV(num_beads2,Q2,"2");}

    Q_prop1 = Q1;
    Q_prop2 = Q2;
 
    double energy = H.get_energy(Q1,Q2);
    double energy_prop = energy;

    double estimator = Esti.get_estimator(Q1,Q2);
    double estimator_prop = estimator;
    double estimator_total = 0;
    int esti_samples = 0;
    vector<double> estimator_t(num_steps/stride,0);

    if (readData){myHelper.read_MC_data(estimator_total,num_steps_total);}

    
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

    for (unsigned long long step=0; step<num_steps; step++) {
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
                    estimator = Esti.get_estimator(Q1,Q2);
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
                    estimator = Esti.get_estimator(Q1,Q2);
                    sys2_steps_accpt += 1;
                }
                else{
                    Q_prop2(mcMove) = Q2(mcMove);
                }
                sys2_steps_total += 1;
            }
        }

        estimator_total += estimator;

        if(step % stride == 0){
            estimator_t[esti_samples] = estimator_total/(num_steps_total + step + 1);
            esti_samples += 1;
        }
    }

    myHelper.print_sys_accpt(sys1_steps_total,sys1_steps_accpt,"1");
    myHelper.print_sys_accpt(sys2_steps_total,sys2_steps_accpt,"2");
    myHelper.print_avg_energy(estimator_total, num_steps + num_steps_total);
    myHelper.write_estimator(estimator_t,stride);

    if(writePSV1){myHelper.write_PSV(num_beads1, Q1, "1");}
    if(writePSV2){myHelper.write_PSV(num_beads2, Q2, "2");}
    if (writeData){myHelper.write_MC_data(estimator_total,num_steps + num_steps_total);}
}

bool two_particle_mc::check_move(double energy, double energy_prop){

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
inline double two_particle_mc::nuc_dist(const double rn, double step_size){
    return (rn * 2.0 * step_size) - step_size;
}
inline int two_particle_mc::rand_bead(const Ullong rn, int num_beads){
    return rn % num_beads;
}
//
//
//void MonteCarlo_MTSastra::set_num_steps(unsigned long long num_steps_In){
//    num_steps = num_steps_In;
//}
//
//void MonteCarlo_MTSastra::set_esti_rate(unsigned long long esti_rate_In){
//    esti_rate = esti_rate_In;
//}
//
//void MonteCarlo_MTSastra::set_write_PSV(bool set_In){
//    writePSV = set_In;
//}
//
//void MonteCarlo_MTSastra::set_read_PSV(bool set_In){
//    readPSV = set_In;
//}
//
//void MonteCarlo_MTSastra::set_read_Data(bool set_In){
//    readData = set_In;
//}
//
//void MonteCarlo_MTSastra::set_write_Data(bool set_In){
//    writeData = set_In;
//}
