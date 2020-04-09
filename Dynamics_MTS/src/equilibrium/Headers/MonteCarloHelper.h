#ifndef MonteCarloHelper_h
#define MonteCarloHelper_h

#include <iostream>
#include <fstream>

// Print system Monte Carlo acceptance ratios after calling runMC.
inline void print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt){
  std::cout << "\t System Acceptance Ratio: " << 100*(double)sys_steps_accpt/sys_steps << std::endl;
}

//* Print electronic Monte Carlo acceptance ratios after calling runMC.  */
inline void print_elec_accpt(unsigned long long elec_steps, unsigned long long elec_steps_accpt){
  std::cout << "\t Electronic Acceptance Ratio: " << 100*(double)elec_steps_accpt/elec_steps << std::endl;
}

/* Print the MV-RPMD average energy after calling runMC from MonteCarlo */
inline void print_avg_energy(double estimator_total, double sgn_total){
  std::cout << "\t Average Energy: " << estimator_total/sgn_total << std::endl;
}

/* Write estimator data to file "estimator" after calling runMC.  */
inline void write_estimator(vector <double> estimator,int interval){

  int size = estimator.size();

  std::ofstream myFile;
  myFile.open("/Users/ellioteklund/Desktop/Equilibrium_MTS_Alpha/MAVARIC_MTS/Results/estimator.txt");
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }

  for(int i=0; i<size; i++){
    myFile << i*interval << " " <<  estimator[i] << std::endl;
  }

  myFile.close();

  std::cout << "\t Successfully wrote energy_estimator file to Results." << std::endl;
}


inline void write_PSV(int nuc_beads,int elec_beads, int num_states,
                      vector<double> Q, matrix<double> x,matrix <double> p){
    
    std::ofstream myFile;
  myFile.open("/Users/ellioteklund/Desktop/Equilibrium_MTS_Alpha/MAVARIC_MTS/Results/PSV.txt");
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    for(int bead=0; bead<nuc_beads; bead++){
        myFile << Q(bead) << std::endl;
    }
    
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile << x(bead,state) << std::endl;
        }
    }
    
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile << p(bead,state) << std::endl;
        }
    }
    
    std::cout << "\t Successfully saved PSV to Results." << std::endl;
    myFile.close();
}

/* Read in PSV.*/
inline void read_PSV(int nuc_beads,int elec_beads, int num_states, vector<double> &Q, matrix<double> &x,
                     matrix<double> &p){
    
    std::ifstream myFile;
  myFile.open("/Users/ellioteklund/Desktop/Equilibrium_MTS_Alpha/MAVARIC_MTS/Results/PSV.txt");
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    for(int bead=0; bead<nuc_beads; bead++){
        myFile >> Q(bead);
    }
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile >> x(bead,state);
        }
    }
    
    for(int bead=0; bead<elec_beads; bead++){
        for (int state=0; state<num_states; state++){
            myFile >> p(bead,state);
        }
    }
    
    myFile.close();
    
    std::cout << "\t Successfully read PSV file from Results." << std::endl;
}

inline void write_MC_data(double sgn_total,double estimator_total){
    
    std::ofstream myFile;
  myFile.open("/Users/ellioteklund/Desktop/Equilibrium_MTS_Alpha/MAVARIC_MTS/Results/mcData.txt");
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile << sgn_total << std::endl;
    myFile << estimator_total << std::endl;
    
    myFile.close();
    
    std::cout << "\t Successfully saved MC data to Results." << std::endl;
    
}

inline void read_MC_data(double &sgn_totalGlobal, double &estimator_total){
    
    std::ifstream myFile;
  myFile.open("/Users/ellioteklund/Desktop/Equilibrium_MTS_Alpha/MAVARIC_MTS/Results/mcData.txt");
    
    if(!myFile.is_open()) {
        std::cout << "Could not open file" << std::endl;
    }
    
    myFile >> sgn_totalGlobal;
    myFile >> estimator_total;
    
    myFile.close();
    
    std::cout << "\t Successfully read MC datat from Results." << std::endl;
}


#endif /* MonteCarloHelper_h */
