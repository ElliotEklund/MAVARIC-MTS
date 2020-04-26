#include "Sampling_MTS.hpp"

Sampling_MTS::Sampling_MTS(int my_id, int num_procs, int root_proc,int num_beads, int elec_beads,double mass,
                           int num_states,double beta_num_beads, double beta_elec_beads, double nuc_ss,
                           double elec_ss, std::string root)
    :/* Initialize parameters*/
     root(root),
     my_id(my_id), num_procs(num_procs), root_proc(root_proc),

     num_beads(num_beads), num_states(num_states), mass(mass),
     beta_num_beads(beta_num_beads), elec_beads(elec_beads),

     /* Initialize phase space variable vectors*/
     Q(num_beads), x(elec_beads,num_states), p(elec_beads,num_states),
     Q_prop(num_beads), x_prop(elec_beads,num_states), p_prop(elec_beads,num_states),
    
     /* Initialize random number generators */
     mt(myRand()), nuc_dist(-nuc_ss,nuc_ss), elec_dist(-elec_ss,elec_ss),
     unif_1(0.0,1.0), rand_bead(0,num_beads-1), rand_elec_bead(0,elec_beads-1),

     /* Initialize Hamiltonian parts and Estimator*/
     V_spring(num_beads, mass, beta_num_beads),
     V0(num_beads, mass),
     G(elec_beads, num_states),
     C(elec_beads,num_states),
     M(num_states,num_beads,beta_elec_beads),
     M_MTS(num_beads,elec_beads,num_states,M),
     thetaMTS(num_states, elec_beads,C,M_MTS),
     H_MTS(beta_num_beads,V_spring,V0,G,thetaMTS)
{}

void Sampling_MTS::runSimulation(){

    readPSV();

    Q_prop = Q;
    x_prop = x;
    p_prop = p;
    
    Q_samp.resize(num_beads*num_samples);
    x_samp.resize(elec_beads*num_samples,num_states);
    p_samp.resize(elec_beads*num_samples,num_states);

    sys_steps_accpt = 0;
    sys_steps = 0;

    elec_steps_accpt = 0;
    elec_steps = 0;

    energy = H_MTS.get_energy(Q,x,p);
        
    for (int sample=0; sample<num_samples; sample++) {
        for (int step=0; step<decor_len; step++) {
            for (int m_type=0; m_type<2; m_type++) {
                if(unif_1(mt) > 0.5){
                    //Sample nuclear coordinates
                    sample_nuc();
                }
                else{
                    //Sample electronic coordinates
                    sample_elec();
                }
            }
        }

        stash_Q(sample);
        stash_x(sample);
        stash_p(sample);
    }
    
    saveQ();
    savex();
    savep();
    saveP();

    if (my_id==root_proc) {
        print_sys_accpt(sys_steps, sys_steps_accpt);
        print_elec_accpt(elec_steps, elec_steps_accpt);
    }

//    std::string fileName = root + "/Results/hist.txt";
//
//    int bins = 300;
//    vector<double> Q_cent (num_samples);
//    double sum = 0;
//    for(int traj=0; traj<num_samples; traj++){
//        sum = 0;
//        for(int bead=0; bead<num_beads; bead++){
//            sum += Q_samp(traj*num_beads + bead);
//        }
//        Q_cent(traj) = sum/num_beads;
//    }
//    histogram(fileName,bins,Q_cent);
}

void Sampling_MTS::sample_nuc(){

    int mcMove = 0;
    
    /* Propose new nuclear moves.*/
    for (int i=0; i<num_beads; i++) {
        mcMove = rand_bead(mt);
        Q_prop(mcMove) = Q(mcMove) + nuc_dist(mt);
        
        double energ_prop = H_MTS.get_energy(Q_prop,x,p);
        
        /* Accept new system moves if energ_prop < energy*/
        if(energ_prop < energy){
            Q = Q_prop;
            energy = energ_prop;
            sys_steps_accpt += 1;
        }
        
        /* Accept new system moves if inequality is met*/
        else if (unif_1(mt) <= exp(-beta_num_beads * (energ_prop - energy))){
            Q = Q_prop;
            energy = energ_prop;
            sys_steps_accpt += 1;
        }
        else{Q_prop = Q;}
    }

    sys_steps += 1;
}

void Sampling_MTS::sample_elec(){

    int mcMove = 0;
    
    //Propose new electronic moves.
    for (int bead=0; bead<elec_beads; bead++) {
        mcMove = rand_elec_bead(mt);
        
        for(int state=0; state<num_states; state++){
            x_prop(mcMove,state) = x(mcMove,state) + elec_dist(mt);
            p_prop(mcMove,state) = p(mcMove,state) + elec_dist(mt);
        }
        
        double energ_prop = H_MTS.get_energy(Q,x_prop,p_prop);
        
        /* Accept new electronic move if energ_prop < energy */
        if(energ_prop < energy){
            x = x_prop;
            p = p_prop;
            energy = energ_prop;
            elec_steps_accpt += 1;
        }
        
        /* Accept new system moves if inequality is met*/
        else if (unif_1(mt) <= exp(-beta_num_beads * (energ_prop - energy))){
            x = x_prop;
            p = p_prop;
            energy = energ_prop;
            elec_steps_accpt += 1;
        }
        
        else{
            x_prop = x;
            p_prop = p;
        }
    }
    
    elec_steps += 1;
}

void Sampling_MTS::set_decor_len(unsigned long long set_In){
    decor_len = set_In;
}

void Sampling_MTS::set_num_samples(int set_In){
    num_samples = set_In;
}

void Sampling_MTS::readPSV(){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "/Results/PSV" + quickConvert.str();
    
    std::ifstream myFile;
    myFile.open(fileName);
    
    if(!myFile.is_open()) {
          std::cout << "Could not open file" << std::endl;
      }
    
    for (int bead=0; bead<num_beads; bead++) {
        myFile >> Q(bead);
    }
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            myFile >> x(bead,state);
        }
    }
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            myFile >> p(bead,state);
        }
    }
}

// Print system Monte Carlo acceptance ratios after calling runMC.
void Sampling_MTS::print_sys_accpt(unsigned long long sys_steps,unsigned long long sys_steps_accpt){
  std::cout << "\t System Acceptance Ratio: " << 100*(double)sys_steps_accpt/sys_steps << std::endl;
}

//* Print electronic Monte Carlo acceptance ratios after calling runMC.  */
void Sampling_MTS::print_elec_accpt(unsigned long long elec_steps, unsigned long long elec_steps_accpt){
  std::cout << "\t Electronic Acceptance Ratio: " << 100*(double)elec_steps_accpt/elec_steps << std::endl;
}

void Sampling_MTS::stash_Q(int sample){
    
    int sample_stride = sample*num_beads;
    
    for (int bead=0; bead<num_beads; bead++) {
        Q_samp(sample_stride + bead) = Q(bead);
    }
    
}

void Sampling_MTS::stash_x(int sample){
    
    int sample_stride = sample*elec_beads;
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            x_samp(sample_stride+bead,state) = x(bead,state);
        }
    }
}

void Sampling_MTS::stash_p(int sample){
    
    int sample_stride = sample*elec_beads;
    
    for (int bead=0; bead<elec_beads; bead++) {
        for (int state=0; state<num_states; state++) {
            p_samp(sample_stride+bead,state) = p(bead,state);
        }
    }
}

void Sampling_MTS::saveQ(){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "/Results/Trajectories/Q" + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName);

    for (int sample=0; sample<num_samples; sample++) {
        for (int bead=0; bead<num_beads; bead++) {
            myFile << Q_samp(sample*num_beads + bead) << std::endl;
        }
    }

    myFile.close();
}

void Sampling_MTS::savex(){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "/Results/Trajectories/xelec" + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName);

    for (int sample=0; sample<num_samples; sample++) {
        for (int bead=0; bead<elec_beads; bead++) {
            for (int state=0; state<num_states; state++) {
                myFile << x_samp(sample*elec_beads + bead,state) << std::endl;
            }
        }
    }

    myFile.close();
}

void Sampling_MTS::savep(){
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "/Results/Trajectories/pelec" + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName);

    for (int sample=0; sample<num_samples; sample++) {
        for (int bead=0; bead<elec_beads; bead++) {
            for (int state=0; state<num_states; state++) {
                myFile << p_samp(sample*elec_beads + bead,state) << std::endl;
            }
        }
    }

    myFile.close();
}

void Sampling_MTS::saveP(){
    
    double stdev = sqrt(mass/beta_num_beads);
    std::normal_distribution<double> boltz_dist(0,stdev);
    
    std::ostringstream quickConvert;
    quickConvert << my_id;
    
    std::string fileName = root + "/Results/Trajectories/P" + quickConvert.str();
    
    std::ofstream myFile;
    myFile.open(fileName);
    
    for(int i=0; i<num_samples*num_beads; i++){
        myFile << boltz_dist(mt) << std::endl;
    }
    
    myFile.close();
}

template <typename T>
void Sampling_MTS::histogram(std::string fileName, int bins, T X){

  std::vector<double> ordered_X;
  int X_size = X.size();

  for(int i=0; i<X_size;i++){
    ordered_X.push_back(X[i]);
  }

  sort(ordered_X.begin(),ordered_X.end());

  double min = ordered_X[0];
  double max = ordered_X[X_size-1];
  double width = (max - min)/bins;

  double average = 0;

  for(int i=0; i<X_size; i++){
    average = average + ordered_X[i];
  }

  average = average/(X_size);

  double std = 0;

  for(int i=0; i<X_size; i++){
    std = std + (ordered_X[i] - average)*(ordered_X[i] - average);
  }

  std = sqrt(std/(X_size-1));

  double skew = 0;
  for(int i=0; i<X_size; i++){
    skew = skew + (ordered_X[i] - average)*(ordered_X[i] - average)
      *(ordered_X[i] - average)/(std*std*std);
  }

  skew = skew/(X_size);

  std::vector<double> finalHistogram (bins,0);

  int j = 1;
  for(int i=0; i<X_size; i++){
    if(ordered_X[i] <= (min + j*width)){
      finalHistogram[j-1] = finalHistogram[j-1] + 1;
    }
    else{j = j+1;}

  }
  std::fstream histogramOut;
  histogramOut.open(fileName.c_str());

  if(!histogramOut.is_open()){
      std::cout << "Failed to open histogram file." << std::endl;
  }
  else { /* no statement */ }

  double mode = finalHistogram[0];
  double largestValue = finalHistogram[0];

  for(int i=0; i<bins; i++){
    histogramOut << double(min + i*width) << " " << finalHistogram[i]/(X_size) << std::endl;

    if(finalHistogram[i] > largestValue){
      largestValue = finalHistogram[i];
      mode = (min + i*width);
    }
    else { /* no statement */ }
  }

  histogramOut.close();

  /* Print Statistics about the distribution. */

  std::cout << std::endl;
  std::cout << "\t--- Histogram Statistics ---" << std::endl;
  std::cout << "\t    Mode: " << mode << std::endl;
  std::cout << "\t    Average: " << average << std::endl;
  std::cout << "\t    Standard Deviation: " << std << std::endl;
  std::cout << "\t    Skew: " << skew << std::endl << std::endl;
}

