#include "PopulationEstimator.hpp"

PopulationEstimator::PopulationEstimator(int elec_beads,int num_states)
    :elec_beads(elec_beads), num_states(num_states),
     xx(elec_beads,num_states,0), pp(elec_beads,num_states,0),
     g(elec_beads,num_states,0), g_sum(elec_beads,0),
     populations(num_states,0),
     ones(num_states,1.0),
     coeff(pow(2.0,num_states+1)/elec_beads)

{}

void PopulationEstimator::update_g(const matrix<double> &x, const matrix<double> &p){
    
    xx = element_prod(x,x);
    pp = element_prod(p,p);

    noalias(g) = xx + pp;
}

void PopulationEstimator::update_populations(const matrix<double> &x,const matrix<double> &p){
  
    update_g(x,p);
    g_sum = prod(g, ones);
    
    double sum = 0;
    for (int state=0; state <num_states; state++) {
        sum = 0;
        for (int bead=0; bead<elec_beads; bead++) {
            sum += exp(-g_sum(bead))*(xx(bead,state) + pp(bead,state) - 0.5);
        }
        populations(state) = coeff * sum;
    }
}

double PopulationEstimator::get_pop(const int state){
    return populations(state);
}

vector<double> PopulationEstimator::get_pop(){
    return populations;
}

void PopulationEstimator::write_populations(matrix<double> final_pops, double dt, int data_steps, int rate, std::string root){
    
    std::string fileName = root + "/Results/PopAC";
    std::ofstream myFile;
    myFile.open(fileName.c_str());
    
    if (!myFile.is_open()) {
        std::cout << "ERROR: Could not open file " << fileName << std::endl;
    }
    
    
    for (int i=0; i<data_steps; i++) {
        myFile << i*dt*rate << " ";
        
        for (int state=0; state<num_states; state++) {
            myFile << final_pops(i,state) << " ";
        }
        myFile << std::endl;
    }
    
    myFile.close();
    
}
