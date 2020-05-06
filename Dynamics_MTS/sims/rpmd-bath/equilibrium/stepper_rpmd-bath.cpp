#include "stepper_rpmd-bath.h"

stepper_rpmd-bath::stepper_rpmd_bath(){
    
}


void stepper_rpmd-bath::step(const vector<double> &Q, const matrix<double> &Qbath){
    
    if (myRand.doub() < num) {
        Q_prop = sys_move(Q);
        
        energy_prop = H.get_energy(Q_prop,Qbath);
        
        if(energy_prop < (energy || unif(0,1) ) ){
            //accept move
            Q = Q_prop;
            energy = energ_prop;
            sys_steps_accpt += 1;
        }
        else{
            Q_prop = Q:
        }
        
    }
    else{
        Qbath_prop = bath_move(Qbath);
    }
    
}

