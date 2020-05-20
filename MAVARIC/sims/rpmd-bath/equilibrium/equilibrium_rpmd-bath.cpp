


equilibrium_rpmd::equilibrium_rpmd(int my_id, int root_proc, int num_procs){
    :/* Initialize parameters*/
    my_id(my_id), root_proc(root_proc), num_procs(num_procs),
    
    /* Initialize Hamiltonian parts and Estimator*/
    V_spring(num_beads, mass, beta_num_beads),
    V0(num_beads, mass),
    //TODO: Add all parts for hamiltonian
    

    H(beta_num_beads,V_spring,V0,G,thetaMTS),

    
}
 
void equilibrium_rpmd-bath::compute(){
    
    
    for (int step=0; step<num_steps; step++){
        
        my_stepper.step(Q,Qbath);
        estimator = Esti.get_estimator(Q,Qbath);
        
        //call stepper
        //compute estimator
        //compute requested calculations
        
    }
    
}

