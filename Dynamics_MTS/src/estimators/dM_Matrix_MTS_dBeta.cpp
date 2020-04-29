#include "dM_Matrix_MTS_dBeta.hpp"

dM_Matrix_MTS_dBeta::dM_Matrix_MTS_dBeta(int nuc_beads, int elec_beads, int num_states,
                                         double beta_nuc_beads, M_Matrix &M_In)
    :dM_dBeta(nuc_beads,elec_beads,num_states, beta_nuc_beads, M_In),
     dM_MTS_dBeta_vec(elec_beads), elec_beads(elec_beads), num_states(num_states),
     gamma(nuc_beads/elec_beads)
{}

void dM_Matrix_MTS_dBeta::update_dM_MTS_dBeta_vec(){
    
    dM_dBeta.update_dM_dBeta_vec();
    matrix<std::complex<double> > temp_sum (num_states,num_states);
      
      for (int bead=0; bead<elec_beads; bead++) {
          std::fill(temp_sum.data().begin(), temp_sum.data().end(), 0);
          for (int slice=0; slice<gamma; slice++) {
              temp_sum += dM_dBeta.get_dM_dBeta_alpha(gamma*bead + slice);
          }
          dM_MTS_dBeta_vec(bead) = (1.0/gamma) * temp_sum;
      }
}

matrix<std::complex<double> >& dM_Matrix_MTS_dBeta::get_dM_MTS_dBeta_alpha(int alpha){
    return dM_MTS_dBeta_vec(alpha);
}
