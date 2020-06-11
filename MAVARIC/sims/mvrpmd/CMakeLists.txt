cmake_minimum_required(VERSION 3.10)

set(MACHINE "MAC")

if(${MACHINE} STREQUAL "MAC")
  set(ROOT "/Users/ellioteklund/Desktop/MAVARIC_v2.0/MAVARIC")
  set(BOOST "/Users/ellioteklund/Downloads")
elseif(${MACHINE} STREQUAL "ASTRA")
  set(ROOT "/home/fs01/ece52/MAVARIC-MTS/MAVARIC")
  set(BOOST "/home/fs01/ece52/MAVARIC_MTS")
endif()

# set the project name
project(mavaric)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_COMPILER mpic++)

#if(${MACHINE} STREQUAL "MAC")
    set(BOOST_DIR "${BOOST}/boost_1_72_0")
    set(EQUILIBRIUM_DIR "${ROOT}/sims/mvrpmd/equilibrium")
    set(SAMPLING_DIR "${ROOT}/sims/mvrpmd/sampling")
    set(DYNAMICS_DIR "${ROOT}/sims/mvrpmd/dynamics")
    set(HAMILTONIAN_DIR "${ROOT}/src/hamiltonian")
    set(ESTIMATORS_DIR "${ROOT}/src/estimators")
    set(FORCES_DIR "${ROOT}/src/forces")
    set(INTEGRATORS_DIR "${ROOT}/src/integrators")
    set(MISC_DIR "${ROOT}/src/misc")
    set(MONTECARLO_DIR "${ROOT}/src/monte_carlo")
    set(EXECUTABLE_OUTPUT_PATH "${ROOT}/bin")
    set(CMAKE_CXX_FLAGS "-DMAC")

    #elseif(${MACHINE} STREQUAL "ASTRA")
    #    set(BOOST_DIR "${BOOST_ASTRA}/boost_1_72_0")
    #    set(DYNAMICS_DIR "${ROOT_ASTRA}/sims/mvrpmd/dynamics")
    #    set(SAMPLING_DIR "${ROOT_ASTRA}/sims/mvrpmd/sampling")
    #    set(EQUILIBRIUM_DIR "${ROOT_ASTRA}/sims/mvrpmd/equilibrium")
    #    set(FORCES_DIR "${ROOT_ASTRA}/src/forces")
    #    set(HAMILTONIAN_DIR "${ROOT_ASTRA}/src/hamiltonian")
    #    set(ESTIMATORS_DIR "${ROOT_ASTRA}/src/estimators")
    #    set(INTEGRATORS_DIR "${ROOT_ASTRA}/src/integrators")
    #    set(MISC_DIR "${ROOT_ASTRA}/src/misc")
    #    set(MONTECARLO_DIR "${ROOT_ASTRA}/src/monte_carlo")
    #    set(EXECUTABLE_OUTPUT_PATH "${ROOT_ASTRA}/bin")
    #    set(CMAKE_CXX_FLAGS "-ASTRA")
    #endif()


add_definitions( "-DNDEBUG -O3" )

# add the executable

add_executable(mvrpmd mvrpmd.cpp

    #equilibrium directory
    ${EQUILIBRIUM_DIR}/MonteCarloHelper.cpp
    ${EQUILIBRIUM_DIR}/equilib_mvrpmd.cpp
    ${EQUILIBRIUM_DIR}/input_mvrpmd.cpp

    #sampling directory
    #${SAMPLING_DIR}/Sampling_MTS.cpp
    ${SAMPLING_DIR}/sampling_mvrpmd.cpp
    ${SAMPLING_DIR}/SamplingHelper.cpp

    #dynamics directory
    ${DYNAMICS_DIR}/Dynamics.cpp

    #hamiltonian directory
    ${HAMILTONIAN_DIR}/C_Matrix.cpp
    ${HAMILTONIAN_DIR}/GTerm.cpp
    ${HAMILTONIAN_DIR}/M_Matrix_MTS.cpp
    ${HAMILTONIAN_DIR}/M_Matrix.cpp
    ${HAMILTONIAN_DIR}/MVRPMD_Hamiltonian.cpp
    ${HAMILTONIAN_DIR}/MVRPMD_MTS_Hamiltonian.cpp
    ${HAMILTONIAN_DIR}/SpringEnergy.cpp
    ${HAMILTONIAN_DIR}/StateDepPots.cpp
    ${HAMILTONIAN_DIR}/StateIndepPot.cpp
    ${HAMILTONIAN_DIR}/Theta_MTS.cpp
    ${HAMILTONIAN_DIR}/Theta.cpp
    ${HAMILTONIAN_DIR}/theta_mixed.cpp
    ${HAMILTONIAN_DIR}/mvrpmd_mixed_ham.cpp

    #estimators directory
    ${ESTIMATORS_DIR}/dM_Matrix_dBeta.cpp
    ${ESTIMATORS_DIR}/dM_Matrix_MTS_dBeta.cpp
    ${ESTIMATORS_DIR}/dTheta_dBeta.cpp
    ${ESTIMATORS_DIR}/dTheta_MTS_dBeta.cpp
    ${ESTIMATORS_DIR}/MVRPMD_Estimator.cpp
    ${ESTIMATORS_DIR}/MVRPMD_MTS_Estimator.cpp
    ${ESTIMATORS_DIR}/theta_mixed_dBeta.cpp
    ${ESTIMATORS_DIR}/mvrpmd_mixed_esti.cpp

    #forces directory
    ${FORCES_DIR}/dCdelec.cpp
    ${FORCES_DIR}/dM_Matrix_dQ.cpp
    ${FORCES_DIR}/dM_Matrix_MTS_dQ.cpp
    ${FORCES_DIR}/dSpring_dQ.cpp
    ${FORCES_DIR}/dStateDep_dQ.cpp
    ${FORCES_DIR}/dStateIndep_dQ.cpp
    ${FORCES_DIR}/dTheta_MTS_dElec.cpp
    ${FORCES_DIR}/dTheta_MTS_dQ.cpp
    ${FORCES_DIR}/Forces_MTS.cpp
    ${FORCES_DIR}/theta_mixed_dQ.cpp
    ${FORCES_DIR}/mvrpmd_mixed_forces.cpp
    ${FORCES_DIR}/theta_mixed_dElec.cpp

    #integrators directory
    ${INTEGRATORS_DIR}/ABM_MVRPMD.cpp
    ${INTEGRATORS_DIR}/RK4_MVRPMD.cpp

    #misc directory
    ${MISC_DIR}/MainHlpr.cpp
    ${MISC_DIR}/mpi_wrapper.cpp
    ${MISC_DIR}/PopulationEstimator.cpp
    ${MISC_DIR}/decorrelation.cpp
    ${MISC_DIR}/pop_estimators.cpp
    ${MISC_DIR}/aggregate.cpp
    ${MISC_DIR}/init_PAC.cpp

    #monte carlo directory
    ${MONTECARLO_DIR}/system_step.cpp
    ${MONTECARLO_DIR}/elec_step.cpp)

include_directories(${BOOST_DIR})
include_directories(${EQUILIBRIUM_DIR}/Headers)
include_directories(${SAMPLING_DIR}/Headers)
include_directories(${DYNAMICS_DIR}/Headers)
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${ESTIMATORS_DIR}/Headers)
include_directories(${FORCES_DIR}/Headers)
include_directories(${INTEGRATORS_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)
include_directories(${MONTECARLO_DIR}/Headers)

