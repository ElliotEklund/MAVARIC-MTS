cmake_minimum_required(VERSION 3.10)

set(MACHINE "MAC")

set(ROOT_MAC "/Users/ellioteklund/Desktop/Dynamics_MTS_git/Dynamics_MTS")
set(BOOST_MAC "/Users/ellioteklund/Downloads")

# set the project name
project(mavaric)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_COMPILER mpic++)

if(${MACHINE} STREQUAL "MAC")
    set(BOOST_DIR "${BOOST_MAC}/boost_1_72_0")
    set(EQUILIBRIUM_DIR "${ROOT_MAC}/sims/rpmd-bath/equilibrium")
    #set(SAMPLING_DIR "${ROOT_MAC}/sims/mvrpmd/sampling")
    #set(DYNAMICS_DIR "${ROOT_MAC}/sims/mvrpmd/dynamics")
    #set(HAMILTONIAN_DIR "${ROOT_MAC}/src/hamiltonian")
    #set(ESTIMATORS_DIR "${ROOT_MAC}/src/estimators")
    #set(FORCES_DIR "${ROOT_MAC}/src/forces")
    #set(INTEGRATORS_DIR "${ROOT_MAC}/src/integrators")
    set(MISC_DIR "${ROOT_MAC}/src/misc")
    set(EXECUTABLE_OUTPUT_PATH "${ROOT_MAC}/bin")
    set(CMAKE_CXX_FLAGS "-DMAC")
endif()

add_definitions( "-DNDEBUG -O3" )

# add the executable

if(${MACHINE} STREQUAL "MAC")
    add_executable(rpmd-bath rpmd-bath.cpp

        #equilibrium directory
        # ${EQUILIBRIUM_DIR}/MonteCarlo.cpp
        # ${EQUILIBRIUM_DIR}/MonteCarloHelper.cpp
        # ${EQUILIBRIUM_DIR}/MonteCarlo_MTS.cpp
        ${EQUILIBRIUM_DIR}/move_system.cpp
        ${EQUILIBRIUM_DIR}/move_bath.cpp
        ${EQUILIBRIUM_DIR}/input_rpmd_bath.cpp

        # #sampling directory
        # ${SAMPLING_DIR}/Sampling_MTS.cpp

        # #dynamics directory
        # ${DYNAMICS_DIR}/Dynamics.cpp

        # #hamiltonian directory
        # ${HAMILTONIAN_DIR}/C_Matrix.cpp
        # ${HAMILTONIAN_DIR}/GTerm.cpp
        # ${HAMILTONIAN_DIR}/M_Matrix_MTS.cpp
        # ${HAMILTONIAN_DIR}/M_Matrix.cpp
        # ${HAMILTONIAN_DIR}/MVRPMD_Hamiltonian.cpp
        # ${HAMILTONIAN_DIR}/MVRPMD_MTS_Hamiltonian.cpp
        # ${HAMILTONIAN_DIR}/SpringEnergy.cpp
        # ${HAMILTONIAN_DIR}/StateDepPots.cpp
        # ${HAMILTONIAN_DIR}/StateIndepPot.cpp
        # ${HAMILTONIAN_DIR}/Theta_MTS.cpp
        # ${HAMILTONIAN_DIR}/Theta.cpp

        # #estimators directory
        # ${ESTIMATORS_DIR}/dM_Matrix_dBeta.cpp
        # ${ESTIMATORS_DIR}/dM_Matrix_MTS_dBeta.cpp
        # ${ESTIMATORS_DIR}/dTheta_dBeta.cpp
        # ${ESTIMATORS_DIR}/dTheta_MTS_dBeta.cpp
        # ${ESTIMATORS_DIR}/MVRPMD_Estimator.cpp
        # ${ESTIMATORS_DIR}/MVRPMD_MTS_Estimator.cpp

        # #forces directory
        # ${FORCES_DIR}/dCdelec.cpp
        # ${FORCES_DIR}/dM_Matrix_dQ.cpp
        # ${FORCES_DIR}/dM_Matrix_MTS_dQ.cpp
        # ${FORCES_DIR}/dSpring_dQ.cpp
        # ${FORCES_DIR}/dStateDep_dQ.cpp
        # ${FORCES_DIR}/dStateIndep_dQ.cpp
        # ${FORCES_DIR}/dTheta_MTS_dElec.cpp
        # ${FORCES_DIR}/dTheta_MTS_dQ.cpp
        # ${FORCES_DIR}/Forces_MTS.cpp

        # #integrators directory
        # ${INTEGRATORS_DIR}/ABM_MVRPMD.cpp
        # ${INTEGRATORS_DIR}/RK4_MVRPMD.cpp

        # #misc directory
        ${MISC_DIR}/aggregate.cpp
        # ${MISC_DIR}/MainHlpr.cpp
        # ${MISC_DIR}/PopulationEstimator.cpp
        # ${MISC_DIR}/decorrelation.cpp
        # ${MISC_DIR}/init_PAC.cpp)
                                )
endif()

include_directories(${BOOST_DIR})
include_directories(${EQUILIBRIUM_DIR}/Headers)
#include_directories(${SAMPLING_DIR}/Headers)
#include_directories(${DYNAMICS_DIR}/Headers)
#include_directories(${HAMILTONIAN_DIR}/Headers)
#include_directories(${ESTIMATORS_DIR}/Headers)
#include_directories(${FORCES_DIR}/Headers)
#include_directories(${INTEGRATORS_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)

