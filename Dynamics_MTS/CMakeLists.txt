cmake_minimum_required(VERSION 3.10)

set(MACHINE "MAC")

set(ROOT_MAC "/Users/ellioteklund/Desktop/Dynamics_MTS_git/Dynamics_MTS")
set(ROOT_LIN "/home/elliot/Desktop")
set(ROOT_ASTRA "/home/fs01/ece52/MAVARIC-MTS/Dynamics_MTS")

set(BOOST_MAC "/Users/ellioteklund/Downloads")
set(BOOST_LIN "/home/elliot")
set(BOOST_ASTRA "/home/fs01/ece52/MAVARIC_MTS")

# set the project name
project(mavaric)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_COMPILER mpic++)

if(${MACHINE} STREQUAL "MAC")
    set(BOOST_DIR "${BOOST_MAC}/boost_1_72_0")
    set(DYNAMICS_DIR "${ROOT_MAC}/src/dynamics")
    set(EQUILIBRIUM_DIR "${ROOT_MAC}/src/equilibrium")
    set(FORCES_DIR "${ROOT_MAC}/src/forces")
    set(HAMILTONIAN_DIR "${ROOT_MAC}/src/hamiltonian")
    set(MISC_DIR "${ROOT_MAC}/src/misc")

elseif(${MACHINE} STREQUAL "LIN")
    set(INCLUDE_DIRS "${BOOST_LIN}/boost_1_72_0")
    set(MY_HEADERS "${ROOT_LIN}/Dynamics_MTS_AlphaDEBUG/Dynamics_MTS/DynamicsParts/Headers")
    set(MY_SOURCE "${ROOT_LIN}/Dynamics_MTS_AlphaDEBUG/Dynamics_MTS/DynamicsParts/Sources")

elseif(${MACHINE} STREQUAL "ASTRA")
    set(INCLUDE_DIRS "${BOOST_ASTRA}/boost_1_72_0")
    set(MY_HEADERS "${ROOT_ASTRA}/DynamicsParts/Headers")
    set(MY_SOURCE "${ROOT_ASTRA}/DynamicsParts/Sources")

endif()

add_definitions( "-DNDEBUG -O3" )
set(EXECUTABLE_OUTPUT_PATH ${ROOT_MAC})

# add the executable
add_executable(mavaric main.cpp

    #dynamics directory
    ${DYNAMICS_DIR}/ABM_MVRPMD.cpp
    ${DYNAMICS_DIR}/Dynamics.cpp
    ${DYNAMICS_DIR}/RK4_MVRPMD.cpp

    #equilibrium directory
    ${EQUILIBRIUM_DIR}/dM_Matrix_dBeta.cpp
    ${EQUILIBRIUM_DIR}/dM_Matrix_MTS_dBeta.cpp
    ${EQUILIBRIUM_DIR}/dTheta_dBeta.cpp
    ${EQUILIBRIUM_DIR}/dTheta_MTS_dBeta.cpp
    ${EQUILIBRIUM_DIR}/MonteCarlo_MTS.cpp
    ${EQUILIBRIUM_DIR}/MonteCarlo.cpp
    ${EQUILIBRIUM_DIR}/MonteCarloHelper.cpp
    ${EQUILIBRIUM_DIR}/MVRPMD_Estimator.cpp
    ${EQUILIBRIUM_DIR}/MVRPMD_MTS_Estimator.cpp
    ${EQUILIBRIUM_DIR}/Sampling_MTS.cpp

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

    #hamiltonian directory
    ${HAMILTONIAN_DIR}/BathSpringEnergy.cpp
    ${HAMILTONIAN_DIR}/C_Matrix.cpp
    ${HAMILTONIAN_DIR}/CouplingEnergy.cpp
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


    #misc directory
    ${MISC_DIR}/MainHlpr.cpp
    ${MISC_DIR}/PopulationEstimator.cpp)


include_directories(${BOOST_DIR})
include_directories(${DYNAMICS_DIR}/Headers)
include_directories(${EQUILIBRIUM_DIR}/Headers)
include_directories(${FORCES_DIR}/Headers)
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)

