cmake_minimum_required(VERSION 3.10)

set(MACHINE "MAC")

set(ROOT_MAC "/Users/ellioteklund/Desktop/MAVARIC_v2.0/MAVARIC")

set(BOOST_MAC "/Users/ellioteklund/Downloads")

# set the project name
project(mavaric)
set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_COMPILER mpic++)

if(${MACHINE} STREQUAL "MAC")
    set(BOOST_DIR "${BOOST_MAC}/boost_1_72_0")
    set(HAMILTONIAN_DIR "${ROOT_MAC}/src/hamiltonian")
    set(ESTIMATORS_DIR "${ROOT_MAC}/src/estimators")
    set(MISC_DIR "${ROOT_MAC}/src/misc")
    set(MONTECARLO_DIR "${ROOT_MAC}/src/monte_carlo")
    set(FORCES_DIR "${ROOT_MAC}/src/forces")
    set(INTEGRATORS_DIR "${ROOT_MAC}/src/integrators")
    set(CMAKE_CXX_FLAGS "-DMAC")
endif()

add_definitions( "-DNDEBUG -O3" )

# add the executable

if(${MACHINE} STREQUAL "MAC")
    add_executable(test test.cpp

    #hamiltonian directory
    ${HAMILTONIAN_DIR}/C_Matrix.cpp
    ${HAMILTONIAN_DIR}/M_Matrix_MTS.cpp
    ${HAMILTONIAN_DIR}/M_Matrix.cpp
    ${HAMILTONIAN_DIR}/StateDepPots.cpp
    ${HAMILTONIAN_DIR}/StateIndepPot.cpp
    ${HAMILTONIAN_DIR}/SpringEnergy.cpp
    ${HAMILTONIAN_DIR}/GTerm.cpp
    ${HAMILTONIAN_DIR}/Theta_MTS.cpp
    ${HAMILTONIAN_DIR}/Theta.cpp
    ${HAMILTONIAN_DIR}/theta_mixed.cpp
    ${HAMILTONIAN_DIR}/MVRPMD_MTS_Hamiltonian.cpp
    ${HAMILTONIAN_DIR}/mvrpmd_ham_mixed.cpp

    #estimators directory
    ${ESTIMATORS_DIR}/dM_Matrix_dBeta.cpp
    ${ESTIMATORS_DIR}/theta_mixed_dBeta.cpp
    ${ESTIMATORS_DIR}/mvrpmd_mixed_esti.cpp
    ${ESTIMATORS_DIR}/dM_Matrix_MTS_dBeta.cpp
    ${ESTIMATORS_DIR}/dTheta_dBeta.cpp
    ${ESTIMATORS_DIR}/dTheta_MTS_dBeta.cpp
    ${ESTIMATORS_DIR}/MVRPMD_MTS_Estimator.cpp

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
    ${FORCES_DIR}/Forces_MTS.cpp
    ${FORCES_DIR}/theta_mixed_dElec.cpp

    #monte Carlo directory
    ${MONTECARLO_DIR}/system_step.cpp

    #integrators directory
    ${INTEGRATORS_DIR}/RK4_MVRPMD.cpp
    ${INTEGRATORS_DIR}/ABM_MVRPMD.cpp

    #misc directory
    ${MISC_DIR}/MainHlpr.cpp
    ${MISC_DIR}/pop_estimators.cpp
    ${MISC_DIR}/aggregate.cpp
    ${MISC_DIR}/PopulationEstimator.cpp)
endif()

include_directories(${BOOST_DIR})
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${ESTIMATORS_DIR}/Headers)
include_directories(${FORCES_DIR}/Headers)
include_directories(${INTEGRATORS_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)
include_directories(${MONTECARLO_DIR}/Headers)

