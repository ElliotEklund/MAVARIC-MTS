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
    set(EQUILIBRIUM_DIR "${ROOT_MAC}/sims/two-particle/equilibrium")
    set(SAMPLING_DIR "${ROOT_MAC}/sims/two-particle/sampling")
    set(DYNAMICS_DIR "${ROOT_MAC}/sims/two-particle/dynamics")
    set(HAMILTONIAN_DIR "${ROOT_MAC}/src/hamiltonian")
    set(FORCES_DIR "${ROOT_MAC}/src/forces")
    set(ESTIMATORS_DIR "${ROOT_MAC}/src/estimators")
    set(INTEGRATORS_DIR "${ROOT_MAC}/src/integrators")
    set(MISC_DIR "${ROOT_MAC}/src/misc")
    set(CMAKE_CXX_FLAGS "-DMAC")
endif()

add_definitions( "-DNDEBUG -O3" )

# add the executable

if(${MACHINE} STREQUAL "MAC")
    add_executable(two-particle two-particle.cpp

        #equilibrium directory
        ${EQUILIBRIUM_DIR}/two_particle_mc.cpp
        ${EQUILIBRIUM_DIR}/MonteCarloHelper.cpp
        ${EQUILIBRIUM_DIR}/input.cpp

        #sampling directory
        ${SAMPLING_DIR}/two_particle_sampling.cpp
        ${SAMPLING_DIR}/SamplingHelper.cpp

        #dynamics directory
        ${DYNAMICS_DIR}/two_particle_Dynamics.cpp

        #hamiltonian directory
        ${HAMILTONIAN_DIR}/two_particle_Hamiltonian.cpp
        ${HAMILTONIAN_DIR}/SpringEnergy.cpp

        #estimator directory
        ${ESTIMATORS_DIR}/two_particle_Estimator.cpp

        #forces directory
        ${FORCES_DIR}/Forces_two_particles.cpp
        ${FORCES_DIR}/dSpring_dQ.cpp

        #integrators directory
        ${INTEGRATORS_DIR}/VV_two_particle.cpp

        #misc directory
        ${MISC_DIR}/MainHlpr.cpp
        ${MISC_DIR}/pop_estimators.cpp
        ${MISC_DIR}/aggregate.cpp
        ${MISC_DIR}/PopulationEstimator.cpp)
endif()

include_directories(${BOOST_DIR})
include_directories(${EQUILIBRIUM_DIR}/Headers)
include_directories(${SAMPLING_DIR}/Headers)
include_directories(${DYNAMICS_DIR}/Headers)
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${FORCES_DIR}/Headers)
include_directories(${ESTIMATORS_DIR}/Headers)
include_directories(${INTEGRATORS_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)

