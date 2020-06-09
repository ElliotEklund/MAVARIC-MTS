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
    set(EQUILIBRIUM_DIR "${ROOT}/sims/csrpmd/equilibrium")
    set(SAMPLING_DIR "${ROOT}/sims/csrpmd/sampling")
    set(DYNAMICS_DIR "${ROOT}/sims/csrpmd/dynamics")
    set(HAMILTONIAN_DIR "${ROOT}/src/hamiltonian")
    set(ESTIMATORS_DIR "${ROOT}/src/estimators")
    set(FORCES_DIR "${ROOT}/src/forces")
    set(INTEGRATORS_DIR "${ROOT}/src/integrators")
    set(MISC_DIR "${ROOT}/src/misc")
    set(EXECUTABLE_OUTPUT_PATH "${ROOT}/bin")
    set(CMAKE_CXX_FLAGS "-DMAC")

add_definitions( "-DNDEBUG -O3" )

# add the executable

file(GLOB   HAMILTONIAN_SRC ${HAMILTONIAN_DIR}/*.cpp)
file(GLOB   FORCES_SRC ${FORCES_DIR}/*.cpp)
file(GLOB   INTEGRATORS_SRC ${INTEGRATORS_DIR}/*.cpp)
file(GLOB   MISC_SRC ${MISC_DIR}/*.cpp)

add_executable(csrpmd csrpmd.cpp

	${HAMILTONIAN_SRC}
    ${FORCES_SRC}
    ${INTEGRATORS_SRC}
    ${MISC_SRC}

    #sampling directory
    ${SAMPLING_DIR}/csrpmd_sampling.cpp
)


include_directories(${BOOST_DIR})
include_directories(${EQUILIBRIUM_DIR}/Headers)
include_directories(${SAMPLING_DIR}/Headers)
include_directories(${DYNAMICS_DIR}/Headers)
include_directories(${HAMILTONIAN_DIR}/Headers)
include_directories(${ESTIMATORS_DIR}/Headers)
include_directories(${FORCES_DIR}/Headers)
include_directories(${INTEGRATORS_DIR}/Headers)
include_directories(${MISC_DIR}/Headers)

