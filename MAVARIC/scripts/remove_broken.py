import numpy as np
import sys

def clean_me(stride_len, file_root, file_tag, broken):

    file_name = file_root + file_tag
    X = np.loadtxt(file_name)
    num_traj = int(np.size(X)/stride_len)

    myFile = open(file_root + file_tag + "_clean","w")
    counter = 0

    for i in range(num_traj):

        if counter < num_broken:
            if i == broken[counter]:
                counter += 1

            else:
                stride = i * stride_len
                for j in range(stride_len):
                    myFile.write(str(X[stride + j]) + "\n")

        else:
            stride = i * stride_len
            for j in range(stride_len):
                myFile.write(str(X[stride + j]) + "\n")


#number of processors passed in
num_procs = int(sys.argv[1])

nuc_beads = 4
elec_beads = 4
num_states = 2

file_root = "broken"
X_final = []

for i in range(num_procs):
    file_name = file_root + str(i)
    X = np.loadtxt(file_name)

    if np.size(X) > 1:
        for j in range(1,np.size(X)):
            X_final.append(int(X[j]))

X_final.sort()
num_broken = len(X_final)

file_root = "/Users/ellioteklund/Desktop/Dynamics_MTS_git/Dynamics_MTS/Results/Trajectories/"

clean_me(nuc_beads,file_root,"Q",X_final)
clean_me(nuc_beads,file_root,"P",X_final)
clean_me(elec_beads*num_states,file_root,"xelec",X_final)
clean_me(elec_beads*num_states,file_root,"pelec",X_final)
