import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import linecache

dir_path = os.path.dirname(os.path.realpath(__file__))
root = os.path.abspath(os.path.join(dir_path, os.pardir))
traj_path = root + "/Results/Trajectories/Q"
bead_path = root + "/InputFiles/SystemParameters.txt"

bead_content = linecache.getline(bead_path,1)
test = bead_content.split(":",1)[1].rstrip("\n")

#number of processors passed in
num_procs = int(sys.argv[1])
nuc_beads = int(1)

Q_total = np.zeros(1)

for i in range(num_procs):
    file_name = traj_path + str(i)
    X = np.loadtxt(file_name)
    X_temp = np.mean(X.reshape(-1,nuc_beads),axis=1)
    Q_total = np.concatenate((Q_total,X_temp))

hist, bin_edges = np.histogram(Q_total,density=True,bins=300)

num_bins = np.size(bin_edges)
bin_centers = np.zeros(num_bins-1)

for i in range(num_bins-1):
    bin_centers[i] = (bin_edges[i] + bin_edges[i+1])/2.0


to_file = np.zeros((np.size(bin_centers),2))
to_file[:,0] = bin_centers
to_file[:,1] = hist
np.savetxt("hist.txt",to_file,fmt='%8e')


plt.plot(bin_centers,hist)
plt.show()
