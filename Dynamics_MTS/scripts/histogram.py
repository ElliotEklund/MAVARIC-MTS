import numpy as np
import matplotlib.pyplot as plt
import sys

root = "/Users/ellioteklund/Desktop/Dynamics_MTS_git/Dynamics_MTS/Results/Trajectories/Q"

#number of processors passed in
num_procs = int(sys.argv[1])
nuc_beads = int(1)

Q_total = np.zeros(1)
test = np.linspace(1,10,10)

test2 = np.concatenate((test,Q_total))

for i in range(num_procs):
    file_name = root + str(i)
    X = np.loadtxt(file_name)
    X_temp = np.mean(X.reshape(-1,nuc_beads),axis=1)
    Q_total = np.concatenate((Q_total,X_temp))

hist, bin_edges = np.histogram(Q_total,density=True,bins=100)

num_bins = np.size(bin_edges)
bin_centers = np.zeros(num_bins-1)

for i in range(num_bins-1):
    bin_centers[i] = (bin_edges[i] + bin_edges[i+1])/2.0

plt.plot(bin_centers,hist)
plt.show()
