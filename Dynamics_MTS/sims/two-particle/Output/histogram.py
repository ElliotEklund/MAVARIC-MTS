import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import linecache

#dir_path = os.path.dirname(os.path.realpath(__file__))
#root = os.path.abspath(os.path.join(dir_path, os.pardir))
#traj_path = root + "/Results/Trajectories/Q"
#bead_path = root + "/InputFiles/SystemParameters.txt"
#
#bead_content = linecache.getline(bead_path,1)
#test = bead_content.split(":",1)[1].rstrip("\n")
#
##number of processors passed in
#num_procs = int(sys.argv[1])
#nuc_beads = int(1)

file_name1 = "P1"
file_name2 = "P2"


#Q_total = np.zeros(1)

#for i in range(num_procs):
#    file_name = traj_path + str(i)
#    X = np.loadtxt(file_name)
#    X_temp = np.mean(X.reshape(-1,nuc_beads),axis=1)
#    Q_total = np.concatenate((Q_total,X_temp))

Q_total1 = np.loadtxt(file_name1)
Q_total2 = np.loadtxt(file_name2)

hist1, bin_edges1 = np.histogram(Q_total1,density=True,bins=300)
hist2, bin_edges2 = np.histogram(Q_total2,density=True,bins=300)

num_bins1 = np.size(bin_edges1)
num_bins2 = np.size(bin_edges2)

bin_centers1 = np.zeros(num_bins1-1)
bin_centers2 = np.zeros(num_bins2-1)

for i in range(num_bins1-1):
    bin_centers1[i] = (bin_edges1[i] + bin_edges1[i+1])/2.0
for i in range(num_bins2-1):
    bin_centers2[i] = (bin_edges2[i] + bin_edges2[i+1])/2.0


to_file1 = np.zeros((np.size(bin_centers1),2))
to_file2 = np.zeros((np.size(bin_centers2),2))

to_file1[:,0] = bin_centers1
to_file2[:,0] = bin_centers2

to_file1[:,1] = hist1
to_file2[:,1] = hist2

np.savetxt("hist1.txt",to_file1,fmt='%8e')
np.savetxt("hist2.txt",to_file2,fmt='%8e')

plt.plot(bin_centers1,hist1)
plt.plot(bin_centers2,hist2)

plt.show()
