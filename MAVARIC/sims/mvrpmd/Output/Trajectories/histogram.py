import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import linecache

Q_total = np.loadtxt("Q")

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
