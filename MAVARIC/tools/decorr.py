import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import linecache

def compute_decorr(k,mu,var,num_samples,x):

    x = x - mu
    x_shift = np.roll(x,-k)
    x_shift[num_samples-k-1:num_samples] = 0.0

    rk = np.dot(x,x_shift)/var

    return rk


dir_path = os.path.dirname(os.path.realpath(__file__))
root = os.path.abspath(os.path.join(dir_path, os.pardir))
traj_path = root + "/Results/decorrelation"
nuc_bead_path = root + "/InputFiles/SystemParameters.txt"
#elec_bead_path = root + "/InputFiles/SystemParameters.txt"

bead_content = linecache.getline(nuc_bead_path,2)
nuc_bead = bead_content.split(":",1)[1].rstrip("\n")
nuc_bead = int(nuc_bead)

#number of processors passed in
num_procs = int(sys.argv[1])

cent = np.loadtxt(traj_path)
num_samples = np.size(cent)
mu = np.average(cent)
var = np.var(cent)


num_k = 10
interval = 100

t = np.zeros(num_k)
y = np.zeros(num_k)
t[0] = 0

r0 = compute_decorr(0,mu,var,num_samples,cent)
y[0] = r0/r0

for i in range(num_k-1):
    r = compute_decorr((i+1)*interval,mu,var,num_samples,cent)
    t[i+1] = (i+1)*interval
    y[i+1] = r/r0

plt.plot(t,y)
plt.show()
