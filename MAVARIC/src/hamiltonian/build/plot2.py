import numpy as np
import matplotlib.pyplot as plt

def plot_diff(v1,v2,n):

    diff = v1 - v2

    for i in range(n):
        plt.plot(x1[:,0],diff[:,1+i])

    plt.show()

def plot_all(v1,v2,n):

    for i in range(n):
        plt.plot(x1[:,0],v1[:,1+i])
        plt.plot(x1[:,0],v2[:,1+i])

    plt.show()

nuc_beads = 6
elec_beads = 6
num_states = 2


Q1 = np.loadtxt("Q1")
Q2 = np.loadtxt("Q2")
P1 = np.loadtxt("P1")
P2 = np.loadtxt("P2")
x1 = np.loadtxt("xElec1")
x2 = np.loadtxt("xElec2")
p1 = np.loadtxt("pElec1")
p2 = np.loadtxt("pElec2")


#plot_diff(P1,P2,nuc_beads)
plot_diff(x1,x2,nuc_beads*num_states)
#plot_diff(p1,p2,nuc_beads*num_states)
#plot_diff(Q1,Q2,nuc_beads)
#plot_diff(P1,P2,nuc_beads)
