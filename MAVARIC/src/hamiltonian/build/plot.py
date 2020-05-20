import numpy as np
import matplotlib.pyplot as plt

nuc_beads = 6
elec_beads = 6
num_states = 2

Q_file1 = "Q_data_normal"
P_file1 = "P_data_normal"
x_file1 = "xElec_data_normal"
p_file1 = "pElec_data_normal"
Q1 = np.loadtxt(Q_file1)
P1 = np.loadtxt(P_file1)
x1 = np.loadtxt(x_file1)
p1 = np.loadtxt(p_file1)

Q_file2 = "Q_data_wrong"
P_file2 = "P_data_wrong"
x_file2 = "xElec_data_wrong"
p_file2 = "pElec_data_wrong"
Q2 = np.loadtxt(Q_file2)
P2 = np.loadtxt(P_file2)
x2 = np.loadtxt(x_file2)
p2 = np.loadtxt(p_file2)

Q_file3 = "Q_data_correct"
P_file3 = "P_data_correct"
x_file3 = "xElec_data_correct"
p_file3 = "pElec_data_correct"
Q3 = np.loadtxt(Q_file3)
P3 = np.loadtxt(P_file3)
x3 = np.loadtxt(x_file3)
p3 = np.loadtxt(p_file3)


dt1 = 0.005
dt2 = 0.001
time = 10
num_steps1 = time/dt2
num_steps2 = time/dt2
t1 = np.linspace(0,time,num_steps1)
t2 = np.linspace(0,time,num_steps2)

file_name = 'P'
fig = plt.figure()
ax1 = fig.add_subplot(4,1,1)
ax2 = fig.add_subplot(4,1,2)
ax3 = fig.add_subplot(4,1,3)
ax4 = fig.add_subplot(4,1,4)

ax1.set_ylabel(file_name)
ax1.set_title('time a.u.')
for i in range(nuc_beads):
    ax1.plot(t1,P1[:,i])

ax2.set_ylabel(file_name)
ax2.set_title('time a.u.')
for i in range(nuc_beads):
    ax2.plot(t2,P2[:,i])

ax3.set_ylabel(file_name)
ax3.set_title('time a.u.')
for i in range(nuc_beads):
    ax3.plot(t2,P3[:,i])

ax4.set_ylabel(file_name)
ax4.set_title('time a.u.')
for i in range(nuc_beads):
    ax4.plot(t2,P3[:,i]-P1[:,i])

plt.legend()
plt.show()
