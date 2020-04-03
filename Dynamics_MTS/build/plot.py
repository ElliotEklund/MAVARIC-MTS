import numpy as np
import matplotlib.pyplot as plt
import sys

X_mac = np.transpose(np.loadtxt("traj_file_MAC"))
X_linux = np.transpose(np.loadtxt("traj_file_LINUX"))
#X_astra = np.transpose(np.loadtxt("traj_file_ASTRA"))

num = np.size(X_mac[:,0])

t = np.linspace(0,50,num)

traj = int(sys.argv[1])

plt.subplot(2,1,1)
plt.plot(t,X_mac[:,traj],label='mac')
plt.plot(t,X_linux[:,traj],label='linux')
#plt.plot(t,X_astra[:,traj],label='astra')

plt.subplot(2,1,2)
plt.plot(t,X_mac[:,traj]-X_linux[:,traj])
plt.xlabel('Time a.u.')
plt.ylabel('Cqq(t)')

plt.legend()
plt.show()
