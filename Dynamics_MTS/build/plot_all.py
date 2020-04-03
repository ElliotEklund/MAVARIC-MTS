import numpy as np
import matplotlib.pyplot as plt

X_mac = np.loadtxt("pos_auto_corr_MAC")
X_linux = np.loadtxt("pos_auto_corr_LINUX")
#X_astra = np.loadtxt("pos_auto_corr_ASTRA")

num = np.size(X_mac[:,0])
t = np.linspace(0,50,num)

plt.subplot(2,1,1)
plt.plot(t,X_mac[:,1],label='mac')
plt.plot(t,X_linux[:,1],label='linux')
#plt.plot(t,X_astra[:,1],label='astra')

plt.subplot(2,1,2)
plt.plot(t,X_mac[:,1]-X_linux[:,1])

plt.xlabel('Time a.u.')
plt.ylabel('Cqq(t)')

plt.legend()
plt.show()
