import numpy as np
import matplotlib
matplotlib.use('qt4agg')
import matplotlib.pyplot as plt

data=np.loadtxt('PCP-PCN.txt')
plt.figure(1)
plt.plot(data[:,0], data[:,1])
plt.ylim([-10, 10])
plt.xlim([0, 1.2])
plt.show()
plt.get_current_fig_manager().window.raise_()
