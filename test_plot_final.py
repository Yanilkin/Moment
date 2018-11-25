import matplotlib.pyplot as plt
import numpy as np

#dirs = ['base1.05','base1.10','base1.20']
dirs = ['base1.05']
dirs = ['base1.10']

for di in dirs:
	dat = np.loadtxt('./%s/log' % di)
	sizes = dat[0,5:-1]
	dist = dat[-1,5:-1]
	vac_conc = dat[2,1]
	plt.plot(sizes[:-1]**0.333,dist[:-1]*sizes[:-1]**0.6666/vac_conc,marker='o',label='%s' % di)

#plt.axis([1,100,0,0.0000002])
plt.ylabel('Cluster concentration (C(i)/Cv)')
plt.xlabel('Cluster size (n^1/3)')
plt.legend(loc='best')
plt.savefig('final.png')
plt.show()