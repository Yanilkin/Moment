import numpy as np
from scipy.integrate import odeint

# Parameters
import parameters as p
G = p.G
max_n = p.max_n

# Material constants
import materials_constants as mc
Di = mc.Di
Dv =  mc.Dv
k2is = mc.k2is
k2id = mc.k2id
k2vs = mc.k2vs
k2vd = mc.k2vd
k2vp = mc.k2vp # Sink strength of pores
k2vp_05 = mc.k2vp_05 # Sink strength of pores
k2ip = mc.k2ip # Sink strength of pores
k2ip_05 = mc.k2ip_05 # Sink strength of pores
k2iv = mc.k2iv
cveq = mc.cveq
prbias = mc.prbias
sizes = mc.sizes
sizes_av = mc.sizes_av
dsizes = mc.dsizes
sigma = mc.sigma

def Dicif(k2p):
	return G/(k2is+k2id+k2p) # Calculation of SIA concentration from stationary

def f(L,t):
	"""Growth rate equations for pore distribution"""
	#c[0] = 1e-2
	L0 = L[0:max_n] # Zero moment function
	L1 = L[max_n:] # First moment function
	dL0 = np.zeros(max_n) # Diff zero moment function
	dL1 = np.zeros(max_n) # Diff first moment function
	J = np.zeros(max_n) # Flux J(xi)
	J_m05 = np.zeros(max_n) # Flux J(<x>i-1/2)

	c = L0 + L1*(sizes-sizes_av) # Concentation c(xi)
	c_m1 = L0 + L1*(sizes-dsizes+1-sizes_av) # Concentation c(x(i-1) + 1)
	c_m05 = L0 - L1*0.5 # Concentation c(<xi> - 1/2)
	c_p05 = L0 + L1*0.5 # Concentation c(<xi> + 1/2)
	
	#k2v = np.sum(k2vp[:-1]*c[:-1])+np.sum((dsizes[1:]-1.0)*k2vp_05[1:]*c_m05[1:])  # Total sink strength of voids for vacancies: first term from L0, second one from L1
	#k2i = np.sum(k2ip*c)+np.sum((dsizes[1:]-1.0)*k2vp_05[1:]*c_p05[1:])  # Total sink strength of voids for vacancies

	# I. Calculation of SIA concentration from stationary
	#Dici = Dicif(k2i)
	#Dici = 0
	Dvcv = Dv*c[0]
	Dici = 0.95*Dvcv

	# II. Changes in vacancy concentration
	#dL0[0] = G*prbias[0]*0 - ((k2vs+k2vd+k2v+k2vp[0]*c[0])*Dv + k2ip[0]*Dici)*c[0] + k2ip[1]*Dici*c[1] + 2*k2vp[0]*Dv*cveq[1]*c[1] + Dv*(np.sum(k2vp[1:-1]*c_m1[2:]*cveq[2:])+np.sum((dsizes[1:]-1.0)*k2vp_05[1:]*c_p05[1:]*cveq[1:]))
	dL0[0] = 0
	# III. Changes in clusters concentration
	# a. clusters
	J[0:max_n-1] = k2vp[0:max_n-1]*c[0:max_n-1]*Dvcv - k2ip[1:max_n]*c_m1[1:max_n]*Dici - k2vp[0:max_n-1]*c_m1[1:max_n]*cveq[1:max_n]*Dv # J(xi) = P(xi)f(xi) - Q(x(i+1))f(x(i+1))
	J[0] = J[0] + 0.0*k2vp[0]*c[0]*Dvcv
	J_m05[1:max_n] = k2vp_05[1:max_n]*c_m05[1:max_n]*Dvcv - k2ip_05[1:max_n]*c_p05[1:max_n]*Dici - k2vp_05[1:max_n]*c_p05[1:max_n]*cveq[1:max_n]*Dv # J(<x>i-1/2) = P(xi)f(<x>i-1/2) - Q(x(i))f(<x>i+1/2))	

	dL0[1:max_n-1] = 1.0/dsizes[1:max_n-1]*(G*prbias[1:max_n-1] + J[0:max_n-2] - J[1:max_n-1])
	dL1[1:max_n-1] = -(dsizes[1:max_n-1]-1)/(2.0*sigma[1:max_n-1])*(J[0:max_n-2] -2*J_m05[1:max_n-1] + J[1:max_n-1])

	# b. Last one
	dL0[max_n-1] = 1.0/dsizes[max_n-1]*J[max_n-2]
	dL1[max_n-1] = -(dsizes[max_n-1]-1)/(2.0*sigma[max_n-1]*dsizes[max_n-1])*(J[max_n-2] -2*J_m05[max_n-1])
	dL = np.hstack((dL0,dL1))
	#print dL0[0]
	#print dL1[max_n-1]*dsizes[max_n-1]*sigma[max_n-1]

	#print dL1
	#print np.sum((dL0*sizes_av+dL1*sigma)*dsizes)
	return dL

if __name__ == "__main__":

	#print "initil derivitives\n", f(p.co,0)

	# Time integration of grow equation
	t = np.linspace(0,p.max_dose/G,p.integ_steps)
	results = odeint(f,p.L_in,t)
	L0 = results[:,0:max_n] # Zero moment function
	L1 = results[:,max_n:] # First moment function

	print "final distribution L0\n", L0[-1]
	print "final distribution L1\n", L1[-1]
	
	#print "final derivitives\n", f(results[-1],0)
	print "total number of clusters\n", np.sum(L0[-1,1:]*dsizes[1:])

	print "total vacancy concentration\n", np.sum((L0[-1]*sizes_av+L1[-1]*sigma)*dsizes)
	print "total vacancy concentration\n", np.sum((L0[-1]*sizes_av)*dsizes)
	#print "total swelling\n," np.sum(results[-1]*sizes/Va)

	#save results
	fo = open('log','w')
	fo.write('0 ')
	for size in sizes:
		fo.write('%f ' % size)
	fo.write('\n')
	i = 0
	for rows in L0:
		if i%p.discretization ==0: 
			fo.write('%f '%(t[i]*G))
			for cols in rows:
				fo.write('%e ' % cols)
			fo.write('\n')
		i = i + 1
	fo.close()

