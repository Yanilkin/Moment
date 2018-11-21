import numpy as np
import parameters as p

#I. Temperature independent coefficients
a = 3.14e-10 # lattice constant in m
b = np.sqrt(3)/2*a # Burger's vector sqrt(3)/2*a
Va = (a*a*a)/2.0 # Atomic volume
pi = 3.14 # Pui number

#k2s = 1.2e13*0 # Surface sink
k2s = 5e13*0 # Surface sink
k2d = 5e18*0 # Dislocation network sink
zid = 3.0 # Bias factor
zvd = 1.0 # Bias factor
k2is = k2s
k2vs = k2s
k2id = zid*k2d
k2vd = zvd*k2d

k2iv = 4*3.14*(3.1e-10+2.9e-10)/Va # Recombination rate constant


# sink strength of vacancies and pores
max_n = p.max_n
dn1 = p.dn1 
base = p.base
dx = p.dx
# Clusters initialization
sizes = np.arange(1,dn1)
for i in range(dn1,max_n+1):
	#sizes=np.hstack((sizes,(2-base)*i+(base-1)*dn1-1+np.power(base,i-dn1)))
	#sizes=np.hstack((sizes,dn1+(i-dn1)*dx))
	sizes=np.hstack((sizes,sizes[i-2]+1+int(sizes[i-2]*(base-1))))

sizes_av = np.zeros(max_n)+1.0
sizes_av[1:] = (sizes[:-1]+1 + sizes[1:])/2.0
dsizes = np.zeros(max_n)+1
dsizes[1:] = sizes[1:]-sizes[:-1]
sigma = np.zeros(max_n)+1.0
sigma[1:] = dsizes[1:]*dsizes[1:]/12.0
#sigma[-1] = (25+36)/2.0-(5.5)**2
print sizes
#print sigma
#print dsizes
#print dsizes/sizes
#print sizes_av

k2vp = [k2iv]
for i in range(1,max_n):
	#k2vp.append(4*pi/Va*(defpr.set_vac_rad(sizes[i])+defpr.set_vac_rad(1))) # Last one is the correction of rate constant
	k2vp.append(k2iv*sizes[i]**0.33333)

#k2vp = np.zeros(max_n)+1.0
#k2vp[0] = 1e-6

k2ip = [k2iv]
for i in range(1,max_n):
	#k2ip.append(4*pi/Va*(defpr.set_vac_rad(sizes[i])+defpr.set_vac_rad(1))) # Last one is the correction of rate constant
	k2ip.append(k2iv*sizes[i]**0.33333)

#k2ip = np.zeros(max_n)+1.0

k2vp_05 = [k2iv]
for i in range(1,max_n):
	#k2vp.append(4*pi/Va*(defpr.set_vac_rad(sizes[i])+defpr.set_vac_rad(1))) # Last one is the correction of rate constant
	k2vp_05.append(k2iv*sizes_av[i]**0.33333)

#k2vp_05 = np.zeros(max_n)+1.0
#k2vp_05[0] = 1e-6

k2ip_05 = [k2iv]
for i in range(1,max_n):
	#k2ip.append(4*pi/Va*(defpr.set_vac_rad(sizes[i])+defpr.set_vac_rad(1))) # Last one is the correction of rate constant
	k2ip_05.append(k2iv*sizes_av[i]**0.33333)

#k2ip_05 = np.zeros(max_n)+1.0

#II. Temperature dependent coefficients
T = p.T
# Diffusion coefficients
Di = 1e-8*np.exp(-0.4*11600/T)
#Dv = 1.5e-6*np.exp(-1.3*11600/T)
Dv = 3.7e-6*np.exp(-1.54*11600/T) #From Ivan Novoselov

#Di = 1.0
#Dv = 1.0


# Equlibrium concentration for dissociation rate
#from Novoselov
Ef = 3.0 #eV
a = 5.4387 # From MD Ivan Novoselov
b = 2.9587 # From MD Ivan Novoselov
cveq = np.exp(-(Ef+a*np.power(sizes-1,0.6667)-a*np.power(sizes,0.6667))*11600/T)
cveq[0] = np.exp(-Ef*11600/T)
cveq[1] = np.exp(-0.2*11600/T) #Binding energy for divacancy
cveq = 1e-10*cveq

#Production bias
prbias = np.zeros(max_n)
#prbias = np.hstack((np.array([1.55,0.227,0.0688,0.0414,0.0236,0.0119,0.005]),prbias))
#prbias = np.hstack((np.array([1.0,0.0,0.0,0.0,0.0,0.0,0.0]),prbias))
prbias[0] = 1.0
prbias = prbias/np.sum(prbias*sizes)