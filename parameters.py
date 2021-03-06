#Parameters
G = 1e-4
max_dose = 1e2
T = 900

# size discretization
max_n = 200
dn1 = 11
base = 1.10
dx = 5

# integration and output
integ_steps = 10000
nframes = 100


# initial distribution
import materials_constants as mc
import numpy as np
L_in = np.zeros(2*max_n)
L_in[0] = np.sqrt(G/mc.k2iv/mc.Di)
#L_in[0] = 1.5e-8
#L_in[0] = 1e-6
scale = 1.5e-8
#L_in[450:470] = scale*np.exp(-np.power((mc.sizes[450:470]-mc.sizes[460]+1.0)/(mc.sizes[460]-mc.sizes[455]),2))

#L_in[0] = 1e0
#L_in[21] = 1e-3
#co[248] = 1e-4
#co[133] = 1e-4

k2i = mc.k2is + mc.k2id + np.sum(mc.k2ip[1:-1]*L_in[1:max_n-1])
k2v = mc.k2vs + mc.k2vd + np.sum(mc.k2vp[1:-1]*L_in[1:max_n-1])
#co[0] = np.sqrt(np.power(k2i*mc.Di/(2*mc.k2iv*(mc.Di+mc.Dv)),2)+G*k2i*mc.Di/(mc.k2iv*(mc.Di+mc.Dv)*k2v*mc.Dv))-k2i*mc.Di/(2*mc.k2iv*(mc.Di+mc.Dv))
#print mc.Dv*co[0]
#print mc.Di*(np.sqrt(np.power(k2v*mc.Dv/(2*mc.k2iv*(mc.Di+mc.Dv)),2)+G*k2v*mc.Dv/(mc.k2iv*(mc.Di+mc.Dv)*k2i*mc.Di))-k2v*mc.Dv/(2*mc.k2iv*(mc.Di+mc.Dv)))

