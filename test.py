from jones_calculus import *
import numpy as np
import matplotlib.pyplot as plt

Ei = np.asarray([0,1])
r0 = 0.5j
dr = 0.1j
gamma = 0
delta = 0
phix = 0.5
phiy = 0.1
r0 = (1/2)*(np.exp(1j*phix) + (1/2)*np.exp(1j*phiy))
dr = (1/2)*(np.exp(1j*phix) - (1/2)*np.exp(1j*phiy))
crystal = GeneralSample(0, r0, dr, gamma, delta)

hwp1 = HWP(0)
hwp2 = HWP(0)

angles = np.linspace(0,np.pi/2,1000)
fig, [ax1, ax2] = plt.subplots(1,2)

## tests

## reflectivity

S_reflectivity = np.zeros(len(angles))
for ii, a in enumerate(angles):
    osi = corotate(a, [hwp1, crystal, hwp2], [0, 2])
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_reflectivity[ii] = S

ax1.plot(angles*2, S_reflectivity)
ax1.set_ylim(0,2)

## birefringence
hwp2.rotate(np.pi/8)

S_bf = np.zeros(len(angles))
for ii, a in enumerate(angles):
    osi = corotate(a, [hwp1, crystal, hwp2], [0, 2])
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_bf[ii] = S

ax2.plot(angles*2, S_bf)
#ax2.set_ylim(0,1.5)

plt.show()
