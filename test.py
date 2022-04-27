from jones_calculus import *
import numpy as np
import matplotlib.pyplot as plt

Ei = np.asarray([0,1])
crystal = Retarder(np.pi/16, 0.5, 0.1)

hwp1 = HWP(0)
hwp2 = HWP(0)

def corotate(angle, optic_list, ind1, ind2):
    optic_list[ind1].rotate(angle)
    optic_list[ind2].rotate(angle)
    ol = optic_list
    system = compute_optical_system(ol)
    optic_list[ind1].rotate(-angle)
    optic_list[ind2].rotate(-angle)
    return system

angles = np.linspace(0,np.pi,1000)
fig, [ax1, ax2] = plt.subplots(1,2)

## tests

## reflectivity

S_reflectivity = np.zeros(len(angles))
for ii, a in enumerate(angles):
    osi = corotate(a, [hwp1, crystal, hwp2], 0, 2)
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_reflectivity[ii] = S

ax1.plot(angles, S_reflectivity)
ax1.set_ylim(0,1.5)

## birefringence
hwp2.rotate(np.pi/8)

S_bf = np.zeros(len(angles))
for ii, a in enumerate(angles):
    osi = corotate(a, [hwp1, crystal, hwp2], 0, 2)
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_bf[ii] = S

ax2.plot(angles, S_bf)
#ax2.set_ylim(0,1.5)

plt.show()
