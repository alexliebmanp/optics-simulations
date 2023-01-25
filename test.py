from optics_simulations.jones_calculus import *
import numpy as np
import matplotlib.pyplot as plt

# hi

Ei = np.asarray([0,1])
r0 = 1
dr = 0.1
gamma = 0.1
delta = 0
phix = 0.5
phiy = 0.1
#r0 = (1/2)*(np.exp(1j*phix) + (1/2)*np.exp(1j*phiy))
#dr = (1/2)*(np.exp(1j*phix) - (1/2)*np.exp(1j*phiy))
'''
            |    r0 + dr    gamma - delta |
      T/R = |                             |
            | gamma + delta    r0 - dr    |
'''

#crystal = Retarder(phix, phiy)
#crystal = OpticalComponent()

crystal = GeneralSample(r0, dr, gamma, delta)
hwp1 = HWP(0)
hwp2 = HWP(0)

angles = np.linspace(0,np.pi/2,1000)
fig, [[ax1, ax2], [ax3,ax4]] = plt.subplots(2,2)

## tests

## reflectivity
hwp1.set_angle(0)
hwp2.set_angle(0)
S_reflectivity = np.zeros(len(angles))
for ii, a in enumerate(angles):
    osi = corotate(a, [hwp1, crystal, hwp2], [0, 2])
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_reflectivity[ii] = S

ax1.plot(angles*2, S_reflectivity)
#ax1.set_ylim(0,None)

## birefringence

hwp1.set_angle(0)
hwp2.set_angle(np.pi/8)
S_bf = np.zeros(len(angles))
for ii, a in enumerate(angles):
    osi = corotate(a, [hwp1, crystal, hwp2], [0, 2])
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_bf[ii] = S

ax2.plot(angles*2, S_bf)
#ax2.set_ylim(0,1.5)

## single axis 1

hwp1.set_angle(0)
hwp2.set_angle(0)
S_sa1 = np.zeros(len(angles))
for ii, a in enumerate(angles):
    hwp1.set_angle(a)
    osi = compute_optical_system([hwp1, crystal, hwp2])
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_sa1[ii] = S

ax3.plot(angles*2, S_sa1)

## single axis 2

hwp1.set_angle(0)
hwp2.set_angle(0)
S_sa2 = np.zeros(len(angles))
for ii, a in enumerate(angles):
    hwp2.set_angle(a)
    osi = compute_optical_system([hwp1, crystal, hwp2])
    Ef = osi.jones @ Ei
    S = abs(Ef[1])**2 - abs(Ef[0])**2
    S_sa2[ii] = S

ax4.plot(angles*2, S_sa2)

ax1.set_xlabel('Polarization Angle (radians)')
ax2.set_xlabel('Polarization Angle (radians)')
ax3.set_xlabel('Polarization Angle (radians)')
ax4.set_xlabel('Polarization Angle (radians)')
ax1.set_ylabel('Reflectivity (symmetric corotation)')
ax2.set_ylabel('Birefringence')
ax3.set_ylabel('Intensity (single axis 1)')
ax4.set_ylabel('Intensity (single axis 2)')

plt.tight_layout(w_pad=1)
plt.show()
