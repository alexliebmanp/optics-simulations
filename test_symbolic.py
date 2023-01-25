from optics_simulations.jones_calculus import *
import numpy as np
import matplotlib.pyplot as plt
sp.init_printing()

# <codecell>
# setup system and initial condition

r11, r12, r21, r22, phi1, phi2 = sp.symbols('r_{11}, r_{12}, r_{21]}, r_{22}, phi_1, phi_2')

material_jones = sp.matrices.Matrix([[r11, r12], [r12, -r11]])
material = OpticalComponent(material_jones, symbolic_flag=True)
hwp1 = HWP(0, symbolic_flag=True)
hwp2 = HWP(0, symbolic_flag=True)
system = [hwp1, material, hwp2]
polarization_0 = sp.matrices.Matrix([0,1])
material.jones
hwp1.jones


# <codecell>
# single rotation
hwp1.set_angle(0)
hwp2.set_angle(phi1)
single_rotate = compute_optical_system(system)
single_rotate.jones
polarization_f = single_rotate.jones*polarization_0
polarization_f
signal = sp.trigsimp(polarization_f[1]**2 - polarization_f[0]**2)
signal

# <codecell>
# symmetric rotation
hwp1.set_angle(phi1)
hwp2.set_angle(phi1)
symmetric_rotate = compute_optical_system(system)
symmetric_rotate.jones
polarization_f = symmetric_rotate.jones*polarization_0
polarization_f
signal = sp.trigsimp(sp.simplify(polarization_f[1]**2 - polarization_f[0]**2))
signal

# <codecell>
# balanced rotation
hwp1.set_angle(phi1)
hwp2.set_angle(phi1 + sp.pi/8)
hwp2.jones
balanced_rotate = compute_optical_system(system)
balanced_rotate.jones
polarization_f = balanced_rotate.jones*polarization_0
polarization_f
signal = sp.trigsimp(polarization_f[1]**2 - polarization_f[0]**2)
signal
