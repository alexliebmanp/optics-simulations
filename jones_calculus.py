"""
    jones_calculus.py - this file contains classes for defining Jones matrices
                        and as well as tools for simulating the affect of
                        an optical systems on a polarized light beam.

                        all angles measured relative to y axis ('vertical direction')

                        based off of numpy.ndarray

                        Todo:
                        (1) come up with better way to act an optical system on a polarization vector, as opposed to doing it by calling the Jones matrix and using @ functional.
                        (2) incorporate formalism for getting R and T from dielectric tensor directly.
                        (3) separate into two files, one for class definitions and one for various functions. alternatively, can make an optical system a class on it's own with its own methods such as compute_optical_system() and corotate()
                        (4) Handle 3D polarization vectors

    @author: oxide
"""
import numpy as np


## methods

def compute_optical_system(optical_components):
    """
    returns the net OpticalComponent object of an optical system.

    Arg:
        - optical_components (list): list of optical components in order in which they appear on from beam direction.
    """

    optical_system = np.identity(2)
    optical_components.reverse()
    for oc in optical_components:
        optical_system = optical_system @ oc.jones

    return OpticalComponent(optical_system)

def corotate(angle_rad, optical_system, ind_list):
    '''
        corotates a set of optical components in a system and returns the rotated systems. Original optical system is left unmodified.

    Args:
        - angle_rad (float): angle to rotate in radians
        - optical_system (list): optical system (list of components) by which to rotate
        - ind_list (list): list of indices that specify which components to rotate
    '''
    for i in ind_list:
        optical_system[i].rotate(angle_rad)
    rotated_system = compute_optical_system(optical_system)
    for i in ind_list:
        optical_system[i].rotate(-angle_rad)
    return rotated_system

## optical components

class OpticalComponent:

    """
    parent class for Jones matrices, which is an identity matrix along with universal functions
    """

    def __init__(self, jones_matrix=np.identity(2)):

        self.jones = jones_matrix
        self.angle = 0

    ## Methods

    def rotate(self, angle_rad):
        '''
        rotates OpticalComponent by angle_rad relative to current angle.

        Args:
            - angle_rad (float): angle of rotation in radians

        '''

        self.jones = self.get_rotated(angle_rad)
        self.angle = self.angle + angle_rad

    def set_angle(self, angle_rad):
        '''
        set absolute angle of OpticalComponent to angle_rad

        Args:
            - angle_rad (float): angle of rotation in radians
        '''

        current_angle = self.angle
        rotation_angle = angle_rad - current_angle
        self.rotate(rotation_angle)

    def get_rotated(self, angle_rad):
        rotation_matrix = np.asarray([[np.cos(angle_rad), -np.sin(angle_rad)],[np.sin(angle_rad), np.cos(angle_rad)]])
        return rotation_matrix.transpose() @ self.jones @ rotation_matrix

    def get_jones(self):
        return self.jones

    def get_angle(self):
        return self.angle

class Polarizer(OpticalComponent):

    """
    a linear polarizer which projects a Jones vector onto an axis defined by the angle.
    """

    def __init__(self, angle_rad=0):

        self.jones = np.asarray([[0,0],[0,1]])
        self.angle = 0
        self.set_angle(angle_rad)

class Retarder(OpticalComponent):
    """
    a general retarder such as a birefringent material
    Args:
        - phix (float): phase shift along the x axis
        - phiy (float): phase shift along y axis

    Todo: how to relate dielectric tensor directly to Jones matrix.
    """

    def __init__(self, phix, phiy, angle_rad=0):

        self.jones = np.asarray([[np.exp(1j*phix), 0], [0, np.exp(1j*phiy)]])
        self.angle = 0
        self.set_angle(angle_rad)

class QWP(Retarder):
    """
    a quarter wave plate (valid only for one wavelength), with fast axis vertical
    """

    def __init__(self, angle_rad=0):

        self.jones = np.asarray([[1,0],[0,-1j]])
        self.angle = 0
        self.set_angle(angle_rad)

class HWP(Retarder):
    """
    half wave plate, fast axis vertical
    """

    def __init__(self, angle_rad=0):

        self.jones = np.asarray([[-1,0],[0,1]])
        self.angle = 0
        self.set_angle(angle_rad)

class GeneralSample(OpticalComponent):
    '''
    a general sample with transmission/reflection matrix of the form:

            |    r0 + dr    gamma - delta |
      T/R = |                             |
            | gamma + delta    r0 - dr    |

    note that all components can be real or imaginary

    ** Is there a better name for this class??**

    Args:
        - r0 (float):       dominant diagonal component, representing reflectivity
        - dr (float):       difference in diagonal, representing birefringence
        - gamma (float):    symmetric off-diagonal
        - delta (float):    antisymmetric off-diagonal
    '''

    def __init__(self, r0, dr, gamma, delta, angle_rad=0):

        self.jones = np.asarray([[r0 + dr, gamma - delta],[gamma + delta, r0 - dr]])
        self.angle = 0
        self.set_angle(angle_rad)
