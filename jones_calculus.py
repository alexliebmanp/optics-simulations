"""
    jones_calculus.py - this file contains classes for defining Jones matrices
                        and as well as tools for simulating the affect of
                        an optical systems on a polarized light beam.

                        all angles measured relative to y axis

                        based off of numpy.ndarray
    @author: oxide
"""
import numpy as np


## methods

def compute_optical_system(optical_components):
    """
    returns the net OpticalComponent object of an optical system.

    Arg:
        - optical_components: list of optical components in order in which they appear on from beam direction.
    """

    optical_system = np.identity(2)
    optical_components.reverse()
    for oc in optical_components:
        optical_system = optical_system @ oc.jones

    return OpticalComponent(optical_system)


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
            - angle_rad: angle of rotation in radians

        '''

        self.jones = self.get_rotated(angle_rad)
        self.angle = self.angle + angle_rad

    def set_angle(self, angle_rad):
        '''
        set absolute angle of OpticalComponent to angle_rad

        Args:
            - angle_rad: angle of rotation in radians
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

    def __init__(self, angle_rad):

        self.jones = np.asarray([[0,0],[0,1]])
        self.angle = 0
        self.set_angle(angle_rad)

class Retarder(OpticalComponent):
    """
    a general retarder such as a birefringent materials
    Args:
        - phix: phase shift along the x axis
        - phiy: phase shift along y axis

    Todo: how to relate dielectric tensor directly to Jones matrix.
    """

    def __init__(self, angle_rad, phix, phiy):

        self.jones = np.asarray([[np.exp(1j*phix), 0], [0, np.exp(1j*phiy)]])
        self.angle = 0
        self.set_angle(angle_rad)

class QWP(Retarder):
    """
    a quarter wave plate (valid only for one wavelength), with fast axis vertical
    """

    def __init__(self, angle_rad):

        self.jones = np.asarray([[1,0],[0,-1j]])
        self.angle = 0
        self.set_angle(angle_rad)

class HWP(Retarder):
    """
    half wave plate, fast axis vertical
    """

    def __init__(self, angle_rad):

        self.jones = np.asarray([[-1,0],[0,1]])
        self.angle = 0
        self.set_angle(angle_rad)
