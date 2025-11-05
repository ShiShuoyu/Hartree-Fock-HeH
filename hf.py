import re
import numpy as np
from numpy.linalg import norm
import mpire
from numpy.linalg import norm
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import tplquad
from scipy.integrate import quad

#######################################
## get sto3g orbit of H and He         
#######################################
def gaussian(alpha):
    a = (2*alpha/np.pi)**0.75
    return lambda x: a*np.exp(-alpha * x**2)
He_data = np.array([
     [0.6362421394E+01,       0.1543289673E+00],
     [0.1158922999E+01,       0.5353281423E+00],
     [0.3136497915E+00,       0.4446345422E+00],
])
H_data = np.array([
      [0.3425250914E+01,       0.1543289673E+00],
      [0.6239137298E+00,      0.5353281423E+00],
      [0.1688554040E+00,       0.4446345422E+00],
])
def He_orbit(r):
    result = 0
    for i in He_data:
        result += i[1]*gaussian(i[0])(r)
    return result
def H_orbit(r):
    result = 0
    for i in H_data:
        result += i[1]*gaussian(i[0])(r)
    return result

def He_wavefunction(x,y,z):
    return He_orbit(norm(np.array([x,y,z])))

def H_wavefunction(x,y,z):
    return H_orbit(norm(np.array([x,y,z-1.4632])))


########################################
## varify normalization                 
########################################
# print(quad(lambda x: 4*np.pi* x**2 * H_orbit(x)**2, 0, np.inf)[0])
# print(quad(lambda x: 4*np.pi* x**2 * He_orbit(x)**2, 0, np.inf)[0])



'''不需要这一段'''
####################################################
## get the linear combination of the above two orbit
####################################################
# 原点设为He，H位于(0,0,1.4632)
# def wave_function(He_const, H_const):
#     def wave_fun(x):
#         H_position = np.array([0, 0, 1.4632])
#         r_He = norm(x)
#         r_H = norm(x-H_position)
#         return He_const*He_orbit(r_He) + H_const*H_orbit(r_H)
#     return wave_fun



######################
## calc overlap matrix
######################
def inner_product(phi1, phi2):
    return tplquad(lambda x,y,z: phi1(x,y,z)*phi2(x,y,z), -np.inf, np.inf, -np.inf, np.inf, -np.inf, np.inf)
    # return tplquad(lambda x,y,z: phi1(x,y,z)*phi2(x,y,z), -10, 10, -10, 10, -10, 10)
    # return tplquad(lambda x,y,z: phi1(x,y,z)*phi2(x,y,z), -1, 1, -1, 1, -1, 1)

# print(inner_product(He_wavefunction, He_wavefunction)) # = 1

''' here to out put'''
# print(inner_product(He_wavefunction, H_wavefunction)[0]) # 计算overlap

################################
## calc single particle H matrix
################################
"""some code"""
H = np.array([[0,0],[0,0]])


##########################################################################################
## a function to generate hatree matrix, depends on current a and b (a for He and b for H)
##########################################################################################
def genHatree(a, b):
    '''some code here'''
    return np.array([0,0],[0,0])
