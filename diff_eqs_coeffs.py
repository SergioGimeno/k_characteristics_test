import numpy as np
import F1ast
import F2ast
import Gast
import drlogut
import dthetalogut
import drF2ast
import dthetaGast
import dthetaF1ast


def cart_to_polar(coords):
	return [np.sqrt(coords[0]**2+coords[1]**2), np.arctan2(coords[0],coords[1])]

def polar_to_cart(coords):
	return [coords[0]*np.sin(coords[1]), coords[0]*np.cos(coords[1])]

def F0(epsilon, kappa, gamma):
	return (1+kappa*epsilon)/(-kappa*gamma)

def depsilonF0(gamma):
	return 1./(-gamma)

def rterm_k(M, a, r, theta, l, s0, kappa, gamma):
	return F2ast.f2ast_f(M, a, r, theta, l, s0, kappa, gamma)

def thetaterm_k(M, a, r, theta, l, s0, kappa, gamma):
	return r*(-F1ast.f1ast_f(M, a, r, theta, l, s0, kappa, gamma) -Gast.gast_f(M, a, r, theta, l, s0, kappa, gamma))

def indeptermpolar_k(M, a, r, theta, l, s0, kappa, gamma, k): 
	return -(drF2ast.drf2ast_f(M, a, r, theta, l, s0, kappa, gamma) - dthetaGast.dthetagast_f(M, a, r, theta, l, s0, kappa, gamma) - dthetaF1ast.dthetaf1ast_f(M, a, r, theta, l, s0, kappa, gamma) + depsilonF0(gamma)*(-F2ast.f2ast_f(M, a, r, theta, l, s0, kappa, gamma)*drlogut.drlogut_f(M, a, r, theta, l) + F1ast.f1ast_f(M, a, r, theta, l, s0, kappa, gamma)*dthetalogut.dthetalogut_f(M, a, r, theta, l) + Gast.gast_f(M, a, r, theta, l, s0, kappa, gamma)*dthetalogut.dthetalogut_f(M, a, r, theta, l)))*k

def rterm_epsilon(theta):
	return np.cos(theta)

def thetaterm_epsilon(theta):
	return -np.sin(theta)

def indeptermpolar_epsilon(M, a, r, theta, l, s0, kappa, gamma, k, epsilon):
	return F0(epsilon, kappa, gamma)*(rterm_epsilon(theta) * drlogut.drlogut_f(M, a, r, theta, l) + (thetaterm_epsilon(theta)/r)*dthetalogut.dthetalogut_f(M, a, r, theta, l)) + k*(-rterm_epsilon(theta) * F1ast.f1ast_f(M, a, r, theta, l, s0, kappa, gamma) - rterm_epsilon(theta) * Gast.gast_f(M, a, r, theta, l, s0, kappa, gamma) - (thetaterm_epsilon(theta)/r)*F2ast.f2ast_f(M, a, r, theta, l, s0, kappa, gamma))

################### cartesian terms

def xterm_k(M, a, x, y, l, s0, kappa, gamma):
	polar_coords = cart_to_polar([x, y])
	return rterm_k(M, a, polar_coords[0], polar_coords[1], l, s0, kappa, gamma)*np.sin(polar_coords[1]) + thetaterm_k(M, a, polar_coords[0], polar_coords[1], l, s0, kappa, gamma)*np.cos(polar_coords[1])

def yterm_k(M, a, x, y, l, s0, kappa, gamma):
	polar_coords = cart_to_polar([x, y])
	return rterm_k(M, a, polar_coords[0], polar_coords[1], l, s0, kappa, gamma)*np.cos(polar_coords[1]) - thetaterm_k(M, a, polar_coords[0], polar_coords[1], l, s0, kappa, gamma)*np.sin(polar_coords[1])

def indeptermcart_k(M, a, x, y, l, s0, kappa, gamma, k):
	polar_coords = cart_to_polar([x, y])
	return indeptermpolar_k(M, a, polar_coords[0], polar_coords[1], l, s0, kappa, gamma, k)

def xterm_epsilon(x, y):
	# polar_coords = cart_to_polar([x, y])
	# return rterm_epsilon(polar_coords[1])*np.sin(polar_coords[1]) + thetaterm_epsilon(polar_coords[1])*np.cos(polar_coords[1])
	return 0.

def yterm_epsilon(x, y):
	# polar_coords = cart_to_polar([x, y])
	# return rterm_epsilon(polar_coords[1])*np.cos(polar_coords[1]) - thetaterm_epsilon(polar_coords[1])*np.sin(polar_coords[1])
	return 1.

def indeptermcart_epsilon(M, a, x, y, l, s0, kappa, gamma, k, epsilon):
	polar_coords = cart_to_polar([x, y])
	return indeptermpolar_epsilon(M, a, polar_coords[0], polar_coords[1], l, s0, kappa, gamma, k, epsilon)

################### 'normalized' cartesian terms

def Ak(M, a, x, y, l, s0, kappa, gamma):
	return xterm_k(M, a, x, y, l, s0, kappa, gamma)/yterm_k(M, a, x, y, l, s0, kappa, gamma)

def Bk():
	return 1.

def Ck(M, a, x, y, l, s0, kappa, gamma, k):
	return indeptermcart_k(M, a, x, y, l, s0, kappa, gamma, k)/yterm_k(M, a, x, y, l, s0, kappa, gamma)


def Aepsilon(x, y):
	return 0.

def Bepsilon(x, y):
	return 1.

def Cepsilon(M, a, x, y, l, s0, kappa, gamma, k, epsilon):
	return indeptermcart_epsilon(M, a, x, y, l, s0, kappa, gamma, k, epsilon)/yterm_epsilon(x, y)

# print(polar_to_cart([4,np.pi/2]))
# print(polar_to_cart([4,np.pi/4]))

# print(cart_to_polar([4, 0]))
# print(cart_to_polar([2.82842712474619, 2.82842712474619]))


# print(F1ast.f1ast_f(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(F1ast.f1ast_f(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))

# print(F2ast.f2ast_f(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(F2ast.f2ast_f(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))

# print(Gast.gast_f(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(Gast.gast_f(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))

# print(drlogut.drlogut_f(1, 0.5, 10, np.pi/2,3.8))
# print(drlogut.drlogut_f(1, 0.5, 10, np.pi/3,3.8))

# print(dthetalogut.dthetalogut_f(1, 0.5, 10, np.pi/2,3.8))
# print(dthetalogut.dthetalogut_f(1, 0.5, 10, np.pi/3,3.8))

# print(drF2ast.drf2ast_f(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(drF2ast.drf2ast_f(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))

# print(dthetaGast.dthetagast_f(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(dthetaGast.dthetagast_f(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))

# print(dthetaF1ast.dthetaf1ast_f(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(dthetaF1ast.dthetaf1ast_f(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))