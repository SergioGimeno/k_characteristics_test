import math
import numpy as np
from scipy.optimize import bisect

def l_kepler(r, M):
	return (M**(1/2.))*(r**(2))/(r**(3/2.)-2*M*r**(1/2.))

def dr_W(r, theta, a, M, l):
	return  ((l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))*(((-16*a**2*M**2*r**3*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**3 + (8*a**2*M**2*r*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(2*r - (4*a**2*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a**2*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)) - ((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/(l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))) - ((l**2*((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2)) + (8*a*l*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 - (4*a*l*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(2*r - (4*a**2*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a**2*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))*((4*a**2*M**2*r**2*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))))/(l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))**2))/(2.*((4*a**2*M**2*r**2*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))))

def find_rmb():

	return 4.

def find_rms():

	return 6.

def find_W_border_eq(l, W_function, M, W_in):

	def find_W_eq(r):

		return W_function(l, r, np.pi/2, M) - W_in

	return find_W_eq

def find_r_in_and_out(find_W_eq, a, b):

	return bisection(find_W_eq, a, b, 10000)

def W_function(l, r, theta, M):

	return np.log(((-2 * M + r)*r**2*np.sin(theta)**2)/(-((l**2)*(-2 * M + r)) + (r**3)*np.sin(theta)**2))/2.

def W_eq(r, l, M):
	theta = np.pi/2
	return np.log(((-2 * M + r)*r**2*np.sin(theta)**2)/(-((l**2)*(-2 * M + r)) + (r**3)*np.sin(theta)**2))/2.

def W_eq_in_out(r, l, M, W_s):
	theta = np.pi/2
	return np.log(((-2 * M + r)*r**2*np.sin(theta)**2)/(-((l**2)*(-2 * M + r)) + (r**3)*np.sin(theta)**2))/2. - W_s

def cart_to_polar(x, y):
	r = np.sqrt(x**2 + y**2)
	theta = np.arctan2(x, y)
	return r, theta

def polar_to_cart(r, theta):
	x = r * np.sin(theta)
	y = r * np.cos(theta)
	return x, y

def W_function_cart(l, x, y, M):
	r, theta = cart_to_polar(x,y)
	return np.log(((-2 * M + r)*r**2*np.sin(theta)**2)/(-((l**2)*(-2 * M + r)) + (r**3)*np.sin(theta)**2))/2.

def obtain_r_c(dr_W, l, M):

	return bisect(dr_W, find_rms(), 200, args=(np.pi/2, 0, M, l))

def beta_m(beta_mc, rc, r, theta, gamma):
	return beta_mc*((rc*(rc - 2))/(r*(r - 2)*(np.sin(theta)**2)))**(gamma - 1)

def beta_m_cart(beta_mc, rc, x, y, gamma):
	r, theta = cart_to_polar(x,y)
	return beta_mc*((rc*(rc - 2))/(r*(r - 2)*(np.sin(theta)**2)))**(gamma - 1)

def p0(K, beta_m, W, Win, gamma):
	return K**(1/(1 - gamma))*((gamma*(1 + beta_m))/((1 - gamma)*(-Win + W)*beta_m))**(gamma/(1 - gamma))

def epsilon0(p0, K, gamma):
	return (p0/K)**(1/gamma)

def A_tilde(r, theta, l0):
	return -2*(-3 + r)*((-8*l0**2*(-3 + r)*(-2 + r) + 4*r**3*(21 - 14*r + 2*r**2))*np.cos(2*theta) + r**3*(-21 + 10*r)*np.cos(4*theta))

def A_tilde_cart(x, y, l0):
	r, theta = cart_to_polar(x,y)

	return -2*(-3 + r)*((-8*l0**2*(-3 + r)*(-2 + r) + 4*r**3*(21 - 14*r + 2*r**2))*np.cos(2*theta) + r**3*(-21 + 10*r)*np.cos(4*theta))

def B_tilde(r, theta, l0):
	return 4*(2*l0**2*(-2 + r)**2*r - (r**3*(-3 + 2*r) + 4*l0**2*(-2 + r)*(-9 + 4*r) - r**3*(-9 + 4*r)*np.cos(2*theta))*np.sin(theta)**2)

def B_tilde_cart(x, y, l0):
	r, theta = cart_to_polar(x,y)

	return 4*(2*l0**2*(-2 + r)**2*r - (r**3*(-3 + 2*r) + 4*l0**2*(-2 + r)*(-9 + 4*r) - r**3*(-9 + 4*r)*np.cos(2*theta))*np.sin(theta)**2)

def k1_func(r, theta, l0, beta_m, epsilon0, K, gamma):
	return -2*(r**3 + l0**2*(2 - r)*((1)/(np.sin(theta)))**2)*np.sin(theta)**2*beta_m*epsilon0 + K*gamma*(-2*r**4 - r**3*(-2 + beta_m) + 8*l0**2*(1 + beta_m) + 2*l0**2*r**2*(2 + beta_m) - 4*l0**2*r*(3 + beta_m) + r**3*np.cos(2*theta)*(-2 + 2*r + beta_m))*epsilon0**gamma

def k1_cart_func(x, y, l0, beta_m, epsilon0, K, gamma):
	r, theta = cart_to_polar(x,y)

	return -2*(r**3 + l0**2*(2 - r)*((1)/(np.sin(theta)))**2)*np.sin(theta)**2*beta_m*epsilon0 + K*gamma*(-2*r**4 - r**3*(-2 + beta_m) + 8*l0**2*(1 + beta_m) + 2*l0**2*r**2*(2 + beta_m) - 4*l0**2*r*(3 + beta_m) + r**3*np.cos(2*theta)*(-2 + 2*r + beta_m))*epsilon0**gamma

def Ccoef(r, theta, l0, beta_m, epsilon0, K, gamma):
	return l0**2*(-2 + r)*beta_m*epsilon0 + K*gamma*(-2*r**3*np.sin(theta)**2 + l0**3*(-2 + r)*(2 + beta_m))*epsilon0**gamma

def Ccoef_cart(x, y, l0, beta_m, epsilon0, K, gamma):
	r, theta = cart_to_polar(x,y)

	return l0**2*(-2 + r)*beta_m*epsilon0 + K*gamma*(-2*r**3*np.sin(theta)**2 + l0**3*(-2 + r)*(2 + beta_m))*epsilon0**gamma

def f1_tilde(r, theta, l0):

	return 32*l0**4*(-3 + r)*(-2 + r)**2 + 16*l0**6*(1 - 2/r)**4*r + 10*r**6 - 6*l0**2*(-2 + r)*r**3*(-12 + 5*r) + (-32*l0**4*(-3 + r)*(-2 + r)**2 - 15*r**6 + 8*l0**2*r**3*(-12 + 5*r))*np.cos(2*theta)

def f1_tilde_cart(x, y, l0):
	r, theta = cart_to_polar(x,y)

	return 32*l0**4*(-3 + r)*(-2 + r)**2 + 16*l0**6*(1 - 2/r)**4*r + 10*r**6 - 6*l0**2*(-2 + r)*r**3*(-12 + 5*r) + (-32*l0**4*(-3 + r)*(-2 + r)**2 - 15*r**6 + 8*l0**2*r**3*(-12 + 5*r))*np.cos(2*theta)

def f2_tilde(r, theta, l0):

	return 4*(1 - 2/r)*r**3*(4*l0**4*(-2 + r)**2 - 2*l0**2*(-2 + r)*r**3 + 3*r**6 + 2*r**3*(-(l0**2*(-2 + r)) + 2*r**3)*np.cos(2*theta) + r**6*np.cos(4*theta))

def f2_tilde_cart(x, y, l0):
	r, theta = cart_to_polar(x,y)

	return 4*(1 - 2/r)*r**3*(4*l0**4*(-2 + r)**2 - 2*l0**2*(-2 + r)*r**3 + 3*r**6 + 2*r**3*(-(l0**2*(-2 + r)) + 2*r**3)*np.cos(2*theta) + r**6*np.cos(4*theta))

def partial_r_term(r, theta, l0, beta_m):
	return ((-2 + r)*r*np.cos(theta)*(r**3 + l0**2*(2 - r)*(1/(np.sin(theta)))**2)**3*np.sin(theta)**5*(1 + beta_m))/beta_m

def partial_theta_term(r, theta, l0, beta_m, epsilon0, K, gamma):
	return r * (-4*(r**3 + l0**2*(2 - r)*(1/(np.sin(theta)))**2)**3*k1_func(r, theta, l0, beta_m, epsilon0, K, gamma)*np.sin(theta)**6*(1 + beta_m))/(Ccoef(r, theta, l0, beta_m, epsilon0, K, gamma)*beta_m)

def partial_x_term(x, y, l0, beta_m, epsilon0, K, gamma):
	r, theta = cart_to_polar(x,y)
	return partial_r_term(r, theta, l0, beta_m) * np.sin(theta) + partial_theta_term(r, theta, l0, beta_m, epsilon0, K, gamma) * np.cos(theta)

def partial_y_term(x, y, l0, beta_m, epsilon0, K, gamma):
	r, theta = cart_to_polar(x,y)
	return partial_r_term(r, theta, l0, beta_m) * np.cos(theta) - partial_theta_term(r, theta, l0, beta_m, epsilon0, K, gamma) * np.sin(theta)

def indep_term(r, theta, l0, beta_m, epsilon0, K, gamma, m1, m2, tau2):
	return 2*l0**2*m1*tau2*(1/(np.tan(theta)))*(A_tilde(r, theta, l0) + (B_tilde(r, theta, l0)*k1_func(r, theta, l0, beta_m, epsilon0, K, gamma))/Ccoef(r, theta, l0, beta_m, epsilon0, K, gamma)) + (m2*(1/(np.tan(theta)))*(f1_tilde(r, theta, l0) - (f2_tilde(r, theta, l0)*k1_func(r, theta, l0, beta_m, epsilon0, K, gamma))/Ccoef(r, theta, l0, beta_m, epsilon0, K, gamma)))/2.

def indep_term_cart(x, y, l0, beta_m, epsilon0, K, gamma, m1, m2, tau2):
	r, theta = cart_to_polar(x,y)
	return 2*l0**2*m1*tau2*(1/(np.tan(theta)))*(A_tilde(r, theta, l0) + (B_tilde(r, theta, l0)*k1_func(r, theta, l0, beta_m, epsilon0, K, gamma))/Ccoef(r, theta, l0, beta_m, epsilon0, K, gamma)) + (m2*(1/(np.tan(theta)))*(f1_tilde(r, theta, l0) - (f2_tilde(r, theta, l0)*k1_func(r, theta, l0, beta_m, epsilon0, K, gamma))/Ccoef(r, theta, l0, beta_m, epsilon0, K, gamma)))/2.

def RK_RHS(x, y, l0, beta_m, epsilon0, K, gamma):
	return (partial_x_term(x, y, l0, beta_m, epsilon0, K, gamma))/(partial_y_term(x, y, l0, beta_m, epsilon0, K, gamma))

def RK_RHS_polar(r, theta, l0, beta_m, epsilon0, K, gamma):
	return (partial_r_term(r, theta, l0, beta_m))/(partial_theta_term(r, theta, l0, beta_m, epsilon0, K, gamma))

def RK_RHS_full_ev(beta_mc, r_c, x, y, gamma, l, K, M, W_in):
	beta_m = beta_m_cart(beta_mc, r_c, x, y, gamma)
	W = W_function_cart(l, x, y, M)
	p0_val = p0(K, beta_m, W, W_in, gamma)
	epsilon0_val = epsilon0(p0_val, K, gamma)
	return RK_RHS(x, y, l, beta_m, epsilon0_val, K, gamma)

def RK_RHS_full_ev_polar(beta_mc, r_c, r, theta, gamma, l, K, M, W_in):
	beta_m_val = beta_m(beta_mc, r_c, r, theta, gamma)
	W = W_function(l, r, theta, M)
	p0_val = p0(K, beta_m_val, W, W_in, gamma)
	epsilon0_val = epsilon0(p0_val, K, gamma)
	return RK_RHS_polar(r, theta, l, beta_m_val, epsilon0_val, K, gamma)

def p_integrand(x, y, l, rc, K, gamma, m1, m2, tau2, W_in, M, beta_mc):
	beta_m_val = beta_m_cart(beta_mc, rc, x, y, gamma)
	W = W_function_cart(l, x, y, M)
	p0_val = p0(K, beta_m_val, W, W_in, gamma)
	epsilon0_val = epsilon0(p0_val, K, gamma)
	f = -indep_term_cart(x, y, l, beta_m_val, epsilon0_val, K, gamma, m1, m2, tau2)
	return f / partial_y_term(x, y, l, beta_m_val, epsilon0_val, K, gamma)

def obtain_p0(l, x, y, M, rc, beta_mc, gamma, K, Win):
	W_val = W_function_cart(l, x, y, M)
	beta_m_val = beta_m_cart(beta_mc, rc, x, y, gamma)
	return p0(K, beta_m_val, W_val, Win, gamma)

def caligraphic_l_function_eq(x, M):
	return (x**2 - 2*M*x)

def g_tt(r, theta, M):
	return - (1 - (2 * M)/(r))

def g_tphi(r, theta, M):
	return 0.

def g_phiphi(r, theta):
	return (r ** 2) * (np.sin(theta) ** 2)

def Omega(g_tt, g_tphi, g_phiphi, l):
	return - (g_tphi + g_tt * l)/(g_phiphi + g_tphi * l)

def utcov(g_tt, g_tphi, g_phiphi, l):
	return -np.sqrt(np.abs((g_tphi ** 2 - g_tt * g_phiphi)/(g_phiphi + 2 * g_tphi * l + g_tt * (l) ** 2)))

def utcon(utcov, Omega, l):
	return - (1)/(utcov * (1 - Omega * l))

def uphicon(utcon, Omega):
	return Omega * utcon

def calligraphic_l_function(g_tt, g_tphi, g_phiphi):
	return g_tphi ** 2 - g_tt * g_phiphi

def Bphi(p0_val, p1, current_beta_m, current_g_tt, calligraphic_l):
	return float(np.sqrt(-current_g_tt)) * float(np.sqrt((2 * (p0_val + p1))/(current_beta_m * calligraphic_l)))

def compute_potential(u_tcov):
	return np.log(- u_tcov)

def compute_e0(p0_val, K, gamma):
	return (p0_val / K) ** (1/gamma)

def compute_e1(p1, e0, gamma, p0_val):
	return (p1 * e0) / (gamma * p0_val)

def compute_rho_atm(rho_min, r):
	return rho_min * (r ** (-3./2.))

def compute_p_atm(p_min, r):
	return p_min * (r ** (-5./2.))


