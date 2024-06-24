import numpy as np

def g_tt(M, a, r, theta):
	return -1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)

def g_tphi(M, a, r, theta):
	return (-2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)

def g_phiphi(M, a, r, theta):
	return np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))

def g_rr(M, a, r, theta):
	return (r**2 + a**2*np.cos(theta)**2)/(a**2 - 2*M*r + r**2)

def g_thetatheta(M, a, r, theta):
	return r**2 + a**2*np.cos(theta)**2

def g_ttcon(M, a, r, theta):
	return g_phiphi(M, a, r, theta)/(-g_tphi(M, a, r, theta)**2 + g_tt(M, a, r, theta)*g_phiphi(M, a, r, theta))

def g_tphicon(M, a, r, theta):
	return g_tphi(M, a, r, theta)/(g_tphi(M, a, r, theta)**2 - g_tt(M, a, r, theta)*g_phiphi(M, a, r, theta))

def g_phiphicon(M, a, r, theta):
	return g_tt(M, a, r, theta)/(-g_tphi(M, a, r, theta)**2 + g_tt(M, a, r, theta)*g_phiphi(M, a, r, theta))

def Christoffel_rtt(M, a, r, theta):
	return (M*(a**2 + r*(-2*M + r))*(r**2 - a**2*np.cos(theta)**2))/(r**2 + a**2*np.cos(theta)**2)**3

def Christoffel_rtphi(M, a, r, theta):
	return (a*M*(a**2 + r*(-2*M + r))*(-r**2 + a**2*np.cos(theta)**2)*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**3

def Christoffel_rphiphi(M, a, r, theta):
	return ((a**2 + r*(-2*M + r))*(-2*r*np.sin(theta)**2 + (2*a**2*M*(r**2 - a**2*np.cos(theta)**2)*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2))/(2*(r**2 + a**2*np.cos(theta)**2))

def Christoffel_thetatt(M, a, r, theta):
	return (-2*a**2*M*r*np.cos(theta)*np.sin(theta))/(r**2 + a**2*np.cos(theta)**2)**3

def Christoffel_thetatphi(M, a, r, theta):
	return (a*M*r*(a**2 + r**2)*np.sin(2*theta))/(r**2 + a**2*np.cos(theta)**2)**3

def Christoffel_thetaphiphi(M, a, r, theta):
	return (np.cos(theta)*np.sin(theta)*(-a**2 - r**2 - (4*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) - (2*a**4*M*r*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2))/(r**2 + a**2*np.cos(theta)**2)

def Christoffel_ttr(M, a, r, theta):
	return (-2*M*(a**2 + r**2)*(a**2 - 2*r**2 + a**2*np.cos(2*theta)))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)

def Christoffel_tttheta(M, a, r, theta):
	return (-4*a**2*M*r*np.sin(2*theta))/(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2

def Christoffel_tphir(M, a, r, theta):
	return (2*a*M*(a**4 - 3*a**2*r**2 - 6*r**4 + a**2*(a**2 - r**2)*np.cos(2*theta))*np.sin(theta)**2)/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)

def Christoffel_tphitheta(M, a, r, theta):
	return (8*a**3*M*r*np.cos(theta)*np.sin(theta)**3)/(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2

def Christoffel_phitr(M, a, r, theta):
	return (4*a*M*(r**2 - a**2*np.cos(theta)**2))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)

def Christoffel_phittheta(M, a, r, theta):
	return (-8*a*M*r*(np.cos(theta)/np.sin(theta)))/(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2

def Christoffel_phiphir(M, a, r, theta):
	return (4*(r**4*(-2*M + r) + a**4*r*np.cos(theta)**4 - a**2*M*r**2*np.sin(theta)**2 + np.cos(theta)**2*(2*a**2*r**2*(-M + r) + a**4*M*np.sin(theta)**2)))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)

def Christoffel_phiphitheta(M, a, r, theta):
	return ((3*a**4 + 8*a**2*M*r + 8*a**2*r**2 + 8*r**4 + 4*a**2*(a**2 + 2*r*(-M + r))*np.cos(2*theta) + a**4*np.cos(4*theta))*(np.cos(theta)/np.sin(theta)))/(2*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)

def R_trphir(M, a, r, theta):
	return (12*a*M*r*(a**2 + r**2)*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(3*a**2 - 2*r**2 + 3*a**2*np.cos(2*theta))*np.sin(theta)**2)/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (4*a*M**2*r**2*(-9*a**4 + 2*a**2*r**2 + 4*r**4 - 2*a**2*(3*a**2 + 5*r**2)*np.cos(2*theta) + 3*a**4*np.cos(4*theta))*np.sin(theta)**2)/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)

def R_trphitheta(M, a, r, theta):
	return (a*M*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)

def R_rphiphir(M, a, r, theta):
	return (-8*M*r*(r**2 + a**2*np.cos(theta)**2)*(3*a**2 - 2*r**2 + 3*a**2*np.cos(2*theta))*(2*a**4 + r**4 + a**2*r*(-2*M + 3*r) - a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(theta)**2)/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4)

def R_rphiphitheta(M, a, r, theta):
	return (24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3)/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4)

def R_trtr(M, a, r, theta):
	return (M*r*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(-27*a**4 - 4*a**2*r**2 + 16*r**4 - 4*a**2*(6*a**2 + 7*r**2)*np.cos(2*theta) + 3*a**4*np.cos(4*theta)))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (24*a**2*M**2*r**2*(3*a**2 - 2*r**2 + 3*a**2*np.cos(2*theta))*np.sin(theta)**2)/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)

def R_trttheta(M, a, r, theta):
	return (8*a**2*M**2*r*np.cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*np.cos(2*theta))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (2*a**2*M*(3*a**2 + r*(-2*M + 3*r))*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)

def R_rphitr(M, a, r, theta):
	return (a*M*r*(3*a**2 + r*(-4*M + 3*r))*(3*a**2 - 2*r**2 + 3*a**2*np.cos(2*theta))*np.sin(theta)**2)/(2*(a**2 - 2*M*r + r**2)*(r**2 + a**2*np.cos(theta)**2)**3)

def R_rphittheta(M, a, r, theta):
	return (2*a*M*(a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*np.cos(2*theta) + a**4*np.cos(4*theta))*np.sin(2*theta))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4)

def detg(M, a, r, theta):
	return -1/4*(((a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)*np.sin(theta)**2)

def utcov(M, a, r, theta, l):
	return -np.sqrt(np.abs(((g_tphi(M, a, r, theta))**2 - g_phiphi(M, a, r, theta)*g_tt(M, a, r, theta))/(g_phiphi(M, a, r, theta) + 2*g_tphi(M, a, r, theta)*l + g_tt(M, a, r, theta)*(l**2))))

def uphicov(M, a, r, theta, l):
	return - l * utcov(M, a, r, theta, l)

def utcon(M, a, r, theta, l):
	return g_ttcon(M, a, r, theta)*utcov(M, a, r, theta, l)+g_tphicon(M, a, r, theta)*uphicov(M, a, r, theta, l)

def uphicon(M, a, r, theta, l):
	return g_tphicon(M, a, r, theta)*utcov(M, a, r, theta, l)+g_phiphicon(M, a, r, theta)*uphicov(M, a, r, theta, l)

def Omega_func(M, a, r, theta, l):
	return -(g_tphi(M, a, r, theta)+g_tt(M, a, r, theta)*l)/(g_phiphi(M, a, r, theta)+g_tphi(M, a, r, theta)*l)

def Potential_func(M, a, r, theta, l):
	return np.log(-utcov(M, a, r, theta, l))

def a_rcon(M, a, r, theta, l):
	return (utcon(M, a, r, theta, l)**2)*Christoffel_rtt(M, a, r, theta) + 2*utcon(M, a, r, theta, l)*uphicon(M, a, r, theta, l)*Christoffel_rtphi(M, a, r, theta)+(uphicon(M, a, r, theta, l)**2)*Christoffel_rphiphi(M, a, r, theta)

def a_thetacon(M, a, r, theta, l):
	return (utcon(M, a, r, theta, l)**2)*Christoffel_thetatt(M, a, r, theta) + 2*utcon(M, a, r, theta, l)*uphicon(M, a, r, theta, l)*Christoffel_thetatphi(M, a, r, theta)+(uphicon(M, a, r, theta, l)**2)*Christoffel_thetaphiphi(M, a, r, theta)

def Datcon(M, a, r, theta, l):
	return utcon(M, a, r, theta, l)*a_rcon(M, a, r, theta, l)*Christoffel_ttr(M, a, r, theta)+utcon(M, a, r, theta, l)*a_thetacon(M, a, r, theta, l)*Christoffel_tttheta(M, a, r, theta)+uphicon(M, a, r, theta, l)*a_rcon(M, a, r, theta, l)*Christoffel_tphir(M, a, r, theta)+uphicon(M, a, r, theta, l)*a_thetacon(M, a, r, theta, l)*Christoffel_tphitheta(M, a, r, theta)

def Daphicon(M, a, r, theta, l):
	return utcon(M, a, r, theta, l)*a_rcon(M, a, r, theta, l)*Christoffel_phitr(M, a, r, theta)+utcon(M, a, r, theta, l)*a_thetacon(M, a, r, theta, l)*Christoffel_phittheta(M, a, r, theta)+uphicon(M, a, r, theta, l)*a_rcon(M, a, r, theta, l)*Christoffel_phiphir(M, a, r, theta)+uphicon(M, a, r, theta, l)*a_thetacon(M, a, r, theta, l)*Christoffel_phiphitheta(M, a, r, theta)

def F1(M, a, r, theta, l):
	return -Omega_func(M, a, r, theta, l)*l*R_trphir(M, a, r, theta)+Omega_func(M, a, r, theta, l)*R_rphiphir(M, a, r, theta)-l*R_trtr(M, a, r, theta)+R_rphitr(M, a, r, theta)

def F2(M, a, r, theta, l):
	return -Omega_func(M, a, r, theta, l)*l*R_trphitheta(M, a, r, theta)+Omega_func(M, a, r, theta, l)*R_rphiphitheta(M, a, r, theta)-l*R_trttheta(M, a, r, theta)+R_rphittheta(M, a, r, theta)

def drlogut(M, a, r, theta, l):
	return ((l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))*(((-16*a**2*M**2*r**3*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**3 + (8*a**2*M**2*r*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(2*r - (4*a**2*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a**2*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)) - ((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/(l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))) - ((l**2*((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2)) + (8*a*l*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 - (4*a*l*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(2*r - (4*a**2*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a**2*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))*((4*a**2*M**2*r**2*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))))/(l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))**2))/(2*((4*a**2*M**2*r**2*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))))

def dthetalogut(M, a, r, theta, l):
	return ((l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))*(-((((4*a**2*M**2*r**2*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))*((4*a**2*l**2*M*r*np.cos(theta)*np.sin(theta))/(r**2 + a**2*np.cos(theta)**2)**2 - (8*a*l*M*r*np.cos(theta)*np.sin(theta))/(r**2 + a**2*np.cos(theta)**2) - (8*a**3*l*M*r*np.cos(theta)*np.sin(theta)**3)/(r**2 + a**2*np.cos(theta)**2)**2 + 2*np.cos(theta)*np.sin(theta)*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)) + np.sin(theta)**2*((4*a**2*M*r*np.cos(theta)*np.sin(theta))/(r**2 + a**2*np.cos(theta)**2) + (4*a**4*M*r*np.cos(theta)*np.sin(theta)**3)/(r**2 + a**2*np.cos(theta)**2)**2)))/(l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))**2) + ((16*a**2*M**2*r**2*np.cos(theta)*np.sin(theta)**3)/(r**2 + a**2*np.cos(theta)**2)**2 + (16*a**4*M**2*r**2*np.cos(theta)*np.sin(theta)**5)/(r**2 + a**2*np.cos(theta)**2)**3 - 2*np.cos(theta)*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)) - (4*a**2*M*r*np.cos(theta)*np.sin(theta)**3*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*((4*a**2*M*r*np.cos(theta)*np.sin(theta))/(r**2 + a**2*np.cos(theta)**2) + (4*a**4*M*r*np.cos(theta)*np.sin(theta)**3)/(r**2 + a**2*np.cos(theta)**2)**2))/(l**2*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2)) - (4*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))))/(2*((4*a**2*M**2*r**2*np.sin(theta)**4)/(r**2 + a**2*np.cos(theta)**2)**2 - (-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))))

def keq():
	return 1.0

def G1ast(M, a, r, theta, l, s0, kappa, gamma):
	return ((2*g_rr(M, a, r, theta)*g_tt(M, a, r, theta)*s0*(1-l*Omega_func(M, a, r, theta, l)))/(kappa*gamma))*np.sqrt(g_thetatheta(M, a, r, theta)/(-detg(M, a, r, theta)))

def F0_0(kappa, gamma):
	return 1./(-kappa*gamma)

def F1ast(M, a, r, theta, l, s0, kappa, gamma):
	return ((s0)/(kappa*gamma))*np.sqrt((g_thetatheta(M, a, r, theta))/(-detg(M, a, r, theta)))*F1(M, a, r, theta, l)

def F2ast(M, a, r, theta, l, s0, kappa, gamma):
	return ((s0)/(kappa*gamma))*np.sqrt((g_thetatheta(M, a, r, theta))/(-detg(M, a, r, theta)))*F2(M, a, r, theta, l)

def G2(M, a, r, theta, l):
	return ((uphicov(M, a, r, theta, l)-(g_tphi(M, a, r, theta)/g_tt(M, a, r, theta))*utcov(M, a, r, theta, l))*Datcon(M, a, r, theta, l) + ((g_tphi(M, a, r, theta)/g_tt(M, a, r, theta))*uphicov(M, a, r, theta, l)-(g_phiphi(M, a, r, theta)/g_tt(M, a, r, theta))*utcov(M, a, r, theta, l))*Daphicon(M, a, r, theta, l))

def Gast(M, a, r, theta, l, s0, kappa, gamma):
	return G1ast(M, a, r, theta, l, s0, kappa, gamma)*G2(M, a, r, theta, l)

def F0(epsilon, kappa, gamma):
	return (1+kappa*epsilon)/(-kappa*gamma)

def d_epsilon_F0ast(gamma):
	return -1./gamma

def d_rF2ast(M, a, r, theta, l, s0, kappa, gamma):
	return (2*s0*np.sqrt(((r**2 + a**2*np.cos(theta)**2)*(1./(np.sin(theta)))**2)/(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)*((-24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))*((4*a*l*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 - (2*a*l*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(2*r - (4*a**2*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a**2*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))**2) + (24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2))) - (4*a*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) - (384*a**2*M*r*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**5*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) - (288*a**2*M*r*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) + (48*a**2*M*r*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) + (24*a**2*M*(-2*M + 2*r)*(a**2 + r**2)*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) - (24*a**2*M*(-2*M + 2*r)*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)**2*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) + (48*a**2*M*r*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) + (2*a*M*(a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(44*a**2*r + 48*r**3 - 28*a**2*r*np.cos(2*theta))*np.sin(2*theta))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - (32*a*M*r*(a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*np.cos(2*theta) + a**4*np.cos(4*theta))*np.sin(2*theta))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**5) + (4*a*M*r*(a**2 + r*(-2*M + r))*(-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*np.cos(2*theta) + a**4*np.cos(4*theta))*np.sin(2*theta))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) + (2*a*M*(-2*M + 2*r)*(r**2 + a**2*np.cos(theta)**2)*(-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*np.cos(2*theta) + a**4*np.cos(4*theta))*np.sin(2*theta))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - (2*a*M*(-2*M + 2*r)*(a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*np.cos(2*theta) + a**4*np.cos(4*theta))*np.sin(2*theta))/((a**2 - 2*M*r + r**2)**2*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - l*((-96*a**2*M**2*r**2*np.cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*np.cos(2*theta))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - (96*a**2*M**2*r**2*np.cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (8*a**2*M**2*(2*M - 2*r)*r*np.cos(theta)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (16*a**2*M**2*r**2*np.cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*np.cos(2*theta))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)**2*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (8*a**2*M**2*r*(-2*M + 2*r)*np.cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*np.cos(2*theta))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))**2*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (8*a**2*M**2*np.cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*np.cos(2*theta))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (24*a**2*M*r*(3*a**2 + r*(-2*M + 3*r))*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - (24*a**2*M*r*(3*a**2 + r*(-2*M + 3*r))*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (2*a**2*M*(3*a**2 + r*(-2*M + 3*r))*((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)+ (2*a**2*M*(-2*M + 6*r)*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (2*a**2*M*(-2*M + 2*r)*(3*a**2 + r*(-2*M + 3*r))*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)) - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))*((a*M*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(a**2*(2*M - 9*r) - 9*a**2*r - 16*r**3 + a**2*(-2*M + 2*r)*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (12*a*M*r*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - (12*a*M*r*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (a*M*((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (a*M*(-2*M + 2*r)*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))**2*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (48*a**3*M**2*r**2*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - (48*a**3*M**2*r**2*(3*a**2 + r*(-2*M + 3*r))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (8*a**3*M**2*r**2*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)**2*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (4*a**3*M**2*r*(-2*M + 6*r)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) - (4*a**3*M**2*r*(-2*M + 2*r)*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))**2*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (4*a**3*M**2*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)))/((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))+ (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))*((4*a*l*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 - (2*a*l*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(2*r - (4*a**2*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a**2*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))*((a*M*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)))/((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))**2 - (l*(-(l*((-4*M*r**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*M)/(r**2 + a**2*np.cos(theta)**2))) - (4*a*M*r**2*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)**2 + (2*a*M*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))*((a*M*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)))/((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))))/(gamma*kappa) + (s0*((-8*r*(r**2 + a**2*np.cos(theta)**2)*(1./(np.sin(theta)))**2)/(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3 + (2*r*(1./(np.sin(theta)))**2)/(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2)*((24*a**2*M*(a**2 + r**2)*(a**2 + r*(-2*M + r))*np.cos(theta)*(r**2 + a**2*np.cos(theta)**2)*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**3*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4*((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))) + (2*a*M*(a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(-3*a**4 + 22*a**2*r**2 + 12*r**4 - 2*a**2*(a**2 + 7*r**2)*np.cos(2*theta) + a**4*np.cos(4*theta))*np.sin(2*theta))/((a**2 - 2*M*r + r**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**4) - l*((8*a**2*M**2*r*np.cos(theta)*(-2*a**2 + (2*M - r)*r + a**2*np.cos(2*theta))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (2*a**2*M*(3*a**2 + r*(-2*M + 3*r))*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)) - (l*(-(l*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))) + (2*a*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2))*((a*M*(-1 + (2*M*r)/(r**2 + a**2*np.cos(theta)**2))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*(-5*a**4 + a**2*(2*M - 9*r)*r - 4*r**4 + a**2*(a**2 + r*(-2*M + r))*np.cos(2*theta))*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3) + (4*a**3*M**2*r*(3*a**2 + r*(-2*M + 3*r))*(a**2 - 6*r**2 + a**2*np.cos(2*theta))*np.sin(theta)**2*np.sin(2*theta))/((a**2 + r*(-2*M + r))*(r**2 + a**2*np.cos(theta)**2)*(a**2 + 2*r**2 + a**2*np.cos(2*theta))**3)))/((-2*a*l*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2) + np.sin(theta)**2*(a**2 + r**2 + (2*a**2*M*r*np.sin(theta)**2)/(r**2 + a**2*np.cos(theta)**2)))))/(gamma*kappa*np.sqrt(((r**2 + a**2*np.cos(theta)**2)*(1./(np.sin(theta)))**2)/(a**2 + 2*r**2 + a**2*np.cos(2*theta))**2))


# print(R_trphir(1,0.5,10,np.pi/3))
# print(R_trphir(1,0.5,10,np.pi/2))

# print(Gast(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(Gast(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))

# print(G2(1, 0.5, 10, np.pi/2,3.8))
# print(G2(1, 0.5, 10, np.pi/3,3.8))

# print(d_rF2ast(1, 0.5, 10, np.pi/2,3.8, 0.25, 1, 2))
# print(d_rF2ast(1, 0.5, 10, np.pi/3,3.8, 0.25, 1, 2))













