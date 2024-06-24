import metric
import matplotlib.pyplot as plt
from scipy.optimize import newton
from diff_eqs_coeffs import *
from scipy.interpolate import interp1d
# from scipy.interpolate import bisplrep
# from scipy.interpolate import bisplev
# from scipy.interpolate import RBFInterpolator
# from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import RectBivariateSpline
# from min_theta_calc import theta_min_calculator
from scipy.interpolate import griddata
# from scipy.spatial import KDTree
import json
from epsilon_computer import *



def to_obtain_cusp(r, M, a, l, s0, gamma, kappa):
	theta = np.pi/2
	return metric.F0_0(kappa, gamma)*drlogut.drlogut_f(M, a, r, theta, l)-F1ast.f1ast_f(M, a, r, theta, l, s0, kappa, gamma)-Gast.gast_f(M, a, r, theta, l, s0, kappa, gamma)


def r_mb(M, a):
	return 2*M*(1-a/2+np.sqrt(1-a))

def radial_epsilon(r, epsilon, M, a, l, s0, kappa, gamma):
	return F0(epsilon, kappa, gamma)*drlogut.drlogut_f(M, a, r, np.pi/2, l) - F1ast.f1ast_f(M, a, r, np.pi/2, l, s0, kappa, gamma)-Gast.gast_f(M, a, r, np.pi/2, l, s0, kappa, gamma)

with open("parameters.json") as param_file:
    parameters = json.load(param_file)

M = float(parameters["M"])
a = float(parameters["a"])
l = float(parameters["l"])
s0 = float(parameters["s0"])
kappa = float(parameters["kappa"])
gamma = float(parameters["gamma"])
r_cusp_0 = r_mb(M, a)

r_cusp = newton(to_obtain_cusp, r_cusp_0, args=(M, a, l, s0, gamma, kappa), maxiter=1000)
print(r_cusp)

h = 0.005
imax = 30/h
current_epsilon = 1e-20
current_r = r_cusp
r_array = []
epsilon_array = []
i = 0


while i<imax:
	r_array.append(current_r)
	epsilon_array.append(current_epsilon)
	pre_r = current_r
	pre_epsilon = current_epsilon
	k1 = h*radial_epsilon(pre_r, pre_epsilon, M, a, l, s0, kappa, gamma)
	k2 = h*radial_epsilon(pre_r + h/2., pre_epsilon+k1/2., M, a, l, s0, kappa, gamma)
	k3 = h*radial_epsilon(pre_r + h/2., pre_epsilon+k2/2., M, a, l, s0, kappa, gamma)
	k4 = h*radial_epsilon(pre_r + h, pre_epsilon+k3, M, a, l, s0, kappa, gamma)
	current_epsilon = pre_epsilon + (1./6.)*(k1+2*k2+2*k3+k4)
	current_r = pre_r + h
	i = i + 1
	if (current_epsilon<1e-20):
		break

epsilon_array = np.asarray(epsilon_array)
r_array = np.asarray(r_array)

epsilon_max = epsilon_array.max()
r_max = r_array[np.where(epsilon_array == epsilon_max)]

r_out = r_array[-1]
print(r_out)
print(r_max, epsilon_max)

plt.plot(r_array, epsilon_array, c = 'black', linewidth = 3)
# plt.show()
plt.close()

boundary_array_k1 = epsilon_computer_f(M, a, l, s0, kappa, gamma, True, 0, 0)
boundary_y_array_k1 = boundary_array_k1[:,1]
boundary_x_array_k1 = boundary_array_k1[:,0]
y_max = 1.10*np.max(boundary_y_array_k1)
# y_max = 1.1*max_locations[1]
# theta_min = np.arctan2(r_max, y_max)

# plt.plot(boundary_x_array_k1, boundary_y_array_k1, c = 'black', linewidth = 3)
# plt.show()
# plt.close()

######## with a correction factor for y of 1.10
##This value could affect an analysis based only on the shape of the characteristics of k, be careful!!!!
boundary_func = interp1d(boundary_x_array_k1, 1.10*boundary_y_array_k1, kind='cubic')
########


x_precision = 0.025
y_precision = 0.005

# theta_min = theta_min_calculator(M, a, l, s0, gamma, kappa, x_precision, y_precision, y_max)

# print(theta_min)

k0 = 1.

N_x = int((r_out - r_cusp)/x_precision)
N_y = int((y_max)/y_precision)

x_eq_array = np.linspace(r_cusp + x_precision, r_out, num=N_x, endpoint=False)
print(x_eq_array)
print(len(x_eq_array))
y_search_array = np.linspace(0, y_max, num=30*N_y)

initial_data_list = []

i = 0
i_max = len(x_eq_array)

while i < i_max:
	initial_data_list.append([x_eq_array[i],1e-10])
	i = i + 1


i = 0
h = y_precision

curves_grid = []

print("computing characteristics")


for initial_point in initial_data_list:
	current_x = initial_point[0]
	current_y = initial_point[1]
	current_curve_array = []
	current_curve_array.append([current_x, current_y])
	current_curve_x = []
	current_curve_y = []
	current_curve_x.append(current_x)
	current_curve_y.append(current_y)
	current_theta = np.pi/2

	while (current_y < boundary_func(current_x)):
		try:
			k1 = Ak(M, a, current_x, current_y, l, s0, kappa, gamma)
		except ValueError:
			current_curve_array = np.asarray(current_curve_array)
			curves_grid.append(current_curve_array)
			plt.plot(current_curve_x, current_curve_y, color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
			break
		try:
			k2 = Ak(M, a, current_x + (1./2.)*(k1*h), current_y + (1./2.)*h, l, s0, kappa, gamma)
		except ValueError:
			current_curve_array = np.asarray(current_curve_array)
			curves_grid.append(current_curve_array)
			plt.plot(current_curve_x, current_curve_y, color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
			break
		try:
			k3 = Ak(M, a, current_x + (1./2.)*(k2*h), current_y + (1./2.)*h, l, s0, kappa, gamma)
		except ValueError:
			current_curve_array = np.asarray(current_curve_array)
			curves_grid.append(current_curve_array)
			plt.plot(current_curve_x, current_curve_y, color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
			break
		try:
			k4 = Ak(M, a, current_x + (k3*h), current_y + h, l, s0, kappa, gamma)
		except ValueError:
			current_curve_array = np.asarray(current_curve_array)
			curves_grid.append(current_curve_array)
			plt.plot(current_curve_x, current_curve_y, color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
			break
		if (k1 != k1 or k2 != k2 or k3 != k3 or k4 != k4):
				current_curve_array = np.asarray(current_curve_array)
				curves_grid.append(current_curve_array)
				plt.plot(current_curve_x, current_curve_y, color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
				break

		current_x = current_x + (h/6.)*(k1 + 2*k2 + 2*k3 + k4)
		current_y = current_y + h
		current_theta = np.arctan2(current_x, current_y)
		current_curve_array.append([current_x, current_y])
		current_curve_x.append(current_x)
		current_curve_y.append(current_y)

	else:
		current_curve_array = np.asarray(current_curve_array)
		curves_grid.append(current_curve_array)
		plt.plot(current_curve_x, current_curve_y, color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])

	print(i)
	i = i + 1

curves_grid = np.asarray(curves_grid)

# np.savetxt("curves_grid_0.dat", curves_grid[0])
# np.savetxt("curves_grid_1.dat", curves_grid[1])
# np.savetxt("curves_grid_2.dat", curves_grid[2])


plt.savefig("characteristics.pdf")
plt.close()

# total_data_array = []

# j = 0
# i = 0

# print(len(curves_grid))

# for curve in curves_grid:

# 	# print(curve)
# 	integrand_array = []
# 	k_initial = k0
# 	current_k = k_initial
# 	curve_k_array = []
# 	print(j)
# 	j = j + 1
# 	i = 0
# 	bad_results_flag = False

# 	for point in curve:
# 		integrand_array.append(Ck(M, a, point[0], point[1], l, s0, kappa, gamma, current_k))

# 	current_curve_x = curve[:,0]
# 	current_curve_y = curve[:,1]
# 	current_x_spline = interp1d(current_curve_y, current_curve_x, kind='cubic')
# 	current_x = current_curve_x[0]
# 	current_y = current_curve_y[0]

# 	for current_y in current_curve_y:
# 		curve_k_array.append(current_k)
# 		total_data_array.append([current_x, current_y, current_k])

# 		if (current_y > boundary_func(current_x)):
# 			break

# 		if (current_y == current_curve_y[-1]):
# 			break
# 		else:
# 			try:
# 				k1_2 = Ck(M, a, current_x, current_y, l, s0, kappa, gamma, current_k)
# 			except ValueError:
# 				break
# 			current_y_2 = current_y + h/2.
# 			current_x_2 = current_x_spline(current_y_2)
# 			try:
# 				k2_2 = Ck(M, a, current_x_2, current_y_2, l, s0, kappa, gamma, current_k + h*k1_2/2.)
# 			except ValueError:
# 				break
# 			try:
# 				k3_2 = Ck(M, a, current_x_2, current_y_2, l, s0, kappa, gamma, current_k + h*k2_2/2.)
# 			except ValueError:
# 				break
# 			current_y_4 = current_y + h
# 			current_x_4 = current_x_spline(current_y_4)
# 			try:
# 				k4_2 = Ck(M, a, current_x_4, current_y_4, l, s0, kappa, gamma, current_k + h*k3_2)
# 			except ValueError:
# 				break
# 			if (k1_2 != k1_2 or k2_2 != k2_2 or k3_2 != k3_2 or k4_2 != k4_2):
# 				break

# 			if ((current_k + (h/6.)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)) < 0 or bad_results_flag == True):
# 				bad_results_flag = True
# 				current_x = current_x_spline(current_y)

# 			else:
# 				current_k = current_k + (h/6.)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)
# 				current_x = current_x_spline(current_y)

# 			i = i + 1
		

# total_data_array = np.asarray(total_data_array)

# np.savetxt("output_k.dat", total_data_array)

# r_data_array = np.sqrt(total_data_array[:,0]**2+total_data_array[:,1]**2)
# theta_data_array = np.arctan2(total_data_array[:,0], total_data_array[:,1])

# print(np.shape(r_data_array), np.shape(theta_data_array))

# polar_total_data_grid = np.array([r_data_array, theta_data_array]).T
# # print(polar_total_data_grid)

# # print(np.shape(polar_total_data_grid))

# # k_neighbors_tree = KDTree(polar_total_data_grid)

# # points = cart_to_polar(np.array([9.8, 1.807]))

# # print(k_neighbors_tree.query(np.array(points), k=4))


# r_min_int = np.min(r_data_array)
# r_max_int = np.max(r_data_array)

# theta_min_int = np.min(theta_data_array)
# theta_max_int = np.max(theta_data_array)

# r_step = int(((r_max_int-r_min_int)/(0.001)))*1j + 1j
# theta_step = int(((theta_max_int-theta_min_int)/(0.001)))*1j + 1j

# # print(r_min_int, r_max_int, theta_min_int, theta_max_int)
# # print(r_data_array)
# # print(theta_data_array)
# # print(polar_total_data_grid)
# # print(total_data_array[:,2])

# grid_r, grid_theta = np.mgrid[r_min_int:r_max_int:r_step, theta_min_int:theta_max_int:theta_step]
# # print(grid_r)
# # print(grid_theta)

# data_to_inter = total_data_array[:,2]
# # data_to_inter -= 1.0
# grid_k = griddata(polar_total_data_grid, total_data_array[:,2], (grid_r, grid_theta), method='cubic', rescale=True, fill_value=1.0)
# # print(grid_k)


# r_1d_array = np.linspace(r_min_int, r_max_int, num=int(np.imag(r_step)))
# theta_1d_array = np.linspace(theta_min_int, theta_max_int, num=int(np.imag(theta_step)))

# # k_interpolated = RectBivariateSpline(r_1d_array, theta_1d_array, grid_k, kx=1, ky=1)
# k_interpolated = RectBivariateSpline(r_1d_array, theta_1d_array, grid_k, kx=3, ky=3)

# def k_interpolated_cart(x, y):
# 	polar_coords = cart_to_polar([x, y])
# 	return k_interpolated.ev(polar_coords[0], polar_coords[1])



# ###################### PLOT K #######################
# # x_step = 0.01
# # y_step = 0.01

# # # print(r_min_int, r_max_int, theta_min_int, theta_max_int)
# # # print(r_data_array)
# # # print(theta_data_array)
# # # print(polar_total_data_grid)
# # # print(total_data_array[:,2])

# # grid_r, grid_y = np.mgrid[r_min_int:r_max_int:x_step, theta_min_int:theta_max_int:y_step]
# # # print(grid_r)
# # # print(grid_theta)

# # data_to_inter = np.concatenate((total_data_array[:,2], total_data_array[:,2]), axis=0)
# # # data_to_inter -= 1.0
# # grid_k = griddata(cart_total_data_grid, data_to_inter, (grid_x, grid_y), method='cubic', fill_value=0.0)
# # # print(grid_k)


# # xArray = full_x_coords
# # yArray = full_y_coords

# # # xi = np.linspace(np.min(xArray), np.max(xArray),5000)
# # # yi = np.linspace(np.min(yArray), np.max(yArray),5000)

# # DEM = griddata((xArray, yArray),data_to_inter, (grid_x, grid_y), method='cubic')

# # f = plt.figure(figsize=(13,11.5))

# # im = plt.imshow(DEM.transpose(), cmap='plasma', vmin=np.min(grid_k), vmax=np.min(grid_k),
# #         extent=[np.min(grid_x), np.max(grid_x), np.min(grid_y), np.max(grid_y)],
# #         interpolation='nearest', origin='lower')
# # # f.colorbar(im)


# # # plt.axis([0, x_max, -y_max, y_max])
# # # plt.axis([-1000, 1010, -1005, 1005])
# # # plb.ylim([0,r_out/2])
# # # plb.xlim([0,r_out*1.2])
# # plb.ylim([-7.5,7.5])
# # plb.xlim([0,15])

# # plt.tick_params(labelsize=35, pad = 12)
# # lx = plt.xlabel("$r$ sin $\\theta$", fontsize=40)
# # ly = plt.ylabel("$r$ cos $\\theta$", fontsize=40)
# # if (no_k_flag == True):
# # 	plt.savefig("epsilon_no_k.pdf")
# # else:
# # 	plt.savefig("epsilon.pdf")

# # plt.close()

# ####################################################

# # print(k_interpolated_cart(9.8, 1.806))
# # print(k_interpolated_cart(9.769126992296294887, 1.805000000099983515))
# # print(k_interpolated.ev(10, np.pi/2.1))
# # print(k_interpolated.ev(np.sqrt(15**2+6**2), np.arctan2(15, 6)))


# epsilon_computer_f(M, a, l, s0, kappa, gamma, False, k_interpolated_cart, boundary_func)

# # k_interpolated = SmoothBivariateSpline(r_data_array, theta_data_array, total_data_array[:,2], kx=1, ky=1)




# # print(k_interpolated_cart(9.769126992296294887, 1.805000000099983515))

# # 9.769126992296294887 e + 00 1.805000000099983515 e +00 1.033366028650707147 e + 00

# # y_search_array = np.linspace(0, y_max, num=30*N_y)

# # print(len(r_data_array), len(theta_data_array), len(total_data_array[:,2]))
# # print(r_data_array)
# # print(theta_data_array)
# # print(total_data_array[:,2])




# # # 9.368042994858670980e+00 1.775000000099973274e+00 9.635292310449019570e-01


# # print(len(total_data_array[:,0]), len(total_data_array[:,1]), len(total_data_array[:,2]), len(x_eq_array), len(y_search_array))



# # print(k_interpolated.ev(8.731705734213683456, 6.245000000099889093))
# #8.731705734213683456e+00 6.245000000099889093e+00 1.051271513145822789e-01
# # k_interpolated = RBFInterpolator(total_data_array[:,:2],total_data_array[:,2])
# # print(k_interpolated(np.asarray([8.731705734213683456, 6.245000000099889093])))

# # k_interpolated = SmoothBivariateSpline(total_data_array[:,0], total_data_array[:,1],total_data_array[:,2], s=1)
# # print(k_interpolated.ev(5.882688826865476450, 0.8400000001000006433))
# # k_interpolated = bisplrep(total_data_array[:,0], total_data_array[:,1],total_data_array[:,2])

# # print(bisplev(5.372034947596807442, 3.955000000099937907, k_interpolated))
# # print(bisplev(5.882688826865476450, 0.8400000001000006433, k_interpolated))
# # # 2.194168934260585235e-01
# # 9.767229419719847261e-01
# # print(k_interpolated.ev(5.882688826865476450, 0.8400000001000006433))







