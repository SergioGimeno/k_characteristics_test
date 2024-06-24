import metric
import matplotlib.pyplot as plt
from scipy.optimize import newton
from diff_eqs_coeffs import *
from scipy.interpolate import interp1d
# from scipy.interpolate import bisplrep
# from scipy.interpolate import bisplev
from scipy.interpolate import RBFInterpolator
from scipy.interpolate import SmoothBivariateSpline
from scipy.interpolate import RectBivariateSpline
# from min_theta_calc import theta_min_calculator
from scipy.interpolate import griddata
import matplotlib.pylab as plb
from scipy.spatial import KDTree

def epsilon_computer_f(M, a, l, s0, kappa, gamma, no_k_flag, k_func, boundary_func, y_max=20, theta_min=0.6):


	def to_obtain_cusp(r, M, a, l, s0, gamma, kappa):
		theta = np.pi/2
		return metric.F0_0(kappa, gamma)*drlogut.drlogut_f(M, a, r, theta, l)-F1ast.f1ast_f(M, a, r, theta, l, s0, kappa, gamma)-Gast.gast_f(M, a, r, theta, l, s0, kappa, gamma)


	def r_mb(M, a):
		return 2*M*(1-a/2+np.sqrt(1-a))

	def radial_epsilon(r, epsilon, M, a, l, s0, kappa, gamma):
		return F0(epsilon, kappa, gamma)*drlogut.drlogut_f(M, a, r, np.pi/2, l) - F1ast.f1ast_f(M, a, r, np.pi/2, l, s0, kappa, gamma)-Gast.gast_f(M, a, r, np.pi/2, l, s0, kappa, gamma)

	r_cusp_0 = r_mb(M, a)

	r_cusp = newton(to_obtain_cusp, r_cusp_0, args=(M, a, l, s0, gamma, kappa), maxiter=1000)
	print(r_cusp)

	h = 0.005
	imax =30/h
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

	np.savetxt("r_max_eps_max.dat", np.array([r_max, epsilon_max]))
	np.savetxt("r_cusp_r_out.dat", np.array([r_cusp, r_out]))

	# plt.plot(r_array, epsilon_array, c = 'black', linewidth = 3)
	# plt.show()
	# plt.close()


	x_precision = 0.02
	y_precision = 0.015

	# theta_min = theta_min_calculator(M, a, l, s0, gamma, kappa, x_precision, y_precision, y_max)

	# print(theta_min)

	# k0 = 1

	N_x = int((r_out - r_cusp)/x_precision)
	N_y = int((y_max)/y_precision)

	x_eq_array = np.linspace(r_cusp + x_precision, r_out, num=N_x, endpoint=False)
	# # print(x_eq_array)
	# # print(len(x_eq_array))
	# y_search_array = np.linspace(0, y_max, num=30*N_y)

	# initial_data_list = []

	# i = 0
	# i_max = len(x_eq_array)

	# while i < i_max:
	# 	initial_data_list.append([x_eq_array[i],1e-10])
	# 	i = i + 1


	# i = 0
	# h = y_precision

	# curves_grid = []

	# print("computing characteristics")


	# for initial_point in initial_data_list:
	# 	current_x = initial_point[0]
	# 	current_y = initial_point[1]
	# 	current_curve_array = []
	# 	current_curve_array.append([current_x, current_y])
	# 	current_curve_x = []
	# 	current_curve_y = []
	# 	current_curve_x.append(current_x)
	# 	current_curve_y.append(current_y)
	# 	current_theta = np.pi/2

	# 	while (current_y <= y_max and current_theta >= theta_min):
	# 		try:
	# 			k1 = Aepsilon(current_x, current_y)
	# 			# k1 = Ak(M, a, current_x, current_y, l, s0, kappa, gamma)
	# 		except ValueError:
	# 			current_curve_array = np.asarray(current_curve_array)
	# 			curves_grid.append(current_curve_array)
	# 			plt.plot(np.asarray(current_curve_x), np.asarray(current_curve_y), color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
	# 			break
	# 		try:
	# 			k2 = Aepsilon(current_x + (1./2.)*(k1*h), current_y + (1./2.)*h)
	# 			# k2 = Ak(M, a, current_x + (1./2.)*(k1*h), current_y + (1./2.)*h, l, s0, kappa, gamma)
	# 		except ValueError:
	# 			current_curve_array = np.asarray(current_curve_array)
	# 			curves_grid.append(current_curve_array)
	# 			plt.plot(np.asarray(current_curve_x), np.asarray(current_curve_y), color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
	# 			break
	# 		try:
	# 			k3 = Aepsilon(current_x + (1./2.)*(k2*h), current_y + (1./2.)*h)
	# 			# k3 = Ak(M, a, current_x + (1./2.)*(k2*h), current_y + (1./2.)*h, l, s0, kappa, gamma)
	# 		except ValueError:
	# 			current_curve_array = np.asarray(current_curve_array)
	# 			curves_grid.append(current_curve_array)
	# 			plt.plot(np.asarray(current_curve_x), np.asarray(current_curve_y), color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
	# 			break
	# 		try:
	# 			k3 = Aepsilon(current_x + (k3*h), current_y + h)
	# 			# k4 = Ak(M, a, current_x + (k3*h), current_y + h, l, s0, kappa, gamma)
	# 		except ValueError:
	# 			current_curve_array = np.asarray(current_curve_array)
	# 			curves_grid.append(current_curve_array)
	# 			plt.plot(np.asarray(current_curve_x), np.asarray(current_curve_y), color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
	# 			break
	# 		if (k1 != k1 or k2 != k2 or k3 != k3 or k4 != k4):
	# 				current_curve_array = np.asarray(current_curve_array)
	# 				curves_grid.append(current_curve_array)
	# 				plt.plot(np.asarray(current_curve_x), np.asarray(current_curve_y), color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])
	# 				break

	# 		current_x = current_x + (h/6.)*(k1 + 2*k2 + 2*k3 + k4)
	# 		current_y = current_y + h
	# 		current_theta = np.arctan2(current_x, current_y)
	# 		current_curve_array.append([current_x, current_y])
	# 		current_curve_x.append(current_x)
	# 		current_curve_y.append(current_y)

	# 	else:
	# 		current_curve_array = np.asarray(current_curve_array)
	# 		curves_grid.append(current_curve_array)
	# 		plt.plot(np.asarray(current_curve_x), np.asarray(current_curve_y), color=[np.random.random_sample(),np.random.random_sample(),np.random.random_sample()])

	# 	# print(i)
	# 	i = i + 1

	# curves_grid = np.asarray(curves_grid)

	# # np.savetxt("curves_grid_0.dat", curves_grid[0])
	# # np.savetxt("curves_grid_1.dat", curves_grid[1])
	# # np.savetxt("curves_grid_2.dat", curves_grid[2])


	# plt.savefig("epsilon_streams.pdf")
	# plt.close()

	curves_grid = []

	for current_x in x_eq_array:
		current_curve = []
		for y_i in range(N_y):
			current_y = y_i*y_precision
			current_curve.append([current_x, current_y])
		curves_grid.append(current_curve)

	curves_grid = np.asarray(curves_grid)

	total_data_array = []

	j = 0
	i = 0

	print(len(curves_grid))

	radial_epsilon_spline = interp1d(r_array, epsilon_array, kind='cubic')
	k0=1.0
	boundary_array = []

	for curve in curves_grid:

		# print(curve)
		if j == 0:
			pass
		else:
			boundary_array.append([current_x, current_y])

		k_initial = k0
		current_k = k_initial
		current_epsilon = radial_epsilon_spline(curve[0][0])
		curve_epsilon_array = []
		# print(j)
		j = j + 1
		i = 0

		current_curve_x = curve[:,0]
		current_curve_y = curve[:,1]
		current_x = current_curve_x[0]
		current_y = current_curve_y[0]

		for current_y in current_curve_y:
			curve_epsilon_array.append(current_epsilon)
			total_data_array.append([current_x, current_y, current_epsilon])

			if no_k_flag == False:
				if current_y > boundary_func(current_x):
					break
			if no_k_flag == True:
				pass
			else:
				current_k = k_func(current_x, current_y)

			if (current_y == current_curve_y[-1]):
				break
			if (current_epsilon <= 1e-10):
				break
			else:
				try:
					k1_2 = Cepsilon(M, a, current_x, current_y, l, s0, kappa, gamma, current_k, current_epsilon)
					# k1_2 = Ck(M, a, current_x, current_y, l, s0, kappa, gamma, current_k)
				except ValueError:
					print("Value error 1")
					break
				current_y_2 = current_y + h/2.
				if no_k_flag == True:
					current_k_2 = current_k
				else:
					current_k_2 = k_func(current_x, current_y_2)
				try:
					k2_2 = Cepsilon(M, a, current_x, current_y_2, l, s0, kappa, gamma, current_k_2, current_epsilon + h*k1_2/2.)
					# k2_2 = Ck(M, a, current_x_2, current_y_2, l, s0, kappa, gamma, current_k + h*k1_2/2.)
				except ValueError:
					print("Value error 2")
					break
				try:
					k3_2 = Cepsilon(M, a, current_x, current_y_2, l, s0, kappa, gamma, current_k_2, current_epsilon + h*k2_2/2.)
					# k3_2 = Ck(M, a, current_x_2, current_y_2, l, s0, kappa, gamma, current_k + h*k2_2/2.)
				except ValueError:
					print("Value error 3")
					break
				current_y_4 = current_y + h
				if no_k_flag == True:
					current_k_4 = current_k
				else:
					current_k_4 = k_func(current_x, current_y_4)
				try:
					k4_2 = Cepsilon(M, a, current_x, current_y_4, l, s0, kappa, gamma, current_k_4, current_epsilon + h*k3_2/2.)
					# k4_2 = Ck(M, a, current_x_4, current_y_4, l, s0, kappa, gamma, current_k + h*k3_2)
				except ValueError:
					print("Value error 4")
					break
				if (k1_2 != k1_2 or k2_2 != k2_2 or k3_2 != k3_2 or k4_2 != k4_2):
					print("Value error extra")
					break

				if ((current_epsilon + (h/6.)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)) < 0):
					# print("negative epsilon")
					break
				# elif ((h/6.)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2) > 0):
				# 	print("positive contribution")
				# 	break
				# elif ((current_epsilon > current_epsilon + (h/6.)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2))):
				# 	current_x = current_x_spline(current_y)
				# 	break

				else:
					current_epsilon = current_epsilon + (h/6.)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2)

				i = i + 1
			
	boundary_array.append([current_x, current_y])

	total_data_array = np.asarray(total_data_array)

	# if no_k_flag == True:
	# 	np.savetxt("output_epsilon_no_k.dat", total_data_array)
	# else:
	# 	np.savetxt("output_epsilon.dat", total_data_array)

	# # r_data_array = np.sqrt(total_data_array[:,0]**2+total_data_array[:,1]**2)
	# # theta_data_array = np.arctan2(total_data_array[:,0], total_data_array[:,1])

	# # print(np.shape(r_data_array), np.shape(theta_data_array))

	# # polar_total_data_grid = np.array([r_data_array, theta_data_array]).T

	# full_x_coords = np.concatenate((total_data_array[:,0],total_data_array[:,0]), axis=0)
	# full_y_coords = np.concatenate((total_data_array[:,1],-total_data_array[:,1]), axis=0)
	# cart_total_data_grid = np.array([full_x_coords, full_y_coords]).T

	# # print(polar_total_data_grid)

	# # # print(np.shape(polar_total_data_grid))

	# # k_neighbors_tree = KDTree(polar_total_data_grid)

	# # points = cart_to_polar(np.array([9.8, 1.807]))

	# # print(k_neighbors_tree.query(np.array(points), k=4))

	# x_min_int = np.min(full_x_coords)
	# x_max_int = np.max(full_x_coords)

	# y_min_int = np.min(full_y_coords)
	# y_max_int = np.max(full_y_coords)

	# x_step = 2*int(((x_max_int-x_min_int)/(0.01)))*1j + 1j
	# y_step = 2*int(((y_max_int-y_min_int)/(0.01)))*1j + 1j

	# # print(r_min_int, r_max_int, theta_min_int, theta_max_int)
	# # print(r_data_array)
	# # print(theta_data_array)
	# # print(polar_total_data_grid)
	# # print(total_data_array[:,2])

	# grid_x, grid_y = np.mgrid[x_min_int:x_max_int:x_step, y_min_int:y_max_int:y_step]
	# # print(grid_r)
	# # print(grid_theta)

	# data_to_inter = np.concatenate((total_data_array[:,2], total_data_array[:,2]), axis=0)
	# # data_to_inter -= 1.0
	# grid_k = griddata(cart_total_data_grid, data_to_inter, (grid_x, grid_y), method='cubic', fill_value=0.0)
	# # print(grid_k)


	# xArray = full_x_coords
	# yArray = full_y_coords

	# # xi = np.linspace(np.min(xArray), np.max(xArray),5000)
	# # yi = np.linspace(np.min(yArray), np.max(yArray),5000)

	# DEM = griddata((xArray, yArray),data_to_inter, (grid_x, grid_y), method='cubic')

	# f = plt.figure(figsize=(13,11.5))

	# im = plt.imshow(DEM.transpose(), cmap='plasma', vmin=0.0, vmax=np.max(data_to_inter),
	#         extent=[np.min(grid_x), np.max(grid_x), np.min(grid_y), np.max(grid_y)],
	#         interpolation='nearest', origin='lower')

	# # f.colorbar(im)


	# # plt.axis([0, x_max, -y_max, y_max])
	# # plt.axis([-1000, 1010, -1005, 1005])
	# # plb.ylim([0,r_out/2])
	# # plb.xlim([0,r_out*1.2])
	# plb.ylim([-7.5,7.5])
	# plb.xlim([0,15])

	# plt.tick_params(labelsize=35, pad = 12)
	# lx = plt.xlabel("$r$ sin $\\theta$", fontsize=40)
	# ly = plt.ylabel("$r$ cos $\\theta$", fontsize=40)
	# if (no_k_flag == True):
	# 	plt.savefig("epsilon_no_k.pdf")
	# else:
	# 	plt.savefig("epsilon.pdf")

	# plt.close()
	# # r_1d_array = np.linspace(r_min_int, r_max_int, num=int(np.imag(r_step)))
	# # theta_1d_array = np.linspace(theta_min_int, theta_max_int, num=int(np.imag(theta_step)))

	# # # k_interpolated = RectBivariateSpline(r_1d_array, theta_1d_array, grid_k, kx=1, ky=1)
	# # k_interpolated = RectBivariateSpline(r_1d_array, theta_1d_array, grid_k, kx=3, ky=3)

	# # def k_interpolated_cart(x, y):
	# # 	polar_coords = cart_to_polar([x, y])
	# 	return k_interpolated.ev(polar_coords[0], polar_coords[1])

	# print(k_interpolated_cart(9.8, 1.806))
	# print(k_interpolated_cart(9.769126992296294887, 1.805000000099983515))
	# print(k_interpolated.ev(10, np.pi/2.1))
	# print(k_interpolated.ev(np.sqrt(15**2+6**2), np.arctan2(15, 6)))


	# k_interpolated = SmoothBivariateSpline(r_data_array, theta_data_array, total_data_array[:,2], kx=1, ky=1)




	# print(k_interpolated_cart(9.769126992296294887, 1.805000000099983515))

	# 9.769126992296294887 e + 00 1.805000000099983515 e +00 1.033366028650707147 e + 00

	# y_search_array = np.linspace(0, y_max, num=30*N_y)

	# print(len(r_data_array), len(theta_data_array), len(total_data_array[:,2]))
	# print(r_data_array)
	# print(theta_data_array)
	# print(total_data_array[:,2])




	# # 9.368042994858670980e+00 1.775000000099973274e+00 9.635292310449019570e-01


	# print(len(total_data_array[:,0]), len(total_data_array[:,1]), len(total_data_array[:,2]), len(x_eq_array), len(y_search_array))



	# print(k_interpolated.ev(8.731705734213683456, 6.245000000099889093))
	#8.731705734213683456e+00 6.245000000099889093e+00 1.051271513145822789e-01
	# k_interpolated = RBFInterpolator(total_data_array[:,:2],total_data_array[:,2])
	# print(k_interpolated(np.asarray([8.731705734213683456, 6.245000000099889093])))

	# k_interpolated = SmoothBivariateSpline(total_data_array[:,0], total_data_array[:,1],total_data_array[:,2], s=1)
	# print(k_interpolated.ev(5.882688826865476450, 0.8400000001000006433))
	# k_interpolated = bisplrep(total_data_array[:,0], total_data_array[:,1],total_data_array[:,2])

	# print(bisplev(5.372034947596807442, 3.955000000099937907, k_interpolated))
	# print(bisplev(5.882688826865476450, 0.8400000001000006433, k_interpolated))
	# # 2.194168934260585235e-01
	# 9.767229419719847261e-01
	# print(k_interpolated.ev(5.882688826865476450, 0.8400000001000006433))

	boundary_array = np.array(boundary_array)
	print(np.max(total_data_array[:,2]))
	# print(np.max(yArray))
	if no_k_flag == True:
		return boundary_array
	else:
		pass



