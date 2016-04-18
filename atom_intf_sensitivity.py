import math
import random
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.constants import *
from sympy import *

# conversions into eV
boltzmann = 8.6*pow(10, -5) #in eV/K
hbar = 6.6*pow(10, -16) # in eV-s
hbarc = 197.3*pow(10, -9) # in eV-m
evJ_conv = 6.24*pow(10, 18) # in eV/J
c = 3.0*pow(10, 8) # in m/s
num_shots = pow(10, 7) #total integration time of experiment
nucleon_mass = 8.16*pow(10, 10) # in eV, for strontium
grav_acceleration = 9.8 # in m/s^2

def integration_time(free_fall_time, m, Q):
	# bin_range = pow(10, 3)
	bin_range = pow(10, 1)
	bin_width = math.pi/free_fall_time
	num_bins = bin_range/bin_width
	coherence_time = Q/m

	time = num_shots
	if free_fall_time > math.pi/m: 
		Qeff = m*free_fall_time/math.pi
		norm = num_bins*Qeff
		time = num_shots/norm
	return time


def background(sens, t, m, Q):
	coherence = Q/m
	background = 0
	if t <= coherence:
		background = sens*pow(t, -0.5)
	elif t > coherence:
		background = sens*pow(coherence, -0.5)*pow(coherence/t, 0.25)
	return background

# def coupling_form(n, p, mu, rho, v):
def axion_coupling_form(noise, spin, rho, v, mass, free_fall_time, resvbroad):

	num_periods = floor(mass*(free_fall_time)/(math.pi))
	remainder_time = free_fall_time - num_periods*math.pi/(mass)
	squeezing = pow(10, 0)

	oscillation = 0
	if resvbroad == 1:
		# oscillation = (1 - math.cos(mass*remainder_time))/mass + num_periods
		oscillation = math.sin(mass*remainder_time) + num_periods
	elif resvbroad == 0:	
		# oscillation = (1 - math.cos(mass*remainder_time))/mass
		oscillation = math.sin(mass*remainder_time)
	# coupling = pow(squeezing, -1)*pow(10, -9)*noise*mass/(spin*v*pow(rho, 0.5)*oscillation)*1.0/pow(hbar*pow(c, 3), 0.5)*pow(10, 15.0/2) # in inverse GeV
	signal = spin*v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation*1.0/(mass*hbar)
	coupling = pow(10, 9)*noise/signal
	return coupling

# def axion_force(noise, keff, free_fall_time, mass, distance):
# 	new_dist = distance/c
# 	# new_dist = keff*hbar/mass_nucleon*min(free_fall_time, math.pi/mass)
# 	dist_dep = pow(new_dist, -4)*(pow(mass*new_dist, 3) + 3*pow(mass*new_dist, 2) + 6*pow(mass*new_dist, 3) + 6)
# 	force = pow(hbar, 2)*math.exp(-mass*new_dist)*dist_dep/(4.0*math.pi) #  in ev^2/s^2
# 	acceleration = force/nucleon_mass*hbarc # in eV/s^2*eV*m = eV^2*m/s^2
# 	coupling = pow(10, -9)*pow(noise/acceleration, 0.5) #in sqrt(1/eV^2)*10^(-9) = 1/GeV
# 	# print acceleration*hbarc, noise
# 	return coupling

# def axion_phase(noise, keff, free_fall_time, mass, distance):
# 	new_dist = distance/c
# 	# nucleon_mass = 8.16*pow(10, 10) # in eV
# 	dist_dep = pow(new_dist, -4)*(pow(mass*new_dist, 3) + 3*pow(mass*new_dist, 2) + 6*pow(mass*new_dist, 3) + 6) # in s^-4
# 	force = pow(hbar, 2)*math.exp(-mass*new_dist)*dist_dep/(4.0*math.pi) #  in ev^2/s^2
# 	acceleration = force/nucleon_mass # in eV/s^2
# 	phase = hbar*keff*acceleration*pow(min(free_fall_time, math.pi/mass),2) #eV*s*Hz*eV/s^2*s^2 = eV^2
# 	print phase
# 	# phase = hbar*keff*acceleration*pow(free_fall_time,2) #eV*s*Hz*eV/s^2*s^2 = eV^2
# 	coupling = pow(10, -9)*pow(noise/phase, 0.5)
# 	# return phase
# 	return coupling

# def axion_potential_phase(noise, keff, free_fall_time, mass):
# 	particle_velocity = keff/nucleon_mass #this should be dimensionless
# 	# final_x = mass*particle_velocity*min(free_fall_time, math.pi/mass) # should also be dimensionless
# 	final_x = mass*particle_velocity*free_fall_time # should also be dimensionless
# 	initial_x = mass*pow(10, -8)
# 	phase = pow(hbar*mass, 2)/(4*math.pi*particle_velocity)*((pow(final_x, 2)/8.0 - pow(final_x, -2) - final_x/3.0) - ((pow(initial_x, 2)/8.0 - pow(initial_x, -2) - initial_x/3.0)))
# 	# phase = pow(hbar*mass, 2)/(4*math.pi*particle_velocity)*(- pow(final_x, -2) + pow(initial_x, -2))
# 	coupling = pow(10, -9)*pow(noise/abs(phase), 0.5)
# 	# print particle_velocity
# 	return coupling

def axion_force(noise, spin, rho, v, mass, keff, free_fall_time):
	acceleration = spin*pow(v, 2)*mass/(mass_nucleon/hbar)*pow(2*rho*pow(10, 15)*pow(hbarc, 3), 0.5)*c/hbar #in eV^2*m/s*1/(eV*s) = eV*m/s^2
	coupling = pow(10, 9)*noise/acceleration
	return coupling

def axion_force_phase(noise, spin, rho, v, mass, keff, free_fall_time):
	acceleration = spin*pow(v, 2)*mass/(mass_nucleon/hbar)*pow(2*rho*pow(10, 15)*pow(hbarc, 3), 0.5) #in eV^2
	# phase = keff*acceleration*pow(min(free_fall_time, math.pi/mass), 2)*pow(hbar, -1.0) #Hz*eV^2*s^2*1/(eV*s) = eV
	phase = keff*acceleration*pow(free_fall_time, 2)*pow(hbar, -1) #Hz*eV^2*s^2*1/(eV*s) = eV
	coupling = pow(10, 9)*noise/phase
	return coupling

# def axial_vector_coupling_form(mass, noise, Np, rho):
# 	coupling = 2.0/(Np)*(mass/pow(10, 6))*pow(0.3/(2*rho), 0.5)*noise*1/(2.31*pow(10, 6)) #dimensionless
# 	return coupling


axion_masses = np.linspace(-22.0, -10.0, 600)
axion_result_1 = [0]*len(axion_masses)
axion_result_1_broadband = [0]*len(axion_masses)
axion_result_acceleration = [0]*len(axion_masses)
axion_result_2 = [0]*len(axion_masses)
axion_result_2_broadband = [0]*len(axion_masses)
astrophysics_bounds = [0]*len(axion_masses)

vector_masses = np.linspace(-22.0, -10.0, 600)
vector_result_1 = [0]*len(vector_masses)
vector_result_2 = [0]*len(vector_masses)
vector_result_3 = [0]*len(vector_masses)
vector_astro = [0]*len(vector_masses)

# coupling_form(pow(10, 6), 30*24*60*60.0, 1.0, 0.1, 0.3, pow(10, -6))
# integration_time = pow(10, 8)
# integration_time = pow(10, 7)

root_Ptheta_1 = pow(10, -4) # rad/sqrt(Hz)
root_Ptheta_2 = 3*pow(10, -5) # rad/sqrt(Hz)
root_Paccel_1 = pow(10, -15)*grav_acceleration # m/s^2/sqrt(Hz)

# Np = pow(10, 6) # will this make an actual difference? 
Np = pow(10, 0) # will this make an actual difference? 
spin = 9.0
# spin = 1.0
vel = pow(10, -3) # in c = 1
# rho = 3*pow(10, 14) # eV/m^3
rho = 0.3 # GeV/cm^3
free_fall_time = 1.15 # in seconds
free_fall_time_2 = 100
# free_fall_time = 2.3 # in seconds (twice for the two arms)

keff = 4.35*pow(10, 14) # in Hz

Q = pow(10, 6)
test_mass = pow(10, 6) # in Hz

alpha = 1.0/137
mass_nucleon = 500*pow(10, 6) # in eV
temp = 30*pow(10, 6) # in eV
energy_loss_rate = 2.7*temp*1.8*pow(10, 38)*0.76*pow(temp, 3)/(3*pow(10, 14))
loss_rate_bound = pow(10, 19)*pow(10, -7)*evJ_conv
# astro_coupling = 6.8*pow(10, -11)*(mass_nucleon/temp)
vector_astro_coupling = pow(loss_rate_bound/energy_loss_rate, 0.5)*(mass_nucleon/temp)*pow(alpha, -0.5)
axion_astro_coupling = pow(loss_rate_bound/energy_loss_rate, 0.5)*(mass_nucleon/temp)*pow(alpha, -0.5)

astrophysics_bound = axion_astro_coupling
vector_astro_bound = vector_astro_coupling

for i in range(0, len(axion_masses)):

	mass = pow(10, axion_masses[i])/hbar #put it in terms of Hz
	total_time = integration_time(free_fall_time, mass, Q)

	#run 1
	noise_1_shot = background(root_Ptheta_1, free_fall_time, mass, Q)
	noise_1 = pow(pow(noise_1_shot, 2)/num_shots, 0.5)


	#run 1
	noise_1 = background(root_Ptheta_1, total_time, mass, Q)
	noise_2 = background(root_Ptheta_2, total_time, mass, Q)
	axion_result_1[i] = math.log10(axion_coupling_form(noise_1, Np*spin, rho, vel, mass, free_fall_time, 1))
	axion_result_2[i] = math.log10(axion_coupling_form(noise_2, Np*spin, rho, vel, mass, free_fall_time_2, 1))

	# acceleration_noise_a = background(root_Paccel_1, num_shots, mass, Q)
	# astrophysics_bounds[i] = math.log10(axion_force(acceleration_noise_a, Np*spin, rho, vel, mass, free_fall_time, 0.5))

	# acceleration_noise_b = background(root_Ptheta_1, num_shots, mass, Q)
	# axion_result_acceleration[i] = math.log10(axion_force_phase(acceleration_noise_b, Np*spin, rho, vel, mass, free_fall_time, 0.5))

	noise_1_broadband = background(root_Ptheta_1, num_shots, mass, Q)
	noise_2_broadband = background(root_Ptheta_2, num_shots, mass, Q)
	axion_result_1_broadband[i] = math.log10(axion_coupling_form(noise_1_broadband, Np*spin, rho, vel, mass, free_fall_time, 0))

	astrophysics_bounds[i] = math.log10(astrophysics_bound)
	vector_astro[i] = math.log10(vector_astro_bound)

	#run 2
	# noise_2 = max(twist_noise_2, thermal_noise_2, total_systematic_uncertainty)
	# axion_result_2[i] = math.log10(axion_coupling_form(noise_2, Np, rho))

	# #run 3
	# noise_3 = max(twist_noise_3, thermal_noise_3, total_systematic_uncertainty)
	# axion_result_3[i] = math.log10(axion_coupling_form(noise_3, Np*pend_mass_3/pend_mass_1, rho))

fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

ax.plot(axion_masses, axion_result_1)
ax.plot(axion_masses, axion_result_2)
ax.plot(axion_masses, axion_result_1_broadband)
# ax.plot(axion_masses, axion_result_acceleration)
ax.plot(axion_masses, astrophysics_bounds)
# ax.legend(("Resonant", "Broadband", "Astrophysics"), loc=2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('$g_{aNN}$ (GeV)$^{-1}$ (log scale)');
ax.set_title('axion sensitivity')
# ax.legend(("G", "E", "CalcG", "CalcE"));
plt.savefig("atom_intf_sensitivity.pdf")
plt.show()

# # now for the axial vector sensitivity

# for i in xrange(0, len(vector_masses)):

# 	mass = pow(10, vector_masses[i])/hbar

# 	#run 1
# 	twist_noise_1 = background(pow(noise_twist_power(pend_moment_1, mass, freq_1, pend_Q_1, Ptheta_1), 0.5), integration_time, mass, Q)
# 	thermal_noise_1 = background(pow(noise_thermal_power(temp_1, kappa_1, freq_1, pend_Q_1), 0.5), integration_time, mass, Q)
# 	noise_1 = max(twist_noise_1, thermal_noise_1)
# 	vector_result_1[i] = math.log10(axial_vector_coupling_form(mass, noise_1, Np, rho))

# 	#run 2
# 	twist_noise_2 = background(pow(noise_twist_power(pend_moment_2, mass, freq_2, pend_Q_2, Ptheta_2), 0.5), integration_time, mass, Q)
# 	thermal_noise_2 = background(pow(noise_thermal_power(temp_2, kappa_2, freq_2, pend_Q_2), 0.5), integration_time, mass, Q)
# 	noise_2 = max(twist_noise_2, thermal_noise_2)
# 	vector_result_2[i] = math.log10(axial_vector_coupling_form(mass, noise_2, Np, rho))

# 	#run 3
# 	twist_noise_3 = background(pow(noise_twist_power(pend_moment_3, mass, freq_3, pend_Q_3, Ptheta_3), 0.5), integration_time, mass, Q)
# 	thermal_noise_3 = background(pow(noise_thermal_power(temp_3, kappa_3, freq_3, pend_Q_3), 0.5), integration_time, mass, Q)
# 	noise_3 = max(twist_noise_3, thermal_noise_3)
# 	vector_result_3[i] = math.log10(axial_vector_coupling_form(mass, noise_3, Np*pend_mass_3/pend_mass_1, rho))


fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

# ax.plot(vector_masses, vector_result_1)
# ax.plot(vector_masses, vector_result_2)
# ax.plot(vector_masses, vector_result_3)
ax.plot(axion_masses, axion_result_1)
ax.plot(axion_masses, axion_result_1_broadband)
# ax.plot(axion_masses, axion_result_3)
ax.plot(axion_masses, vector_astro)
# ax.legend(("Run 1", "Run 2", "Run 3", "Astro"));
ax.legend(("Resonant", "Broadband","Astro"), loc = 2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('g_{ANN} (log scale)');
ax.set_title('hidden photon sensitivity')
plt.show()

# fig, ax = plt.subplots()
# # handles, labels = ax.get_legend_handles_labels()
# # ax.legend(handles, labels)

# ax.plot(vector_masses, vector_result_1)
# ax.plot(vector_masses, vector_result_2)
# ax.plot(vector_masses, vector_result_3)
# ax.legend(("Run 1", "Run 2", "Run 3"));
# ax.set_xlabel('mass (eV) (log scale)');
# ax.set_ylabel('g_{aee} (log scale)');
# ax.set_title('axial vector sensitivity')
# # ax.legend(("G", "E", "CalcG", "CalcE"));
# plt.show()