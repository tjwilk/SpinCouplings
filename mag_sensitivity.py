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
evT_conv = 1.44*pow(10, -3) # in T/eV^2

c = 3.0*pow(10, 8) # in m/s
num_shots = pow(10, 7) #total integration time of experiment
nucleon_mass = 8.16*pow(10, 10) # in eV, for strontium
grav_acceleration = 9.8 # in m/s^2


integration_time = num_shots

# def integration_time(free_fall_time, m, Q):
# 	# bin_range = pow(10, 3)
# 	bin_range = pow(10, 1)
# 	bin_width = math.pi/free_fall_time
# 	num_bins = bin_range/bin_width
# 	coherence_time = Q/m

# 	time = num_shots
# 	if free_fall_time > math.pi/m: 
# 		Qeff = m*free_fall_time/math.pi
# 		norm = num_bins*Qeff
# 		time = num_shots/norm
# 	return time

#returns a noise in fT/sqrt(Hz)
def german_noise(floor, freq):
	noise = 0
	if freq > 1: 
		noise = floor
	else: noise = floor/(freq)

def background(sens, t, m, Q):
	coherence = Q/m
	background = 0
	if t <= coherence:
		background = sens*pow(t, -0.5)
	elif t > coherence:
		background = sens*pow(coherence, -0.5)*pow(coherence/t, 0.25)
	return background

# def coupling_form(n, p, mu, rho, v):
def axion_coupling_form(noise, rho, v, mass, Q, mag_moment, relaxation_time, shot_time):

	coherence = Q/mass
	signal_time = min(coherence, pow(relaxation_time, -1.0))
	precession_freq = 10 #something like this
	sideband_freq1 = precession_freq + mass
	sideband_freq2 = abs(precession_freq - mass)

	num_periods = floor(mass*(signal_time)/(math.pi))
	remainder_time = signal_time - num_periods*math.pi/(mass)
	oscillation = abs(math.sin(2.0*math.pi*mass*remainder_time))
	# print oscillation


	signal = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation # in eV^2
	# if mass < precession_freq: print oscillation

	print "blah blah test blah"

	coupling = 1.0
	if mass < precession_freq: coupling = pow(10, 9)*noise/signal # in GeV (since noise is in eV)

	return coupling


# def coupling_form(n, p, mu, rho, v):
def axion_coupling_form_polarization(noise, rho, v, mass, Q, z_pol, relaxation_time, el_mag_moment, nucl_mag_moment, shot_time):

	coherence = Q/mass
	signal_time = min(coherence, pow(relaxation_time, -1.0))
	precession_freq = 10 #something like this (in Hz)
	sideband_freq1 = precession_freq + mass
	sideband_freq2 = abs(precession_freq - mass)

	num_periods = floor(sideband_freq1*(signal_time)/(math.pi))
	remainder_time = signal_time - num_periods*math.pi/(sideband_freq1)
	oscillation = abs(math.sin(2.0*math.pi*sideband_freq1*remainder_time))
	# suppression = 1.0

	grad_axion = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation # in eV^2

	# signal_lowfreq = pow(hbar, -1.0)*z_pol*grad_axion*amplitude*nucl_relax_time/(relaxation_time*nucl_relax_time + precession_freq*relaxation_time*0.02 + nucl_relax_time*10*z_pol)
	# signal_highfreq = z_pol*grad_axion*pow(hbar*relaxation_time, -1.0)*suppression*amplitude

	signal_e = z_pol*grad_axion/(relaxation_time*hbar) # in eV
	signal_nucl = signal_e*(el_mag_moment/nucl_mag_moment)

	coupling = 1.0
	if mass < precession_freq: coupling = pow(10, 9)*noise/signal_nucl # noise should be dimensionless

	# if mass < pow(10, -2): coupling = pow(10, 9)*noise/signal_lowfreq # in GeV (since noise is in T)
	# elif mass > pow(10, -2) and mass < precession_freq: coupling = pow(10, 9)*noise/signal_highfreq

	return coupling

axion_masses = np.linspace(-22.0, -10.0, 600)
axion_result_1 = [0]*len(axion_masses)
axion_result_2 = [0]*len(axion_masses)
axion_result_3 = [0]*len(axion_masses)
astrophysics_bounds = [0]*len(axion_masses)

vector_masses = np.linspace(-22.0, -10.0, 600)
vector_result_1 = [0]*len(vector_masses)
vector_result_2 = [0]*len(vector_masses)
vector_result_3 = [0]*len(vector_masses)
vector_astro = [0]*len(vector_masses)

# coupling_form(pow(10, 6), 30*24*60*60.0, 1.0, 0.1, 0.3, pow(10, -6))
# integration_time = pow(10, 8)
# integration_time = pow(10, 7)

german_noise_floor = 2.3*pow(10, -15) # in T/sqrt(Hz)


nucl_mag_moment = 2.0*math.pi*3.243*pow(10, 3)*hbar*pow(10, 4) # from Hz/G -> eV/T
el_mag_moment = 2.0*math.pi*2.8*pow(10, 6)*hbar*pow(10, 4) # in eV/T

root_Pmag_1 = pow(10, -15)*el_mag_moment # T/sqrt(Hz)*eV/T = eV/sqrt(Hz)
root_Pmag_2 = pow(10, -17)*el_mag_moment # T/sqrt(Hz)*eV/T = eV/sqrt(Hz)

# root_Pmag_1 = pow(10, -15)*abs(nucl_mag_moment - el_mag_moment) # T/sqrt(Hz)*eV/T = eV/sqrt(Hz)
# root_Pmag_2 = pow(10, -17)*abs(nucl_mag_moment - el_mag_moment) # T/sqrt(Hz)*eV/T = eV/sqrt(Hz)

# where did this come from?
root_Ptheta_1 = 4*pow(10, -8) # rad/sqrt(Hz)

vel = pow(10, -3) # in c = 1
rho = 0.3 # GeV/cm^3

ez_pol = 0.5
relaxation_time = 353 # in 1/s
nuclear_relax_time = 2*pow(10, -4) # in 1/s

# relaxation_time = 1 # 1 second
shot_time = 86400 # 1 day
Q = pow(10, 6)




#astro calculation
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
	total_time = integration_time

	#run 1
	noise_1 = background(root_Pmag_1, total_time, mass, Q)
	noise_2 = background(root_Pmag_2, total_time, mass, Q)

	axion_result_1[i] = math.log10(axion_coupling_form(noise_1, rho, vel, mass, Q, nucl_mag_moment, relaxation_time, shot_time))
	axion_result_2[i] = math.log10(axion_coupling_form(noise_2, rho, vel, mass, Q, nucl_mag_moment, relaxation_time, shot_time))

	noise_pol_1 = background(root_Ptheta_1, total_time, mass, Q)
	axion_result_3[i] = math.log10(axion_coupling_form_polarization(noise_pol_1, rho, vel, mass, Q, ez_pol, relaxation_time, el_mag_moment, nucl_mag_moment, shot_time))


	astrophysics_bounds[i] = math.log10(astrophysics_bound)


fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

ax.plot(axion_masses, axion_result_1)
ax.plot(axion_masses, axion_result_2)
ax.plot(axion_masses, axion_result_3)
ax.plot(axion_masses, astrophysics_bounds)
# ax.legend(("Mag Analysis", "Polarization"), loc=2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('$g_{aee} - g_{aNN}$ (GeV)$^{-1}$ (log scale)');
ax.set_title('axion sensitivity')
# ax.legend(("G", "E", "CalcG", "CalcE"));
plt.savefig("magnetometer_sensitivity.pdf")
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


# fig, ax = plt.subplots()
# # handles, labels = ax.get_legend_handles_labels()
# # ax.legend(handles, labels)

# # ax.plot(vector_masses, vector_result_1)
# # ax.plot(vector_masses, vector_result_2)
# # ax.plot(vector_masses, vector_result_3)
# ax.plot(axion_masses, axion_result_1)
# ax.plot(axion_masses, axion_result_broadband)
# # ax.plot(axion_masses, axion_result_3)
# ax.plot(axion_masses, vector_astro)
# # ax.legend(("Run 1", "Run 2", "Run 3", "Astro"));
# ax.legend(("Resonant", "Broadband","Astro"), loc = 2);
# ax.set_xlabel('mass (eV) (log scale)');
# ax.set_ylabel('g_{ANN} (log scale)');
# ax.set_title('hidden photon sensitivity')
# plt.show()

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