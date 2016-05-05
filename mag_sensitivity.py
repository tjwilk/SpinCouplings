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
grav_acceleration = 9.8 # in m/s^2

num_shots = pow(10, 6) #total integration time of experiment

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

	return noise

def background(sens, t, m, Q, relaxation_rate):
	coherence = Q/m
	signal_time = min(coherence, pow(relaxation_rate, -1.0))

	background = 0
	if t <= signal_time:
		background = sens*pow(t, -0.5)
	elif t > signal_time:
		background = sens*pow(signal_time, -0.5)*pow(signal_time/t, 0.25)
	return background

# def coupling_form(n, p, mu, rho, v):
def axion_coupling_form(noise, rho, v, mass, Q, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, coupling_bool):

	coherence = Q/mass
	signal_time = min(coherence, pow(relaxation_rate, -1.0))
	precession_freq = 10 #something like this

	num_periods = floor(mass*(signal_time)/(math.pi))
	remainder_time = signal_time - num_periods*math.pi/(mass)
	oscillation = abs(math.sin(mass*remainder_time))
	# print oscillation

	grad_axion = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation # in eV^2

	# signal = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation # in eV^2
	signal = 0
	signal_e = grad_axion/el_mag_moment # in eV*T
	signal_nucl = grad_axion/nucl_mag_moment # in eV*T

	if coupling_bool == 0: signal = signal_e
	if coupling_bool == 1: signal = signal_nucl

	coupling = pow(10, 9)*noise/signal # in GeV (since noise is in T)

	return coupling


# def coupling_form(n, p, mu, rho, v):
def axion_coupling_form_polarization(noise, rho, v, mass, Q, z_pol, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, coupling_bool):

	coherence = Q/mass
	signal_time = min(coherence, pow(relaxation_rate, -1.0))
	precession_freq = 10 #something like this (in Hz)

	num_periods = floor(mass*(signal_time)/(math.pi))
	remainder_time = signal_time - num_periods*math.pi/(mass)
	# oscillation = abs(math.sin(mass*remainder_time))
	oscillation = abs(math.sin(mass*signal_time))

	grad_axion = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation # in eV^2

	# signal_lowfreq = pow(hbar, -1.0)*z_pol*grad_axion*amplitude*nucl_relax_time/(relaxation_rate*nucl_relax_time + precession_freq*relaxation_rate*0.02 + nucl_relax_time*10*z_pol)
	# signal_highfreq = z_pol*grad_axion*pow(hbar*relaxation_rate, -1.0)*suppression*amplitude

	# phi = 1/2 * l * r_e * c * f * n_a * D(nu)
	length = 2.4 # cm
	r_e = 2.8*pow(10, -13) # cm
	c = 3*pow(10, 10) # cm/s
	n = 7*pow(10, 13) # cm^(-3)
	lambda_0 = 770.1*pow(10, -7) # in cm
	lambda_p = 770.76*pow(10, -7) # in cm
	nu_0 = c/lambda_0
	nu_p = c/lambda_p
	delta_nu = 13.2*pow(10, 10) # in Hz

	D_nu = (nu_0 - nu_p)/(pow((nu_0 - nu_p), 2) + pow(delta_nu*0.5, 2))

	pol_to_angle = 0.5*length*r_e*c*n*D_nu

	# print D_nu

	signal = 0
	signal_e = pol_to_angle*z_pol*grad_axion/(mass*hbar) # in eV (1/mass term comes from integrating the frequency formula)
	signal_nucl = signal_e*(el_mag_moment/nucl_mag_moment)

	if coupling_bool == 0: signal = signal_e
	if coupling_bool == 1: signal = signal_nucl

	coupling = pow(10, 9)*noise/signal # noise should be dimensionless
	# if mass < precession_freq: coupling = pow(10, 9)*noise/signal # noise should be dimensionless

	# if mass < pow(10, -2): coupling = pow(10, 9)*noise/signal_lowfreq # in GeV (since noise is in T)
	# elif mass > pow(10, -2) and mass < precession_freq: coupling = pow(10, 9)*noise/signal_highfreq

	return coupling

def axion_coupling_form_germans(noise, rho, v, mass, Q, relaxation_rate, shot_time, nucl_mag_moment, german_frequency):

	coherence = Q/mass
	signal_time = min(coherence, pow(relaxation_rate, -1.0))
	precession_freq = 0.001 #something like this 
	earth_freq = pow(10, -6) #something like this 


	sideband_1 = german_frequency + mass
	sideband_2 = abs(german_frequency - mass)

	num_periods = floor(mass*(signal_time)/(math.pi))
	remainder_time = signal_time - num_periods*math.pi/(mass)
	oscillation = earth_freq*abs(math.sin(mass*remainder_time))/(mass)
	# oscillation = precession_freq/(mass)
	# print oscillation

	grad_axion = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation # in eV^2
	# signal = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5)*oscillation # in eV^2

	signal_nucl = grad_axion/nucl_mag_moment # in eV*T

	coupling = pow(10, 9)*noise/signal_nucl # in GeV (since noise is in T)

	return coupling

def axion_coupling_form_germans_sample(noise, rho, v, mass, Q, relaxation_rate, shot_time, transverse_mag_amplitude):

	coherence = Q/mass
	signal_time = min(coherence, pow(relaxation_rate, -1.0), shot_time)
	precession_freq = 10 #something like this 


	sideband_1 = precession_freq + mass
	sideband_2 = abs(precession_freq - mass)

	oscillation_1 = abs(math.sin(sideband_1*shot_time))
	oscillation_2 = abs(math.sin(sideband_2*shot_time))

	# print oscillation_1, oscillation_2

	grad_axion = v*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5) # in eV^2
	signal = transverse_mag_amplitude*grad_axion*signal_time*max(oscillation_1, oscillation_2)/hbar # in eV*T

	# freq_offset = pow(10, -9)*grad_axion*math.sin(mass*shot_time)
	# offset = math.sin(precession_freq*shot_time)- math.sin((precession_freq + freq_offset)*shot_time)
	# print freq_offset, offset
	# if abs(offset) > pow(10, -15): print mass


	coupling = pow(10, 9)*noise/signal

	return coupling

axion_masses = np.linspace(-22.0, -10.0, 600)
axion_result_1 = [0]*len(axion_masses)
axion_result_2 = [0]*len(axion_masses)
axion_result_3 = [0]*len(axion_masses)
axion_result_4 = [0]*len(axion_masses)
axion_result_germans = [0]*len(axion_masses)
axion_result_germans_2 = [0]*len(axion_masses)
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
earth_frequency = 11*pow(10, -6) # sidereal frequency in Hz
# german_frequency = pow(10, -3) # sidereal frequency in Hz
german_frequency = 10 #something like this
german_factor = 1 - 13.0/4.7 # 1 - gamma_He/gamma_Xe
german_shot_run = 86400 #they run for ~a day
german_shot_samples = 3 #num seconds per sample
num_german_shots = 100
german_relaxation_rate = 3*pow(10, -5) # in 1/s

transverse_mag_amplitude = 10*pow(10, -9) #10 pT is from paper


nucl_mag_moment = 2.0*math.pi*3.243*pow(10, 3)*hbar*pow(10, 4) # from Hz/G -> eV/T
el_mag_moment = 2.0*math.pi*2.8*pow(10, 6)*hbar*pow(10, 4) # in eV/T

root_Pmag_1 = pow(10, -15) # T/sqrt(Hz)
root_Pmag_2 = pow(10, -17) # T/sqrt(Hz)

# root_Pmag_1 = pow(10, -15)*abs(nucl_mag_moment - el_mag_moment) # T/sqrt(Hz)*eV/T = eV/sqrt(Hz)
# root_Pmag_2 = pow(10, -17)*abs(nucl_mag_moment - el_mag_moment) # T/sqrt(Hz)*eV/T = eV/sqrt(Hz)

root_Ptheta_1 = 4*pow(10, -8) # rad/sqrt(Hz)
root_Ptheta_2 = pow(10, -9) # rad/sqrt(Hz)

vel = pow(10, -3) # in c = 1
rho = 0.3 # GeV/cm^3
Q = pow(10, 6)

# experimental parameters

K_density = 7*pow(10, 13) # in 1/cm^3
volume = 4*pow(2.4, 3) # in cm^3
num_K = K_density*volume
ez_pol = 0.5
relaxation_rate = 353 # in 1/s
nuclear_relax_time = 2*pow(10, -4) # in 1/s
shot_time = 5 # number of seconds they record

# I believe this to be the shot-noise limited polarization sensitivity
test_sensitivity = pow(num_K*relaxation_rate, -0.5)

test_mag_sensitivity_1 = root_Ptheta_1*relaxation_rate*hbar*pow(el_mag_moment, -1.0)
test_mag_sensitivity_2 = root_Ptheta_2*relaxation_rate*hbar*pow(el_mag_moment, -1.0)

# root_Ptheta_1 = test_sensitivity
# print test_mag_sensitivity_1, test_mag_sensitivity_2

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

	sample_resolved = false

	mass = pow(10, axion_masses[i])/hbar #put it in terms of Hz
	sideband_1 = earth_frequency + mass
	sideband_2 = abs(earth_frequency - mass)

	bandwidth = math.pi/german_shot_samples
	DM_bandwidth = min(mass/Q, bandwidth)
	# german_freq_resolution = pow(bandwidth*pow(german_shot_samples, 3), -0.5)

	german_noise_sb1 = german_noise(german_noise_floor, sideband_1)
	german_noise_sb2 = german_noise(german_noise_floor, sideband_2)

	# print(sideband_1, sideband_2)

	noise_1_germans_shot = background(german_noise_sb1, german_shot_run, mass, Q, german_relaxation_rate)
	noise_1_germans = pow(pow(noise_1_germans_shot, 2)/num_german_shots, 0.5)

	noise_2_germans_shot_sb1 = background(german_noise_sb1, german_shot_samples, mass, Q, german_relaxation_rate)
	noise_2_germans_shot_sb2 = background(german_noise_sb1, german_shot_samples, mass, Q, german_relaxation_rate)
	noise_2_germans = pow(pow(min(noise_2_germans_shot_sb1, noise_2_germans_shot_sb2), 2)/(num_german_shots*german_shot_run/german_shot_samples), 0.5)

	if bandwidth/2.0 < (mass + DM_bandwidth/2.0) : 
		sample_resolved = true
		# print "yay ", mass*hbar

	axion_result_germans[i] = math.log10(axion_coupling_form_germans(noise_1_germans, rho, vel, mass, Q, german_relaxation_rate, german_shot_run, nucl_mag_moment, german_frequency))
	if sample_resolved: 
		axion_result_germans_2[i] = math.log10(axion_coupling_form_germans_sample(noise_2_germans, rho, vel, mass, Q, german_relaxation_rate, german_shot_samples, transverse_mag_amplitude))
	else: axion_result_germans_2[i] = 1.0

	# total_time = integration_time
	total_time = shot_time

	#run 1
	noise_1_shot = background(root_Pmag_1, shot_time, mass, Q, relaxation_rate)
	noise_2_shot = background(root_Pmag_2, shot_time, mass, Q, relaxation_rate)

	noise_1 = pow(pow(noise_1_shot, 2)/num_shots, 0.5)
	noise_2 = pow(pow(noise_2_shot, 2)/num_shots, 0.5)

	noise_pol_1_shot = background(root_Ptheta_1, shot_time, mass, Q, relaxation_rate)
	noise_pol_2_shot = background(root_Ptheta_2, shot_time, mass, Q, relaxation_rate)

	noise_pol_1 = pow(pow(noise_pol_1_shot, 2)/num_shots, 0.5)
	noise_pol_2 = pow(pow(noise_pol_2_shot, 2)/num_shots, 0.5)

	# axion_result_1[i] = math.log10(axion_coupling_form(noise_1, rho, vel, mass, Q, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 0))
	# axion_result_2[i] = math.log10(axion_coupling_form(noise_2, rho, vel, mass, Q, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 0))

	# axion_result_3[i] = math.log10(axion_coupling_form(noise_1, rho, vel, mass, Q, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 1))
	# axion_result_4[i] = math.log10(axion_coupling_form(noise_2, rho, vel, mass, Q, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 1))

	axion_result_1[i] = math.log10(axion_coupling_form_polarization(noise_pol_1, rho, vel, mass, Q, ez_pol, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 0))
	axion_result_2[i] = math.log10(axion_coupling_form_polarization(noise_pol_2, rho, vel, mass, Q, ez_pol, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 0))

	axion_result_3[i] = math.log10(axion_coupling_form_polarization(noise_pol_1, rho, vel, mass, Q, ez_pol, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 1))
	axion_result_4[i] = math.log10(axion_coupling_form_polarization(noise_pol_2, rho, vel, mass, Q, ez_pol, relaxation_rate, shot_time, el_mag_moment, nucl_mag_moment, 1))


	astrophysics_bounds[i] = math.log10(astrophysics_bound)


fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

ax.plot(axion_masses, axion_result_1)
ax.plot(axion_masses, axion_result_2)
ax.plot(axion_masses, astrophysics_bounds)
ax.legend(("Current", "Shot-noise Limited"), loc=2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('$g_{aee}$ (GeV)$^{-1}$ (log scale)');
ax.set_title('axion sensitivity (electron coupling)')
# ax.legend(("G", "E", "CalcG", "CalcE"));
plt.savefig("magnetometer_sensitivity_el.pdf")
plt.show()

fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

ax.plot(axion_masses, axion_result_3)
ax.plot(axion_masses, axion_result_4)
ax.plot(axion_masses, axion_result_germans)
ax.plot(axion_masses, axion_result_germans_2)
ax.plot(axion_masses, astrophysics_bounds)
# ax.legend(("Current", "Shot-noise Limited", "Germans", "Germans 2"), loc=2);
ax.legend(("Current", "Shot-noise Limited", "Astrophysics"), loc=2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('$g_{aNN}$ (GeV)$^{-1}$ (log scale)');
ax.set_title('axion sensitivity (nuclear coupling)')
# ax.legend(("G", "E", "CalcG", "CalcE"));
plt.savefig("magnetometer_sensitivity_nucl.pdf")
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