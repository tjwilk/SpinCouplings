import math
import random
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.constants import *
from sympy import *
from qutip import *

# conversions into eV
boltzmann = 8.6*pow(10, -5) #in eV/K
hbar = 6.6*pow(10, -16) # in eV-s
hbarc = 197.3*pow(10, -9) # in eV-m
evJ_conv = 6.24*pow(10, 18) # in eV/J


def noise_twist_power(moment, freq, fund_freq, pend_Q, Ptheta):
	# power = pow(2*math.pi, 4)*pow(moment*evJ_conv, 2)*(pow(pow(freq, 2) - pow(fund_freq, 2), 2) + pow(fund_freq, 4)*pow(pend_Q, -2))*Ptheta # J^2*(eV/J)^2 * rad^2/Hz
	power = pow(2*math.pi, 4)*pow(moment, 2)*(pow(pow(freq, 2) - pow(fund_freq, 2), 2) + pow(fund_freq, 4)*pow(pend_Q, -2))*Ptheta # J^2*(eV/J)^2 * rad^2/Hz
	return power

def noise_thermal_power(temp, kappa, freq, pend_Q):
	power = 48*boltzmann*temp*kappa*1.0/(2*math.pi*freq*pend_Q) #ev/K * K * eV/rad * s => eV^2/Hz
	return power

def background(sens, t, m, Q):
	coherence = Q/m
	background = 0
	if t <= coherence:
		background = sens*pow(t, -0.5)
		# print sens, sens*pow(t, -0.5)
	elif t > coherence:
		background = sens*pow(coherence, -0.5)*pow(coherence/t, 0.25)
		# print sens, sens*pow(coherence, -0.5), sens*pow(coherence, -0.5)*pow(coherence/t, 0.25)
	return background

# def coupling_form(n, p, mu, rho, v):
def axion_coupling_form(noise, Np, rho):
	# coupling = pow(10, -9)*2/(Np) * pow(0.3/(2*rho), 0.5)*noise*1/(pow(10*hbarc, 1.5)*5.4*pow(10, -13)) # in inverse GeV
	signal = 0.5*Np*pow(10, -3)*pow(2*rho*pow(10, 15)*pow(hbarc, 3.0), 0.5) # in eV^2
	coupling = pow(10, 9)*noise/signal # in 
	 # in inverse GeV
	return coupling

def axial_vector_coupling_form(mass, noise, Np, rho):
	# coupling = 1.0/(Np)*(mass/pow(10, 6))*pow(0.3/(2*rho), 0.5)*noise*1/(2.31*pow(10, 6)) #dimensionless
	# coupling = noise/(Np)*(mass*hbar)*pow(2*rho*pow(hbarc, 3)*pow(10, 15), -0.5) #dimensionless
	coupling = noise/(Np)*(mass*hbar)*pow(0.04, -2) #dimensionless
	return coupling


axion_masses = np.linspace(-22.0, -10.0, 100)
axion_result_1 = [0]*len(axion_masses)
axion_result_2 = [0]*len(axion_masses)
axion_result_3 = [0]*len(axion_masses)
axion_result_4 = [0]*len(axion_masses)
astrophysics_bounds = [0]*len(axion_masses)

vector_masses = np.linspace(-22.0, -10.0, 100)
vector_result_1 = [0]*len(vector_masses)
vector_result_2 = [0]*len(vector_masses)
vector_result_3 = [0]*len(vector_masses)
vector_astro = [0]*len(vector_masses)

axial_vector_masses = np.linspace(-22.0, -10.0, 100)
axial_vector_result_1 = [0]*len(axial_vector_masses)
axial_vector_result_2 = [0]*len(axial_vector_masses)
axial_vector_result_3 = [0]*len(axial_vector_masses)
axial_vector_astro = [0]*len(vector_masses)

# coupling_form(pow(10, 6), 30*24*60*60.0, 1.0, 0.1, 0.3, pow(10, -6))
integration_time = pow(10, 8) # 3 years :(
# integration_time = pow(10, 7) #4 months


#run 1 values
Ptheta_1 = pow(10, -18) # rad^2/Hz
# pend_mass_1 = 64*pow(10, -3) # kg
pend_mass_1 = 107*pow(10, -3) # kg
pend_radius_1 = .0545/2 # in m
# kappa_1 = 7.365*pow(10, 10) # eV/rad
kappa_1 = 1.155*pow(10, 11) # eV/rad
# freq_1 = 2.570*pow(10, -3) # in Hz
freq_1 = 5.379*pow(10, -3) # in Hz
# pend_moment_1 = pend_mass_1*pow(pend_radius_1, 2)*evJ_conv # in eV/s^2
pend_moment_1 = kappa_1/(pow(2*math.pi*freq_1, 2)) # in eV/s^2
print pend_mass_1*pow(pend_radius_1, 2)*evJ_conv, pend_moment_1


# pend_moment_2 = kappa/(4*pow(math.pi, 2)*pow(freq, 2)) # in eV*s^2
pend_Q_1 = 2000 # 
temp_1 = 300 # in Kelvin

#run 2 values
Ptheta_2 = pow(10, -24) # rad^2/Hz
# pend_mass_2 = 64*pow(10, -3) # kg
pend_mass_2 = 107*pow(10, -3) # kg
pend_radius_2 = .0545/2 # in m
# kappa_2 = 7.365*pow(10, 10) # eV/rad
kappa_2 = 1.155*pow(10, 11) # eV/rad
# freq_2 = 2.570*pow(10, -3) # in Hz
freq_2 = 5.379*pow(10, -3) # in Hz
pend_moment_2 = kappa_2/(pow(2*math.pi*freq_2, 2)) # in kg*m^2


# pend_moment_2 = kappa/(4*pow(math.pi, 2)*pow(freq, 2)) # in eV*s^2
pend_Q_2 = pow(10, 6) # don't know this for sure
temp_2 = 300 # in Kelvin

#run 3 values
Ptheta_3 = pow(10, -36) # rad^2/Hz
# Ptheta_3 = pow(10, -24) # rad^2/Hz
pend_mass_3 = 20 # kg
pend_radius_3 = .2 # in m
# pend_moment_3 = pend_mass_3*pow(pend_radius_3, 2) # in kg*m^2
kappa_3 = 5*pow(10, 15) # eV/rad
freq_3 = 100*pow(10, -3) # in Hz
# freq_3 = 5*pow(10, -3) # in Hz
# print freq_3, 1.0/(2*math.pi)*pow(kappa_3/(pend_moment_3*evJ_conv), 0.5)
pend_Q_3 = pow(10, 8) # don't know this for sure
temp_3 = 6 # in Kelvin
pend_moment_3 = kappa_3/(pow(2*math.pi*freq_3, 2)) # in kg*m^2
pend_moment_3_r = kappa_2/(pow(2*math.pi*freq_3, 2))


# total_systematic_uncertainty = 0.33 #order of magnitude estimate based on paper
total_systematic_uncertainty = 0.0 #order of magnitude estimate based on paper

Np = 1.0*pow(10, 23)
vel = pow(10, -3) # in c = 1
# rho = 3*pow(10, 14) # eV/m^3
rho = 0.3 # GeV/cm^3

W = 1.0
Q = pow(10, 6)

alpha = 1.0/137
mass_nucleon = 500*pow(10, 6) # in eV
temp = 30*pow(10, 6) # in eV
energy_loss_rate = 2.7*temp*1.8*pow(10, 38)*0.76*pow(temp, 3)/(3.0*pow(10, 14))
loss_rate_bound = pow(10, 19)*pow(10, -7)*evJ_conv
# astro_coupling = 6.8*pow(10, -11)*(mass_nucleon/temp)
vector_astro_coupling = pow(loss_rate_bound/energy_loss_rate, 0.5)*(mass_nucleon/temp)*pow(alpha, -0.5)
axion_astro_coupling = pow(loss_rate_bound/energy_loss_rate, 0.5)*(mass_nucleon/temp)*pow(alpha, -0.5)

astrophysics_bound = axion_astro_coupling


test_mass = 1.0*pow(10, 3) # in Hz
twist_noise_test = background(pow(noise_twist_power(pend_moment_1, test_mass, freq_1, pend_Q_1, Ptheta_1), 0.5), integration_time, test_mass, Q)
thermal_noise_test = background(pow(noise_thermal_power(temp_1, kappa_1, freq_1, pend_Q_1), 0.5), integration_time, test_mass, Q)
noise_test = max(twist_noise_test, thermal_noise_test)
axion_result_test = (axion_coupling_form(noise_test, Np, rho))
print math.log10(test_mass*hbar), axion_result_test



for i in xrange(0, len(axion_masses)):

	mass = pow(10, axion_masses[i])/hbar

	#run 1
	twist_noise_1 = background(pow(noise_twist_power(pend_moment_1, mass, freq_1, pend_Q_1, Ptheta_1), 0.5), integration_time, mass, Q)
	thermal_noise_1 = background(pow(noise_thermal_power(temp_1, kappa_1, freq_1, pend_Q_1), 0.5), integration_time, mass, Q)
	# thermal_noise_1 = background(pow(noise_thermal_power(temp_1, kappa_1, mass, pend_Q_1), 0.5), integration_time, mass, Q)
	noise_1 = max(twist_noise_1, thermal_noise_1, total_systematic_uncertainty)
	axion_result_1[i] = math.log10(axion_coupling_form(noise_1, Np, rho))

	#run 2
	twist_noise_2 = background(pow(noise_twist_power(pend_moment_2, mass, freq_2, pend_Q_2, Ptheta_2), 0.5), integration_time, mass, Q)
	thermal_noise_2 = background(pow(noise_thermal_power(temp_2, kappa_2, freq_2, pend_Q_2), 0.5), integration_time, mass, Q)
	# thermal_noise_2 = background(pow(noise_thermal_power(temp_2, kappa_2, mass, pend_Q_2), 0.5), integration_time, mass, Q)
	noise_2 = max(twist_noise_2, thermal_noise_2, total_systematic_uncertainty)
	axion_result_2[i] = math.log10(axion_coupling_form(noise_2, Np, rho))

	#run 3
	twist_noise_3 = background(pow(noise_twist_power(pend_moment_3, mass, freq_3, pend_Q_3, Ptheta_3), 0.5), integration_time, mass, Q)
	thermal_noise_3 = background(pow(noise_thermal_power(temp_3, kappa_3, freq_3, pend_Q_3), 0.5), integration_time, mass, Q)
	# thermal_noise_3 = background(pow(noise_thermal_power(temp_3, kappa_3, mass, pend_Q_3), 0.5), integration_time, mass, Q)
	noise_3 = max(twist_noise_3, thermal_noise_3, total_systematic_uncertainty)
	axion_result_3[i] = math.log10(axion_coupling_form(noise_3, Np*pend_mass_3/pend_mass_1, rho))

	twist_noise_4 = background(pow(noise_twist_power(pend_moment_3_r, mass, freq_3, pend_Q_3, Ptheta_3), 0.5), integration_time, mass, Q)
	thermal_noise_4 = background(pow(noise_thermal_power(temp_3, kappa_2, freq_3, pend_Q_3), 0.5), integration_time, mass, Q)
	# thermal_noise_3 = background(pow(noise_thermal_power(temp_3, kappa_3, mass, pend_Q_3), 0.5), integration_time, mass, Q)
	noise_4 = max(twist_noise_4, thermal_noise_4, total_systematic_uncertainty)
	axion_result_4[i] = math.log10(axion_coupling_form(noise_4, Np, rho))

	# print noise_1, noise_2, noise_3

	astrophysics_bounds[i] = math.log10(astrophysics_bound)


fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

ax.plot(axion_masses, axion_result_1)
ax.plot(axion_masses, axion_result_2)
ax.plot(axion_masses, axion_result_3)
ax.plot(axion_masses, axion_result_4)
ax.plot(axion_masses, astrophysics_bounds)
ax.legend(("Run 1", "Run 2", "Run 3", "Run 3 (constant R)", "Astrophysics"),  loc=2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('$g_{aee}$ (GeV)$^{-1}$ (log scale)');
ax.set_title('axion sensitivity')
# ax.legend(("G", "E", "CalcG", "CalcE"));
plt.savefig("torsion_pend_sens.pdf")
plt.show()

# now for the axial vector sensitivity
mass_nucleon = 500*pow(10, 6) # in eV
temp = 30*pow(10, 6) # in eV
energy_loss_rate = 2.7*temp*1.8*pow(10, 38)*0.76*pow(temp, 3)/(3*pow(10, 14))
loss_rate_bound = pow(10, 19)*pow(10, -7)*evJ_conv
# astro_coupling = 6.8*pow(10, -11)*(mass_nucleon/temp)
astro_coupling = pow(loss_rate_bound/energy_loss_rate, 0.5)*(mass_nucleon/temp)

vector_astro_bound = astro_coupling

for i in xrange(0, len(vector_masses)):

	mass = pow(10, vector_masses[i])/hbar

	vector_astro[i] = math.log10(vector_astro_bound)

	#run 1
	twist_noise_1 = background(pow(noise_twist_power(pend_moment_1, mass, freq_1, pend_Q_1, Ptheta_1), 0.5), integration_time, mass, Q)
	thermal_noise_1 = background(pow(noise_thermal_power(temp_1, kappa_1, freq_1, pend_Q_1), 0.5), integration_time, mass, Q)
	noise_1 = max(twist_noise_1, thermal_noise_1)
	vector_result_1[i] = math.log10(axion_coupling_form(noise_1, Np, rho))

	#run 2
	twist_noise_2 = background(pow(noise_twist_power(pend_moment_2, mass, freq_2, pend_Q_2, Ptheta_2), 0.5), integration_time, mass, Q)
	thermal_noise_2 = background(pow(noise_thermal_power(temp_2, kappa_2, freq_2, pend_Q_2), 0.5), integration_time, mass, Q)
	noise_2 = max(twist_noise_2, thermal_noise_2)
	vector_result_2[i] = math.log10(axion_coupling_form(noise_2, Np, rho))

	#run 3
	twist_noise_3 = background(pow(noise_twist_power(pend_moment_3, mass, freq_3, pend_Q_3, Ptheta_3), 0.5), integration_time, mass, Q)
	thermal_noise_3 = background(pow(noise_thermal_power(temp_3, kappa_3, freq_3, pend_Q_3), 0.5), integration_time, mass, Q)
	noise_3 = max(twist_noise_3, thermal_noise_3)
	vector_result_3[i] = math.log10(axion_coupling_form(noise_3, Np*pend_mass_3/pend_mass_1, rho))

fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

# ax.plot(vector_masses, vector_result_1)
# ax.plot(vector_masses, vector_result_2)
# ax.plot(vector_masses, vector_result_3)
ax.plot(axion_masses, axion_result_1)
ax.plot(axion_masses, axion_result_2)
ax.plot(axion_masses, axion_result_3)
ax.plot(axion_masses, vector_astro)
ax.legend(("Run 1", "Run 2", "Run 3", "Astro"), loc=2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('$g_{ANN}$ (log scale)');
ax.set_title('hidden photon sensitivity')
# ax.legend(("G", "E", "CalcG", "CalcE"));
plt.show()


for i in xrange(0, len(vector_masses)):

	mass = pow(10, axial_vector_masses[i])/hbar

	vector_astro[i] = math.log10(vector_astro_bound)

	#run 1
	twist_noise_1 = background(pow(noise_twist_power(pend_moment_1, mass, freq_1, pend_Q_1, Ptheta_1), 0.5), integration_time, mass, Q)
	thermal_noise_1 = background(pow(noise_thermal_power(temp_1, kappa_1, freq_1, pend_Q_1), 0.5), integration_time, mass, Q)
	noise_1 = max(twist_noise_1, thermal_noise_1)
	axial_vector_result_1[i] = math.log10(axial_vector_coupling_form(mass, noise_1, Np, rho))

	#run 2
	twist_noise_2 = background(pow(noise_twist_power(pend_moment_2, mass, freq_2, pend_Q_2, Ptheta_2), 0.5), integration_time, mass, Q)
	thermal_noise_2 = background(pow(noise_thermal_power(temp_2, kappa_2, freq_2, pend_Q_2), 0.5), integration_time, mass, Q)
	noise_2 = max(twist_noise_2, thermal_noise_2)
	axial_vector_result_2[i] = math.log10(axial_vector_coupling_form(mass, noise_2, Np, rho))

	#run 3
	twist_noise_3 = background(pow(noise_twist_power(pend_moment_3, mass, freq_3, pend_Q_3, Ptheta_3), 0.5), integration_time, mass, Q)
	thermal_noise_3 = background(pow(noise_thermal_power(temp_3, kappa_3, freq_3, pend_Q_3), 0.5), integration_time, mass, Q)
	noise_3 = max(twist_noise_3, thermal_noise_3)
	axial_vector_result_3[i] = math.log10(axial_vector_coupling_form(mass, noise_3, Np*pend_mass_3/pend_mass_1, rho))

fig, ax = plt.subplots()
# handles, labels = ax.get_legend_handles_labels()
# ax.legend(handles, labels)

# ax.plot(vector_masses, vector_result_1)
# ax.plot(vector_masses, vector_result_2)
# ax.plot(vector_masses, vector_result_3)
ax.plot(axion_masses, axial_vector_result_1)
ax.plot(axion_masses, axial_vector_result_2)
ax.plot(axion_masses, axial_vector_result_3)
# ax.plot(axion_masses, vector_astro)
ax.legend(("Run 1", "Run 2", "Run 3"), loc=2);
ax.set_xlabel('mass (eV) (log scale)');
ax.set_ylabel('$g_{ax}$ (log scale)');
ax.set_title('axial vector sensitivity')
# ax.legend(("G", "E", "CalcG", "CalcE"));
# plt.show()