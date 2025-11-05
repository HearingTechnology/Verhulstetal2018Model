import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy import signal
import os

# This file runs and plots the results of a click simulation performed with
# ExampleSimulation.py
# Load the simulation results
data = sio.loadmat('Simulations.mat')
output = data['output'][0]  # MATLAB struct becomes a structured array

# Define the different conditions
L = [0]
p0 = 2e-5
CF = output['cf'].item().flatten()  # Handle nested array structure

# Create figure for CF plot
plt.figure()
plt.plot(CF)
plt.xlabel('Cochlear Channel Number [-]')
plt.ylabel('Characteristic Frequency [Hz]')
plt.title('Characteristic Frequencies')
plt.grid(True)
plt.show()

# Sampling frequencies
fs_c = output['fs_bm'].item().item()  # BM, OAE and IHC sampling frequency (higher to avoid numerical errors)
fs = output['fs_an'].item().item()    # AN, CN, IC and waves sampling frequency (5 times lower)

# Time vectors
t_c = np.arange(output['v'].item().shape[0]) / fs_c
t = np.arange(output['anfH'].item().shape[0]) / fs

# Frequency vectors
f_c = np.arange(output['v'].item().shape[0]) * fs_c / output['v'].item().shape[0]
f = np.arange(output['anfH'].item().shape[0]) * fs / output['anfH'].item().shape[0]

# Pick a channel number to plot results from
No = 245

# Reorganization of the data for easier processing
num_conditions = len(L)
v_data = output['v'].item()
ihc_data = output['ihc'].item()
e_data = output['e'].item()

oae = np.zeros((e_data.shape[1], num_conditions))
vrms = np.zeros((v_data.shape[1], num_conditions))
ihcrms = np.zeros((ihc_data.shape[1], num_conditions))
v = np.zeros((v_data.shape[0], num_conditions))
ihc = np.zeros((ihc_data.shape[0], num_conditions))

for n in range(num_conditions):
    oae[:, n] = e_data.flatten()
    vrms[:, n] = np.sqrt(np.mean(v_data**2, axis=0))
    ihcrms[:, n] = np.sqrt(np.mean(ihc_data**2, axis=0))
    v[:, n] = v_data[:, No]
    ihc[:, n] = ihc_data[:, No]

# The OAE (Otoacoustic Emission)
# This figure shows the ear-canal pressure:
# To get the reflection-source OAE, do the following simulation:
# OAE_{reflections on} - OAE_{reflections off}
# To get the distortion-source OAE, do a normalization using a low level
# (linear) simulation with reflections off
plt.figure(figsize=(10, 8))

plt.subplot(2, 1, 1)
plt.plot(1000 * t_c, oae)
plt.xlabel('Time [ms]')
plt.ylabel('Ear Canal Pressure [Pa]')
plt.xlim([0, 20])
plt.ylim([-0.02, 0.02])
plt.legend([str(l) for l in L], frameon=False)
plt.title('OAE Time Domain')

# Zero padding to plot OAE spectrum
oae_spectrum = oae.copy()
oae_spectrum[:200, :] = 0

plt.subplot(2, 1, 2)
oae_fft = np.fft.fft(oae_spectrum / p0, axis=0)
plt.plot(f_c / 1000, 20 * np.log10(np.abs(oae_fft)))
plt.xlabel('Frequency [kHz]')
plt.ylabel('EC Magnitude [dB re p0]')
plt.xlim([0, 12])
plt.legend([str(l) for l in L], frameon=False)
plt.title('OAE Frequency Domain')
plt.grid(True)
plt.tight_layout()
plt.show()

# v_bm and V_IHC
plt.figure(figsize=(12, 10))

plt.subplot(2, 2, 1)
plt.plot(1000 * t_c, v)
plt.xlabel('Time [ms]')
plt.ylabel('v_{bm} [m/s]')
plt.xlim([0, 30])
plt.legend([str(l) for l in L], frameon=False)
plt.title(f'CF of {round(CF[No])} Hz')
plt.grid(True)

plt.subplot(2, 2, 2)
plt.plot(CF / 1000, 20 * np.log10(vrms))
plt.xlabel('CF [kHz]')
plt.ylabel('rms of v_{bm} [dB re 1 m/s]')
plt.xlim([0, 14])
plt.ylim([np.max(20 * np.log10(vrms)) - 100, np.max(20 * np.log10(vrms)) + 10])
plt.title('Excitation Pattern')
plt.grid(True)

plt.subplot(2, 2, 3)
plt.plot(1000 * t_c, ihc)
plt.xlabel('Time [ms]')
plt.ylabel('V_{ihc} [V]')
plt.xlim([0, 30])
plt.legend([str(l) for l in L], frameon=False)
plt.title(f'CF of {round(CF[No])} Hz')
plt.grid(True)

plt.subplot(2, 2, 4)
plt.plot(CF / 1000, ihcrms)
plt.xlabel('CF [kHz]')
plt.ylabel('rms of V_{ihc} [V]')
plt.xlim([0, 14])
plt.title('Excitation Pattern')
plt.grid(True)

plt.tight_layout()
plt.show()

# Reorganization of the data for easier processing
anfH_data = output['anfH'].item()
anfM_data = output['anfM'].item()
anfL_data = output['anfL'].item()
an_summed_data = output['an_summed'].item()
cn_data = output['cn'].item()
ic_data = output['ic'].item()
w1_data = output['w1'].item()
w3_data = output['w3'].item()
w5_data = output['w5'].item()

HSR = np.zeros((anfH_data.shape[0], num_conditions))
MSR = np.zeros((anfM_data.shape[0], num_conditions))
LSR = np.zeros((anfL_data.shape[0], num_conditions))
AN = np.zeros((an_summed_data.shape[0], num_conditions))
CN = np.zeros((cn_data.shape[0], num_conditions))
IC = np.zeros((ic_data.shape[0], num_conditions))
W1 = np.zeros((w1_data.shape[1], num_conditions))  # Use shape[1] for wave data
W3 = np.zeros((w3_data.shape[1], num_conditions))
W5 = np.zeros((w5_data.shape[1], num_conditions))
EFR = np.zeros((w1_data.shape[1], num_conditions))

for n in range(num_conditions):
    HSR[:, n] = anfH_data[:, No]
    MSR[:, n] = anfM_data[:, No]
    LSR[:, n] = anfL_data[:, No]
    AN[:, n] = an_summed_data[:, No]
    CN[:, n] = cn_data[:, No]
    IC[:, n] = ic_data[:, No]
    W1[:, n] = w1_data.flatten()
    W3[:, n] = w3_data.flatten()
    W5[:, n] = w5_data.flatten()
    EFR[:, n] = w1_data.flatten() + w3_data.flatten() + w5_data.flatten()

# Single unit responses
plt.figure(figsize=(12, 10))

plt.subplot(3, 2, 1)
plt.plot(1000 * t, HSR)
plt.title(f'CF of {round(CF[No])} Hz')
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('HSR fiber [spikes/s]')
plt.legend([str(l) for l in L], frameon=False)
plt.grid(True)

plt.subplot(3, 2, 3)
plt.plot(1000 * t, MSR)
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('MSR fiber [spikes/s]')
plt.grid(True)

plt.subplot(3, 2, 5)
plt.plot(1000 * t, LSR)
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('LSR fiber [spikes/s]')
plt.grid(True)

plt.subplot(3, 2, 2)
plt.plot(1000 * t, AN)
plt.title(f'CF of {round(CF[No])} Hz')
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('sum AN [spikes/s]')
plt.grid(True)

plt.subplot(3, 2, 4)
plt.plot(1000 * t, CN)
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('CN [spikes/s]')
plt.grid(True)

plt.subplot(3, 2, 6)
plt.plot(1000 * t, IC)
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('IC [spikes/s]')
plt.grid(True)

plt.tight_layout()
plt.show()

# Population responses
plt.figure(figsize=(10, 12))

plt.subplot(4, 1, 1)
plt.plot(1000 * t, 1e6 * W1)
plt.title('Population Responses summed across simulated CFs')
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('W-1 [μV]')
plt.legend([str(l) for l in L], frameon=False)
plt.grid(True)

plt.subplot(4, 1, 2)
plt.plot(1000 * t, 1e6 * W3)
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('W-3 [μV]')
plt.grid(True)

plt.subplot(4, 1, 3)
plt.plot(1000 * t, 1e6 * W5)
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('W-5 [μV]')
plt.grid(True)

plt.subplot(4, 1, 4)
plt.plot(1000 * t, 1e6 * EFR)
plt.xlim([0, 20])
plt.xlabel('Time [ms]')
plt.ylabel('EFR [μV]')
plt.grid(True)

plt.tight_layout()
plt.show()

print("Analysis complete!")