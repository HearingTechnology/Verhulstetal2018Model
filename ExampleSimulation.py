#!/usr/bin/env python3
"""
Easy Model Runner with EFR Calculation

This script provides a simple interface to run the model2018 and calculate EFR (Envelope Following Response).

Usage:
    python easy_model_run.py

Features:
- Generates RAM stimulus automatically
- Loads poles profile (configurable)
- Runs model2018 simulation
- Calculates EFR using FFT analysis
- Displays results clearly
- Saves output to files

Created by: Brent Nissens
Based on: Verhulst et al. 2018 model and FullSimulationRAM.py
"""

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import os
from get_RAM_stims import get_RAM_stims
from model2018 import model2018

def calculate_EFR(output):
    """
    Calculate the EFR from model output using FFT analysis.
    
    This function computes the EFR (Envelope Following Response) harmonics sum 
    from simulation output, similar to the calculation in FullSimulationRAM.py.
    
    Parameters:
    -----------
    output : ModelOutput
        Output from model2018 function
        
    Returns:
    --------
    float
        EFR value in microvolts (μV)
    """
    try:
        # Handle the case where output is a list containing a dictionary
        if isinstance(output, list) and len(output) > 0:
            output = output[0]  # Get the first (and likely only) element
        
        # Extract relevant waveforms and sampling frequency from the output object/struct.
        fs = float(output.fs_an)
        w1 = output.w1.flatten()
        w3 = output.w3.flatten()
        w5 = output.w5.flatten()

        # Combine waves to get EFR
        EFR = w1 + w3 + w5

        # Fourier Transform
        L = len(EFR)
        Y = np.fft.fft(EFR)
        P2 = np.abs(Y / L)
        P1 = P2[:L//2 + 1]
        P1[1:-1] = 2 * P1[1:-1]
        f = fs * np.arange(L//2 + 1) / L

        # Find indices for 4 harmonics of 110 Hz (fundamental modulation frequency)
        fundamental = 110  # Hz
        num_harmonics = 4
        harmonics = np.arange(1, num_harmonics + 1) * fundamental
        idx = []
        for harmonic in harmonics:
            idx.append(np.argmin(np.abs(f - harmonic)))
        
        # Calculate sum of harmonics and convert to microV
        harmonic_sum = np.sum(P1[idx]) * 1e6

        return harmonic_sum, f, P1, harmonics, idx

    except Exception as e:
        print(f"Error calculating EFR: {e}")
        return np.nan, None, None, None, None

def run_easy_model(carrier_freq=4000, poles_profile='Flat00', show_plots=True, save_results=True):
    """
    Easy interface to run the model2018 and calculate EFR.
    
    Parameters:
    -----------
    carrier_freq : float, optional
        Carrier frequency for RAM stimulus in Hz (default: 4000)
    poles_profile : str, optional
        Name of poles profile folder in ./Poles/ directory (default: 'Flat00')
    show_plots : bool, optional
        Whether to show plots of results (default: True)
    save_results : bool, optional
        Whether to save results to MAT files (default: True)
        
    Returns:
    --------
    dict
        Dictionary containing:
        - 'efr_value': EFR value in μV
        - 'output': Full model output
        - 'stimulus': Generated stimulus
        - 'frequency_spectrum': Frequency array for FFT
        - 'power_spectrum': Power spectrum
    """
    
    print("="*60)
    print("Easy Model2018 Runner with EFR Calculation")
    print("="*60)
    
    # Configuration
    fs = 1e5  # Sampling frequency
    fRAM = np.array([carrier_freq])
    
    print(f"Carrier frequency: {carrier_freq} Hz")
    print(f"Sampling frequency: {fs} Hz")
    print(f"Poles profile: {poles_profile}")
    print()
    
    # Generate stimulus
    print("Generating RAM stimulus...")
    try:
        stimulus = get_RAM_stims(fs, fRAM)
        print(f"✓ Generated stimulus: {stimulus.shape[1]} samples, duration: {stimulus.shape[1]/fs:.3f} s")
    except Exception as e:
        print(f"✗ Error generating stimulus: {e}")
        return None
    
    # Load poles profile
    poles_path = f'./Poles/{poles_profile}/StartingPoles.dat'
    print(f"Loading poles from: {poles_path}")
    try:
        sheraP = np.loadtxt(poles_path)
        # Take only the first row if there are multiple rows
        if sheraP.ndim > 1:
            sheraP = sheraP[0, :]
        print(f"✓ Loaded poles profile: {len(sheraP)} values")
    except Exception as e:
        print(f"✗ Error loading poles: {e}")
        print(f"Available poles profiles in ./Poles/:")
        try:
            poles_dirs = [d for d in os.listdir('./Poles/') if os.path.isdir(f'./Poles/{d}')]
            for d in sorted(poles_dirs):
                print(f"  - {d}")
        except:
            print("  Could not list poles directories")
        return None
    
    # Run model
    print("\nRunning model2018...")
    try:
        results = model2018(
            stimulus, 
            fs, 
            fc='abr',                    # Use ABR frequency points
            irregularities=1,            # Enable irregularities
            storeflag='evihmlbw',        # Store all relevant outputs
            subject=1,                   # Subject seed
            sheraPo=sheraP,              # Poles profile
            IrrPct=0.05,                 # 5% irregularities
            non_linear_type='vel',       # Velocity-based nonlinearity
            nH=13,                       # High spontaneous rate fibers
            nM=3,                        # Medium spontaneous rate fibers
            nL=3,                        # Low spontaneous rate fibers
            clean=1,                     # Clean parameter
            data_folder='./'             # Data folder
        )
        print("✓ Model simulation completed successfully!")
    except Exception as e:
        print(f"✗ Error running model: {e}")
        return None
    
    # Calculate EFR
    print("\nCalculating EFR...")
    try:
        efr_value, f, P1, harmonics, harmonic_idx = calculate_EFR(results)
        if not np.isnan(efr_value):
            print(f"✓ EFR calculated: {efr_value:.4f} μV")
        else:
            print("✗ EFR calculation failed")
            return None
    except Exception as e:
        print(f"✗ Error calculating EFR: {e}")
        return None
    
    # Display harmonic details
    print("\nHarmonic Analysis:")
    for i, (harm_freq, idx) in enumerate(zip(harmonics, harmonic_idx)):
        power_uv = P1[idx] * 1e6
        print(f"  {harm_freq} Hz (harmonic {i+1}): {power_uv:.4f} μV")
    
    # Save results if requested
    if save_results:
        print("\nSaving results...")
        try:
            output = results[0]
            
            # Save EFR and waves
            EFR_combined = output.w1 + output.w3 + output.w5
            sio.savemat('easy_model_EFR.mat', {
                'EFR': EFR_combined,
                'w1': output.w1,
                'w3': output.w3,
                'w5': output.w5,
                'efr_value_uV': efr_value,
                'fs': output.fs_an,
                'carrier_freq': carrier_freq,
                'poles_profile': poles_profile
            })
            
            # Save velocity and other outputs
            sio.savemat('easy_model_output.mat', {
                'v': output.v,
                'cf': output.cf,
                'stimulus': stimulus,
                'fs_bm': output.fs_bm
            })
            
            print("✓ Results saved to:")
            print("  - easy_model_EFR.mat (EFR and waves)")
            print("  - easy_model_output.mat (full model output)")
            
        except Exception as e:
            print(f"✗ Error saving results: {e}")
    
    # Create plots if requested
    if show_plots and f is not None:
        print("\nGenerating plots...")
        try:
            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            fig.suptitle(f'Model2018 Results - Carrier: {carrier_freq} Hz, Profile: {poles_profile}', fontsize=14)
            
            # Plot 1: Stimulus
            t_stim = np.arange(len(stimulus[0])) / fs
            axes[0,0].plot(t_stim[:1000], stimulus[0][:1000])  # Show first 1000 samples
            axes[0,0].set_title('RAM Stimulus (first 10 ms)')
            axes[0,0].set_xlabel('Time (s)')
            axes[0,0].set_ylabel('Amplitude')
            axes[0,0].grid(True)
            
            # Plot 2: EFR waveform
            output = results[0]
            EFR_combined = output.w1 + output.w3 + output.w5
            t_efr = np.arange(len(EFR_combined)) / output.fs_an
            axes[0,1].plot(t_efr, EFR_combined)
            axes[0,1].set_title('EFR Waveform (w1+w3+w5)')
            axes[0,1].set_xlabel('Time (s)')
            axes[0,1].set_ylabel('Amplitude')
            axes[0,1].grid(True)
            
            # Plot 3: Frequency spectrum
            axes[1,0].semilogy(f[:2000], P1[:2000])  # Show up to 400 Hz
            for i, idx in enumerate(harmonic_idx):
                axes[1,0].semilogy(f[idx], P1[idx], 'ro', markersize=8, 
                                 label=f'{harmonics[i]} Hz')
            axes[1,0].set_title('FFT Spectrum with Harmonics')
            axes[1,0].set_xlabel('Frequency (Hz)')
            axes[1,0].set_ylabel('Power')
            axes[1,0].grid(True)
            axes[1,0].legend()
            axes[1,0].set_xlim(0, 500)
            
            # Plot 4: Individual waves
            axes[1,1].plot(t_efr, output.w1, label='Wave 1', alpha=0.7)
            axes[1,1].plot(t_efr, output.w3, label='Wave 3', alpha=0.7)
            axes[1,1].plot(t_efr, output.w5, label='Wave 5', alpha=0.7)
            axes[1,1].set_title('Individual ABR Waves')
            axes[1,1].set_xlabel('Time (s)')
            axes[1,1].set_ylabel('Amplitude')
            axes[1,1].legend()
            axes[1,1].grid(True)
            
            plt.tight_layout()
            plt.savefig('easy_model_results.png', dpi=150, bbox_inches='tight')
            plt.show()
            
            print("✓ Plots saved as 'easy_model_results.png'")
            
        except Exception as e:
            print(f"✗ Error creating plots: {e}")
    
    # Prepare return dictionary
    result_dict = {
        'efr_value': efr_value,
        'output': results[0],
        'stimulus': stimulus,
        'frequency_spectrum': f,
        'power_spectrum': P1,
        'harmonics': harmonics,
        'carrier_freq': carrier_freq,
        'poles_profile': poles_profile
    }
    
    print(f"\n" + "="*60)
    print(f"SUMMARY: EFR = {efr_value:.4f} μV")
    print("="*60)
    
    return result_dict

if __name__ == "__main__":
    # Example usage with different configurations
    
    print("Running easy model with default parameters...")
    results = run_easy_model()
    
    if results is not None:
        print(f"\nSuccess! EFR value: {results['efr_value']:.4f} μV")
        
        # You can also run with different parameters:
        # results = run_easy_model(carrier_freq=2000, poles_profile='Normal', show_plots=False)
        # results = run_easy_model(carrier_freq=8000, poles_profile='Flat10', save_results=False)
    else:
        print("\nModel run failed. Please check the error messages above.")
