#!/usr/bin/env python3
"""
RAM EFR Analysis Script - Parallel Processing Pipeline

This script provides a parallel processing pipeline for RAM EFR analysis:
1. Processes multiple subjects with custom audiogram data
2. Converts audiogram data to Poles using OHC_ind function for each subject
3. Runs simulations with user-specified ANF distributions (HSR, MSR, LSR)
4. Calculates EFR values for each simulation
5. Saves results to CSV file with subject names and EFR values

Usage:
    Configure your settings in the main block and run the script.
    The script will process all subjects in parallel using all available CPU cores.

Created by: Brent Nissens
Date: October 22, 2025
"""

import numpy as np
import os
import sys
import pandas as pd
import multiprocessing as mp
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from get_RAM_stims import get_RAM_stims
from model2018 import model2018
import OHC_ind

def load_shera_poles_profile(sheraP_folder):
    """
    Load Shera poles profile from a specific path.
    """
    try:
        sheraP_data = np.loadtxt(sheraP_folder + '/StartingPoles.dat')
        # Take only the first row if there are multiple rows
        if sheraP_data.ndim > 1:
            sheraP = sheraP_data[0, :] 
        else:
            sheraP = sheraP_data
        return sheraP
        
    except Exception as e:
        print(f"Error loading Shera poles profile from {sheraP_folder}: {e}")
        return None

def run_simulation(stim, fs, sheraP, HSR, MSR, LSR):
    """
    Run a simulation with the given parameters.
    """
    try:
        output = model2018(stim, fs, 'abr', 1, 'evihmlbw', 1, sheraP, 0.05, 'vel', HSR, MSR, LSR, 1, os.getcwd())
        return output
    except Exception as e:
        print(f"Error running simulation: {e}")
        return None

def calculate_EFR(output):
    """
    Calculate the EFR from the output.
    """
    # This function computes the EFR (Envelope Following Response) harmonics sum from simulation output.
    # It mirrors get_RAM_EFRS1, but omits plotting and converts to microV.

    try:
        # Handle the case where output is a list containing a dictionary
        if isinstance(output, list) and len(output) > 0:
            output = output[0]  # Get the first (and likely only) element
        
        # Extract relevant waveforms and sampling frequency from the output object/struct.
        fs = float(output.fs_an)
        w1 = output.w1.flatten()
        w3 = output.w3.flatten()
        w5 = output.w5.flatten()

        EFR = w1 + w3 + w5

        # Fourier Transform
        L = len(EFR)
        Y = np.fft.fft(EFR)
        P2 = np.abs(Y / L)
        P1 = P2[:L//2 + 1]
        P1[1:-1] = 2 * P1[1:-1]
        f = fs * np.arange(L//2 + 1) / L

        # Find indices for 4 harmonics of 110 Hz
        fundamental = 110  # Hz
        num_harmonics = 4
        harmonics = np.arange(1, num_harmonics + 1) * fundamental
        idx = []
        for harmonic in harmonics:
            idx.append(np.argmin(np.abs(f - harmonic)))
        # Calculate sum and convert to microV
        harmonic_sum = np.sum(P1[idx]) * 1e6

        return harmonic_sum

    except Exception as e:
        print(f"Error calculating EFR: {e}")
        return np.nan

def process_single_subject(args):
    """
    Process a single subject's simulation.
    
    Parameters:
    -----------
    args : tuple
        (subject_data, HSR, MSR, LSR, stim, fs, poles_output_dir)
        where subject_data is (idx, subject_name, hl_freqs_hz, hl_db)
        HSR, MSR, LSR can be scalars or arrays (frequency-dependent ANF distributions)
    
    Returns:
    --------
    dict
        Result row with subject name and EFR value
    """
    subject_data, HSR, MSR, LSR, stim, fs, poles_output_dir = args
    idx, subject_name, hl_freqs_hz, hl_db = subject_data
    
    try:
        # Create poles using OHC_ind (without showing figures)
        OHC_ind.ohc_ind(
            name=subject_name,
            hl_freqs_hz=hl_freqs_hz,
            hl_db=hl_db,
            base_dir=poles_output_dir,
            show_figs=False
        )
        
        # Load the poles profile
        subject_poles_path = os.path.join(poles_output_dir, 'Poles', subject_name)
        sheraP = load_shera_poles_profile(subject_poles_path)
        if sheraP is None:
            print(f"  Could not load poles for {subject_name}, skipping...")
            return None
        
        # Run simulation with ANF distributions
        # Note: HSR, MSR, LSR can be scalars (constant across frequency) or 
        # arrays (frequency-dependent ANF distributions)
        output = run_simulation(stim, fs, sheraP, HSR, MSR, LSR)
        
        if output is None:
            return None
        
        # Calculate EFR
        efr = calculate_EFR(output)
        
        # Store results
        result_row = {
            'Name': subject_name,
            'EFR': efr
        }
        
        return result_row
        
    except Exception as e:
        print(f"Error processing subject {subject_name}: {e}")
        return None


if __name__ == "__main__":
    # ============================================================================
    # CONFIGURATION - Modify these settings for your use case
    # ============================================================================
    
    # Path to Excel file containing audiogram data
    # Expected format: Excel file with columns 'ID' and 'Audio_XXXHz' (e.g., 'Audio_125Hz', 'Audio_250Hz')
    excel_path = './data/audiograms.xlsx'
    
    # ANF distribution settings
    # HSR, MSR, LSR can be scalars (constant across frequency) or arrays (frequency-dependent)
    # For frequency-dependent distributions, provide arrays with values for each frequency channel
    # Default values: HSR=13, MSR=3, LSR=3 (constant across all frequencies)
    HSR = 13  # High Spontaneous Rate fibers
    MSR = 3   # Medium Spontaneous Rate fibers
    LSR = 3   # Low Spontaneous Rate fibers
    
    # Poles output directory (where OHC_ind will save generated poles)
    poles_output_dir = '.'
    
    # Stimulus parameters
    fs = 1e5  # Sampling frequency in Hz
    fRAM = np.array([4000])  # RAM frequency in Hz
    
    # Output settings
    output_csv = 'EFR_results.csv'  # Output CSV filename
    
    # Parallel processing settings
    num_workers = mp.cpu_count()  # Number of parallel workers (default: all CPU cores)
    
    # ============================================================================
    # MAIN PROCESSING PIPELINE
    # ============================================================================
    
    print("="*60)
    print("Starting RAM EFR analysis pipeline...")
    print("="*60 + "\n")
    
    # Load Excel data
    print("="*60)
    print(f"Loading data from {excel_path}")
    print("="*60 + "\n")

    try:
        df = pd.read_excel(excel_path)
    except Exception as e:
        print(f"Error loading Excel file: {e}")
        sys.exit(1)
    
    # Get audiogram column names (assuming they are Audio_XXXHz)
    audio_columns = [col for col in df.columns if 'Audio_' in col]
    
    if not audio_columns:
        print("Error: No audiogram columns found (expected format: 'Audio_XXXHz')")
        sys.exit(1)
    
    # Initialize stimulus
    stim = get_RAM_stims(fs, fRAM)
    print("="*60)
    print(f"Generated RAM stimulus with fs={fs} Hz, fRAM={fRAM} Hz")
    print("="*60 + "\n")

    # Prepare all subject data for parallel processing
    subject_data_list = []
    for idx, row in df.iterrows():
        subject_name = str(row.get('ID', f'Subject_{idx}'))
        
        # Extract audiogram frequencies and values
        audiogram_data = []
        for col in audio_columns:
            freq_hz = int(col.replace('Audio_', '').replace('Hz', ''))
            hl_db = row[col]
            audiogram_data.append((freq_hz, hl_db))
        
        # Sort by frequency
        audiogram_data.sort(key=lambda x: x[0])
        hl_freqs_hz = [x[0] for x in audiogram_data]
        hl_db = [x[1] for x in audiogram_data]
        
        subject_data_list.append((idx, subject_name, hl_freqs_hz, hl_db))
    
    print("="*60)
    print(f"Processing {len(subject_data_list)} subjects using {num_workers} workers")
    print(f"ANF distributions: HSR={HSR}, MSR={MSR}, LSR={LSR}")
    print("="*60 + "\n")
    
    # Prepare arguments for each worker
    worker_args = [(sd, HSR, MSR, LSR, stim, fs, poles_output_dir) for sd in subject_data_list]
    
    # Process subjects in parallel
    start_time = time.time()
    completed_subjects = 0
    results = []
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {executor.submit(process_single_subject, args): idx for idx, args in enumerate(worker_args)}
        
        # Collect results with progress bar
        pbar = tqdm(total=len(futures), desc="Processing subjects")
        
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                results.append(result)
            
            # Update progress and time estimates
            completed_subjects += 1
            elapsed = time.time() - start_time
            avg_time_per_subject = elapsed / completed_subjects if completed_subjects > 0 else 0
            remaining_subjects = len(futures) - completed_subjects
            estimated_time_left = avg_time_per_subject * remaining_subjects / num_workers if num_workers > 0 else 0
            
            # Format time estimates
            elapsed_str = f"{int(elapsed // 3600)}h {int((elapsed % 3600) // 60)}m {int(elapsed % 60)}s"
            eta_str = f"{int(estimated_time_left // 3600)}h {int((estimated_time_left % 3600) // 60)}m {int(estimated_time_left % 60)}s"
            
            pbar.set_postfix({
                'Elapsed': elapsed_str,
                'ETA': eta_str,
                'Avg/Subj': f"{avg_time_per_subject:.1f}s"
            })
            pbar.update(1)
        
        pbar.close()
    
    # Calculate total processing time
    total_time = time.time() - start_time
    hours = int(total_time // 3600)
    minutes = int((total_time % 3600) // 60)
    seconds = int(total_time % 60)
    
    # Convert results to DataFrame and save to CSV
    if results:
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_csv, index=False)
        print(f"\nResults saved to {output_csv}")
        print(f"\nTotal subjects processed: {len(results)}")
        print(f"Total processing time: {hours}h {minutes}m {seconds}s")
        print("\nResults summary:")
        print(results_df)
    else:
        print("\nWarning: No results were generated. Check your configuration and input data.")