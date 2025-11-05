import numpy as np
from scipy import signal

def get_RAM_stims(fs, fRAM):
    """
    Generating stimulation waveform for RAM (Rectangular Amplitude Modulation)
    Based on Sarineh Keshishzadeh 10.11.2021.
    
    Parameters:
    -----------
    fs : float
        Sampling frequency [Hz]
    fRAM : array_like
        [Mx1] vector of RAM carrier frequencies [Hz]
        M - number of stimuli to generate
    
    Returns:
    --------
    stim : numpy.ndarray
        [MxN] matrix of the generated stimuli
        N - number of samples in the stimuli
    
    Created by:
    -----------
    Brent Nissens
    October 22, 2025
    """
    
    # Reference SAM parameters (RAM stimuli are calibrated with reference of SAM tone)
    level = 70                    # Stimulation intensity [dBSPL]
    p0 = 20e-6                    # Reference sound pressure level [Pa]
    fSAM = 4000                   # SAM carrier frequency [Hz]
    
    # RAM parameters
    fMod = 110                    # modulation frequency [Hz]
    md = 1                        # modulation depth
    phi = 3 * np.pi / 2           # Carrier phase in [rad]
    dur = 0.4                     # Stimulus duration in [s]
    
    # Get the time vector
    n = round(dur * fs)           # number of samples
    tVec = np.arange(n) / fs      # time vector
    
    # Generate the reference SAM tone
    carrierSAM = np.sin(2 * np.pi * fSAM * tVec)
    modSAM = np.sin(2 * np.pi * fMod * tVec + phi)
    SAM = carrierSAM * (1 + md * modSAM)
    
    # Calculate RMS for normalization
    SAM_rms = np.sqrt(np.mean(SAM**2))
    SAM = p0 * 10**(level/20) * SAM / SAM_rms  # adjust level such that 1Pa is the full-scale signal
    
    # RAM stimuli, for modulator 110Hz and carrier frequencies in fRAM
    # Generate square wave modulator (25% duty cycle) - matches MATLAB square(x, 25)
    modRAM = signal.square(2 * np.pi * fMod * tVec + phi, duty=0.25)
    
    # Initialize output matrix
    stim = np.zeros((len(fRAM), n))
    
    # Generate RAM stimuli for each carrier frequency
    for i in range(len(fRAM)):
        carrierRAM = np.sin(2 * np.pi * fRAM[i] * tVec)
        RAM = carrierRAM * (1 + md * modRAM)
        
        # Scale the RAM signal such that it has the same max value as the SAM
        stim[i, :] = np.max(SAM) / np.max(RAM) * RAM
    
    return stim


if __name__ == "__main__":
    # Example usage
    fs = 44100  # Sampling frequency
    fRAM = np.array([4000])  # Carrier frequencies
    
    stim = get_RAM_stims(fs, fRAM)
    print(f"Generated {stim.shape[0]} RAM stimuli with {stim.shape[1]} samples each")
    print(f"Sampling frequency: {fs} Hz")
    print(f"Carrier frequencies: {fRAM} Hz")