"""
Python equivalent of the MATLAB model2018.m function.

Computational model of the auditory periphery (Verhulst, Altoe, Vasilkov, 2018).

Author: Converted from MATLAB by Brent Nissens
Based on original work by Alessandro Altoe and Sarah Verhulst
"""

import numpy as np
import scipy as sp
import scipy.signal
import multiprocessing as mp
import warnings
from typing import Union, Optional, Dict, List, Any
import os
from get_RAM_stims import get_RAM_stims
import scipy.io as sio

# Import the model components
from cochlear_model2018 import cochlea_model
import inner_hair_cell2018 as ihc
import auditory_nerve2018 as anf
import ic_cn2018 as nuclei

# Suppress warnings like in the original run_model2018.py
warnings.filterwarnings("ignore")


class ModelOutput:
    """Class to store model outputs similar to MATLAB struct."""
    
    def __init__(self):
        self.cf = None
        self.v = None           # BM velocity
        self.y = None           # BM displacement  
        self.emission = None    # pressure output from middle ear
        self.ihc = None         # IHC receptor potential
        self.anfH = None        # HSR fiber spike probability
        self.anfM = None        # MSR fiber spike probability
        self.anfL = None        # LSR fiber spike probability
        self.an_summed = None   # summation of HSR, MSR and LSR per channel
        self.cn = None          # cochlear nuclei output
        self.ic = None          # IC output
        self.w1 = None          # wave 1
        self.w3 = None          # wave 3
        self.w5 = None          # wave 5
        self.fs_bm = None       # sampling frequency of BM simulations
        self.fs_ihc = None      # sample rate of IHC output
        self.fs_an = None       # sample rate of AN output
        self.fs_abr = None      # sample rate of IC, CN and W1/3/5 outputs


def model2018(
    sign: np.ndarray,
    fs: float,
    fc: Union[np.ndarray, str, int, float] = 'all',
    irregularities: Union[int, np.ndarray] = 1,
    storeflag: str = 'vihlmeb',
    subject: int = 1,
    sheraPo: Union[float, np.ndarray] = 0.06,
    IrrPct: float = 0.05,
    non_linear_type: str = 'vel',
    nH: Union[int, np.ndarray] = 13,
    nM: Union[int, np.ndarray] = 3,
    nL: Union[int, np.ndarray] = 3,
    clean: int = 1,
    data_folder: str = './') -> List[ModelOutput]:
    """
    Computational model of the auditory periphery (Verhulst, Altoe, Vasilkov, 2018).
    
    Parameters:
    -----------
    sign : np.ndarray
        Stimulus signal
    fs : float
        Sample rate
    fc : Union[np.ndarray, str, int, float], optional
        Probe frequency or alternatively a string 'all' to probe all cochlear
        sections or 'half' to probe half sections, 'abr' to store the
        401 sections used to compute the abr responses in Verhulst et al. 2017
    irregularities : Union[int, np.ndarray], optional
        Decide whether turn on (1) or off (0) irregularities and
        nonlinearities of the cochlear model (default 1)
    storeflag : str, optional
        String that sets what variables to store from the computation, 
        each letter correspond to one desired output variable (e.g., 'avhl' 
        to store acceleration, displacement, high and low spont. rate fibers.) 
        Default: 'vihlmeb'.
    subject : int, optional
        Number representing the seed to generate the random
        irregularities in the cochlear sections (default 1)
    sheraPo : Union[float, np.ndarray], optional
        Starting real part of the poles of the cochlear model
        it can be either an array with one value per BM section, or a
        single value for all sections (default 0.06)
    IrrPct : float, optional
        Magnitude of random perturbations on the BM (irregularities, default 0.05=5%)
    non_linear_type : str, optional
        Select the type of nonlinearity in the BM model.
        Currently implemented:
           'vel'= instantaneous nonlinearity based on local BM velocity (see Verhulst et al. 2012)
           'none'= linear model
    nH, nM, nL : Union[int, np.ndarray], optional
        Number of high, medium and low spont. fibers employed to
        compute the response of cn and ic nuclei. Default 13,3,3. These
        parameters can be passed either as a single value for all sections or
        as an array with each value corresponding to a single CF location
    clean : int, optional
        Not used in Python version (kept for compatibility)
    data_folder : str, optional
        Not used in Python version (kept for compatibility)
        
    Returns:
    --------
    List[ModelOutput]
        List of ModelOutput objects, one per channel, containing:
        - v: BM velocity (store 'v')
        - y: BM displacement (store 'y') 
        - emission: pressure output from the middle ear (store 'e')
        - cf: center frequencies (always stored)
        - ihc: IHC receptor potential (store 'i')
        - anfH: HSR fiber spike probability [0,1] (store 'h')
        - anfM: MSR fiber spike probability [0,1] (store 'm')
        - anfL: LSR fiber spike probability [0,1] (store 'l')
        - an_summed: summation of HSR, MSR and LSR per channel (storeflag 'b')
        - cn: cochlear nuclei output (storeflag 'b')
        - ic: IC (storeflag 'b')
        - w1, w3, w5: wave 1,3 and 5 (storeflag 'w')
        - fs_bm: sampling frequency of the bm simulations
        - fs_ihc: sample rate of the inner hair cell output
        - fs_an: sample rate of the an output
        - fs_abr: sample rate of the IC,CN and W1/3/5 outputs
    """
    
    DECIMATION = 5
    sectionsNo = 1000
    
    # Handle input signal dimensions
    sign = np.atleast_2d(sign)
    if sign.shape[0] > sign.shape[1]:
        sign = sign.T  # Make sure it's channels x samples
    
    channels = sign.shape[0]
    
    # Handle irregularities parameter
    if np.isscalar(irregularities):
        irregularities = irregularities * np.ones(channels)
    
    # Handle probe points (fc parameter)
    if isinstance(fc, str):
        if fc == 'all':
            l = sectionsNo
            probes = 'all'
        elif fc == 'half':
            l = sectionsNo // 2
            probes = 'half'
        elif fc == 'abr':
            l = 401
            probes = 'abr'
        else:
            raise ValueError(f"Unknown fc string option: {fc}")
    else:
        fc = np.atleast_1d(fc)
        l = len(fc)
        probes = np.round(fc).astype(int)
    
    # Handle sheraPo parameter
    if np.isscalar(sheraPo):
        sheraPo_val = sheraPo
    else:
        sheraPo_val = np.atleast_1d(sheraPo)
    
    # Handle fiber numbers
    if np.isscalar(nH):
        numH = nH
    else:
        numH = np.atleast_1d(nH)
        
    if np.isscalar(nM):
        numM = nM
    else:
        numM = np.atleast_1d(nM)
        
    if np.isscalar(nL):
        numL = nL
    else:
        numL = np.atleast_1d(nL)
    
    # Create cochlear models for each channel
    cochlear_list = []
    for i in range(channels):
        coch = cochlea_model()
        cochlear_list.append([coch, sign[i], irregularities[i], i])
    
    def solve_one_cochlea(model_data):
        """Process one cochlear channel."""
        coch, sig, irr_on, channel_idx = model_data
        
        # Initialize model
        coch.init_model(
            sig, fs, sectionsNo, probes,
            Zweig_irregularities=irr_on,
            sheraPo=sheraPo_val,
            subject=subject,
            IrrPct=IrrPct,
            non_linearity_type=non_linear_type
        )
        
        # Solve cochlear model
        coch.solve()
        
        # Create output structure
        output = ModelOutput()
        output.cf = coch.cf
        output.fs_bm = fs
        output.fs_ihc = fs
        output.fs_an = fs // DECIMATION
        output.fs_abr = fs // DECIMATION
        
        # Store BM velocity if requested
        if 'v' in storeflag:
            output.v = coch.Vsolution
            
        # Store BM displacement if requested  
        if 'y' in storeflag:
            output.y = coch.Ysolution
            
        # Store emissions if requested
        if 'e' in storeflag:
            output.emission = coch.oto_emission
            
        # Process IHC if needed
        if any(flag in storeflag for flag in 'ihmlbw'):
            magic_constant = 0.118
            Vm = ihc.inner_hair_cell_potential(coch.Vsolution * magic_constant, fs)
            
            if 'i' in storeflag:
                output.ihc = Vm
                
            # Resample for AN processing
            dec_factor = 5
            Vm_resampled = sp.signal.decimate(Vm, dec_factor, axis=0, n=30, ftype='fir')
            Vm_resampled[0:5, :] = Vm[0, 0]  # resting value to eliminate noise from decimate
            Fs_res = fs / dec_factor
            
            # Process auditory nerve fibers if needed
            if any(flag in storeflag for flag in 'hmlbw'):
                
                # High spontaneous rate fibers
                if 'h' in storeflag or 'b' in storeflag:
                    anfH = anf.auditory_nerve_fiber(Vm_resampled, Fs_res, 2) * Fs_res
                    if 'h' in storeflag:
                        output.anfH = anfH
                        
                # Medium spontaneous rate fibers  
                if 'm' in storeflag or 'b' in storeflag or 'w' in storeflag:
                    anfM = anf.auditory_nerve_fiber(Vm_resampled, Fs_res, 1) * Fs_res
                    if 'm' in storeflag:
                        output.anfM = anfM
                        
                # Low spontaneous rate fibers
                if 'l' in storeflag or 'b' in storeflag or 'w' in storeflag:
                    anfL = anf.auditory_nerve_fiber(Vm_resampled, Fs_res, 0) * Fs_res
                    if 'l' in storeflag:
                        output.anfL = anfL
                
                # Process brainstem nuclei if needed
                if 'b' in storeflag or 'w' in storeflag:
                    cn, anSummed = nuclei.cochlearNuclei(anfH, anfM, anfL, numH, numM, numL, Fs_res)
                    ic = nuclei.inferiorColliculus(cn, Fs_res)
                    
                    if 'b' in storeflag:
                        output.cn = cn
                        output.ic = ic
                        output.an_summed = anSummed
                        
                    # Process ABR waves if needed
                    if 'w' in storeflag:
                        output.w1 = nuclei.M1 * np.sum(anSummed, axis=1)
                        output.w3 = nuclei.M3 * np.sum(cn, axis=1)
                        output.w5 = nuclei.M5 * np.sum(ic, axis=1)
        
        return output
    
    # Process all channels (parallel processing)
    if channels == 1:
        # Single channel - no need for multiprocessing
        results = [solve_one_cochlea(cochlear_list[0])]
    else:
        # Multiple channels - use multiprocessing
        with mp.Pool(mp.cpu_count(), maxtasksperchild=1) as p:
            results = p.map(solve_one_cochlea, cochlear_list)
    
    print("cochlear simulation: done")
    
    return results


if __name__ == "__main__":
    # Create a simple test stimulus4000
    fs = 1e5
    stimulus = get_RAM_stims(fs,np.array([4000]))

    # Load the poles profile
    sheraP = np.loadtxt(f'./Poles/Flat00/StartingPoles.dat')

    # Run model
    print("Running model2018 example...")
    results = model2018(
        stimulus, 
        fs, 
        fc='abr',
        irregularities=0.05,
        storeflag='evihmlbw',
        subject=1,
        sheraPo=sheraP,
        IrrPct=0.05,
        non_linear_type='vel',
        nH=13,
        nM=3,
        nL=3,
        clean=1,
        data_folder='./'
    )
    
    # Display results
    output = results[0]
    print(f"Model completed successfully!")