import numpy as np
from scipy import signal

M1 = 4.2767e-14
M3 = 5.1435e-14
M5 = 13.3093e-14


def get_bilinear_filter(fs, tau):
    '''
    Calculate filter using bilinear transform
    '''
    m = 2 * tau *fs
    a = (m-1) / (m+1)
    b = 1.0/(m+1)**2 * np.array([1,2,1])
    a = np.array([1, -2*a, a**2])
    return b, a


def delay(x, delay, fs):
    # Index where inhibition delay begins
    i_start = int(round(delay * fs))
    i_end = x.shape[-1] - i_start

    # This is a much faster appraoch than using np.zeros_like as verified by
    # profiling.
    delayed = np.empty_like(x)
    delayed[..., i_start:] = x[..., :i_end]
    delayed[..., :i_start] = 0
    return delayed


def cochlearNuclei(summed_an, fs, Acn=1.5, Scn=0.6, Tex=0.5e-3, Tin=2e-3,
                   inhibition_delay=1e-3):
    delayed_inhibition = delay(summed_an, inhibition_delay, fs)
    bEx, aEx = get_bilinear_filter(fs, Tex)
    bIn, aIn = get_bilinear_filter(fs, Tin)
    anEx = signal.lfilter(bEx, aEx, summed_an, axis=-1)
    anIn = signal.lfilter(bIn, aIn, delayed_inhibition, axis=-1)
    return Acn * (anEx - Scn * anIn)


def inferiorColliculus(cn, fs, Aic=1, Sic=1.5, Tex=0.5e-3, Tin=2e-3,
                       inhibition_delay=2e-3):
    delayed_inhibition = delay(cn, inhibition_delay, fs)
    bEx, aEx = get_bilinear_filter(fs, Tex)
    bIn, aIn = get_bilinear_filter(fs, Tin)
    cnEx = signal.lfilter(bEx, aEx, cn, axis=-1)
    cnIn = signal.lfilter(bIn, aIn, delayed_inhibition, axis=-1)
    return Aic * (cnEx - Sic * cnIn)
