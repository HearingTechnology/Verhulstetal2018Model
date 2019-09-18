import numpy as np
from scipy import signal

M1=4.2767e-14 
M3=5.1435e-14 
M5=13.3093e-14 

def cochlearNuclei(anfH,anfM,anfL,numH,numM,numL,fs):
    size=len(anfH[0,:])
    

    TF=19; #total no of fibers on each IHC
    # HSnormal=13;
    # MSnormal=3;
    # LSnormal=3;

    Acn=1.5;
    Scn=0.6;
    inhibition_delay=int(round(1e-3*fs));
    Tex=0.5e-3;
    Tin=2e-3;

    summedAN=numL*anfL+numM*anfM+numH*anfH;

    summedAN=summedAN;
    delayed_inhibition=np.zeros_like(summedAN)
    delayed_inhibition[inhibition_delay:,:]=summedAN[0:len(summedAN)-inhibition_delay,:]

    #filters obtained with bilinear transform
    # # Excitatory filter:
    m = (2*Tex*fs)
    a = (m-1)/(m+1)
    bEx = 1.0/(m+1)**2*np.array([1,2,1]); #numerator
    aEx = np.array([1, -2*a, a**2]); #denominator

    # # Inhibitory filter:
    m = (2*Tin*fs)
    a = (m-1)/(m+1)
    bIncn = 1.0/(m+1)**2*np.array([1, 2, 1]); # numerator
    aIncn = np.array([1, -2*a, a**2]); # denominator

    cn=Acn*(signal.lfilter(bEx,aEx,summedAN,axis=0)-Scn*signal.lfilter(bIncn,aIncn,delayed_inhibition,axis=0))
    return cn,summedAN
             
def inferiorColliculus(cn,fs):
    size=len(cn[0,:])
    Tex=0.5e-3;
    Tin=2e-3;
    Aic=1;
    Sic=1.5;
    inhibition_delay=int(round(2e-3*fs));

    delayed_inhibition=np.zeros_like(cn)
    delayed_inhibition[inhibition_delay:,:]=cn[0:len(cn)-inhibition_delay,:]

    # # Excitatory filter:
    m = (2*Tex*fs)
    a = (m-1)/(m+1)
    bEx=1.0/(m+1)**2*np.array([1,2,1]); #numerator
    aEx=np.array([1, -2*a, a**2]); #denominator

    # # Inhibitory filter:
    m = (2*Tin*fs)
    a = (m-1)/(m+1)
    bIncn = 1.0/(m+1)**2*np.array([1, 2, 1]); # numerator
    aIncn = np.array([1, -2*a, a**2]); # denominator

    ic=Aic*(signal.lfilter(bEx,aEx,cn,axis=0)-Sic*signal.lfilter(bIncn,aIncn,delayed_inhibition,axis=0))
    return ic
