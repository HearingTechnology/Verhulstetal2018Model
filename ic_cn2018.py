#from numpy import *
import numpy as np
from scipy import signal

M1=6.2755e-14
M3=7.2161e-14
M5=3.5200e-20

def cochlearNuclei(anfH,anfM,anfL,numH,numM,numL,fs):
    size=len(anfH[0,:])
    

    TF=19; #total no of fibers on each IHC
#    HSnormal=13;
#    MSnormal=3;
#    LSnormal=3;

    Acn=1.5;
    Scn=0.6;
    inhibition_delay=int(round(1e-3*fs));
    Tex=0.5e-3;
    Tin=2e-3;

    summedAN=numL*anfL+numM*anfM+numH*anfH;

    summedAN=summedAN;
    delayed_inhibition=np.zeros_like(summedAN)
    delayed_inhibition[inhibition_delay:,:]=summedAN[0:len(summedAN)-inhibition_delay,:]

#    print('here OK')
    #filters obtained with bilinear transform
    bEx = np.array([1, 2, 1]); #numerator
    aEx = np.array([1, -2*(1-1/Tex/fs/2)/(1+1/Tex/fs/2), (1-1/Tex/fs/2)**2/(1+1/Tex/fs/2)**2]); #denominator
    bIncn = np.array([1, 2, 1]); #numerator
    aIncn = np.array([1, -2*(1-1/Tin/fs/2)/(1+1/Tin/fs/2), (1-1/Tin/fs/2)**2/(1+1/Tin/fs/2)**2]) #denominator
    cn=Acn*(1/(Tex**2)*signal.lfilter(bEx,aEx,summedAN,axis=0)-Scn*1/(Tin**2)*signal.lfilter(bIncn,aIncn,delayed_inhibition,axis=0))
    cn=1.0/(4*fs**2)*cn;
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
    bEx = [1, 2, 1]; #numerator
    aEx = [1, -2*(1-1/Tex/fs/2)/(1+1/Tex/fs/2), (1-1/Tex/fs/2)**2/(1+1/Tex/fs/2)**2]; #denominator
    bIncn = [1, 2, 1]; #numerator
    aIncn = [1, -2*(1-1/Tin/fs/2)/(1+1/Tin/fs/2), (1-1/Tin/fs/2)**2/(1+1/Tin/fs/2)**2] #denominator
    ic=Aic*(1/(Tex**2)**2*signal.lfilter(bEx,aEx,cn,axis=0)-Sic/(Tin**2)**2*signal.lfilter(bIncn,aIncn,delayed_inhibition,axis=0))
    ic=ic
    return 1.0/(4*fs**2)*ic


