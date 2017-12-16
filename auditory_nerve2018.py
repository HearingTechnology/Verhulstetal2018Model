from numpy import *
from inner_hair_cell2018 import resting_potential, peak_potential

#parameters
r1=220 #reserve pool max. replenishment rate
x=700 #RRP replenishment rate (Pangrisc 2010,Chapocnikov 2014)
M=14 # Max. vesicles in the ready release pool or release sites (reasonable, ref?)
M2=60 # Max. vesicles in the second pool, fitted parameter
# with the paramaters is 250 release/s, because of refractoriness the steady state spike rate goes around 200 spike/s
refr_tail=0.6e-3 # relative refreactory period (from Peterson and Heil 2014)
abs_refractoriness=0.6e-3 # absolute refractory period (from Peterson and Heil 2014)
ss=1.5e-3 #sensitivity of release rate on IHC potential, fitted paramter
tCa=0.2e-3 #time constant of the Ca2 channel
#driven exocytosis rate is a boltzmann version of Vm, low-pass filter with the time constant of the Ca2+ channels. This is a shortcut to tune the nonlinear relationship between Vm and exocytosis rate based on AF data. Ca2+ signaling varies substantialy between synapses of the same IHC (Frank 2009,Ohn 2016). We use a shortcut Vm->Exocytosis rate, otherwise Vm->Ca2+ (at individual synapse)->exocytosis rate would require the tuning of too many parameters (here one free parameter ss, and 2 fitted parameters sp,psr)

def auditory_nerve_fiber(Vm,fs,spont):
    size=len(Vm[0,:])
    dt=1.0/fs
    if(spont==0):
        sp=1 #spont rate
        psr=800 #peak spike rate, based upon  Taberner and Liberman 2005
    elif(spont==1):
        sp=10.0
        psr=1.0e3
    else:
        sp=68.5
        psr=3.0e3
    #parameters for vesicles trafficking
#    psr=2.7e3
#    psr=psr*4;
    r=r1/M2 # replenishment rate per vesicle location
    M2_steady=M2-sp/r #steady state values
    Msteady=M*(M2_steady/M2-sp/x)
    wt=M2_steady+zeros([1,size]) #reserve pool
    qt=Msteady+zeros([1,size]) #RRP
    xdt=x*dt
    rdt=r*dt
    #parameters for refractoriness
    alpha_ref=exp(-dt/refr_tail)
    relRefFraction=zeros([1,size]) #how much the firing probability decreases due to relative refractoriness
    available=1.0+zeros([1,size]) #number of fibers not in a refractory state
    buf_lgt=int(round(abs_refractoriness*fs)) #length of buffer to store the firing history (in order to account for absolute refractoriness)
    buf_ref=zeros([buf_lgt,size],dtype=float) #buffer to store the history of firing
    buf_pointer=0
    pp=psr/Msteady#/(1/(1+exp(-(peak_potential-vh)/ss))) #multiplier relating the activation nonlinearity (between 0 and 1) with actual firing rate
    #parameters to relate Vm with firing rate
    rat=log((psr-sp)/sp) #find half activation voltage of exocytosis rate based on peak and spontaneous rate
    vh=rat*ss+resting_potential
    k0=sqrt(1/(1+exp(-(resting_potential-vh)/ss)))+zeros([1,size])  #driven exocytosis rate at rest
    #take the square root, filter it with a first order filter and then square it. This is equivalent to a second order activation of the ion channels
    alphaCa=exp(-dt/tCa)
 

    zero_time=int(50e-3*fs)
    for i in range(zero_time):
        vesicleReleaseRate=pp*(k0**2)
        releaseProb=vesicleReleaseRate*dt
        qt[qt>M]=M
        wt[wt>M2]=M2
        ejected=releaseProb*qt
        replenishReserve= rdt*(M2-wt)
        replenishRRP=xdt*fmax(wt/M2-qt/M,0)
        qt= qt + replenishRRP - ejected
        wt= wt - replenishRRP + replenishReserve
        firing=(available-relRefFraction)*ejected
        recovered=buf_ref[buf_pointer,:]
        relRefFraction=alpha_ref*relRefFraction+(1-alpha_ref)*recovered
        available=available-firing+recovered
        buf_ref[buf_pointer,:]=firing
        buf_pointer=mod(buf_pointer+1,buf_lgt)
    k=k0
    kin=sqrt(1/(1+exp(-(Vm-vh)/ss)))
    solution=zeros_like(Vm)
    for i in range(len(kin[:,0])):
        k=alphaCa*k+(1-alphaCa)*kin[i,:]
        vesicleReleaseRate=pp*(k**2)
        releaseProb=vesicleReleaseRate*dt
        qt[qt<0]=0
        wt[wt<0]=0
        qt[qt>M]=M
        wt[wt>M2]=M2
        ejected=releaseProb*qt
        replenishReserve= rdt*(M2-wt)
        replenishRRP=xdt*fmax(wt/M2-qt/M,0)
        qt= qt + replenishRRP - ejected
        wt= wt - replenishRRP + replenishReserve
        firing=(available-relRefFraction)*ejected
        recovered=buf_ref[buf_pointer,:]
        relRefFraction=alpha_ref*relRefFraction+(1-alpha_ref)*recovered
        available=available-firing+recovered
        buf_ref[buf_pointer,:]=firing
        buf_pointer=mod(buf_pointer+1,buf_lgt)
        solution[i,:]=firing
    return solution





