from numpy import *
from scipy import signal
#parameters
Cm=12.5e-12;
Gmet=30e-9;
s1=16e-9;
s0=s1*3;
x0=20e-9;
#x0=34e-9
tauMet=50e-6;
Gleak=00.0e-9;
EP=90e-3;
Ekf=-71e-3;
Eks=-78e-3;
Gk=230e-9;
xk=-31e-3;
sk=10.5e-3;
tkf1=0.3e-3;
tkf2=0.1e-3;
tks1=8e-3;
tks2=2e-3;
GcaM=4.1e-9;
tauCa=0.25e-3;
xca=-30.0e-3;
sca=7.5e-3;
Eca=45e-3;
tauInS=0.5;
tauInF=50e-3;
CaInMax=0.4;
xCaIn=-43e-3;
sCaIn=6e-3;
Ileak=0.0e-9;
resting_potential=-0.05703; #resting_potential at equilibrium
peak_potential=-0.04; #peak resting potential at 100 dB 4 kHz (where nerve fibers saturate)

def inner_hair_cell_potential(mu,fs):
    size=len(mu[0,:])
    fs=fs
    dt=1.0/fs
    alphaMet=exp(-dt/tauMet);
    a=zeros([size,2]);
    mt=zeros(size)
    alphakf1=exp(-dt/tkf1);
    alphaks1=exp(-dt/tks1);
    Vm=zeros([1,size])+resting_potential
    mt=zeros([1,size])+1/(1+exp(x0/s0)*(1+exp(x0/s1)));
#    mt=zeros([1,size])+1/(1+exp(x0/s1));
    Vsol=zeros_like(mu,dtype=float)
    mtIn=1/(1+exp((x0-mu)/s0)*(1+exp((x0-mu)/s1)));
    mkf1=1/(1+exp(-(Vm-xk)/sk));
    mks1=1/(1+exp(-(Vm-xk)/sk));
    zero_time=int(fs*50e-3);
    for i in range(zero_time):
        Imet=(Gmet*mt)*(Vm-EP);
        Ileak=Gleak*(Vm-Eca);
        mk=1/(1+exp(-(Vm-xk)/sk));
        mkf1=(1-alphakf1)*mk+alphakf1*mkf1;
        mks1=(1-alphaks1)*mk+alphaks1*mks1;
        Gkf=Gk*mkf1;
        Gks=Gk*mks1;
        Ikf=Gkf*(Vm-Ekf);
        Iks=Gks*(Vm-Eks);
        dV=-(Ileak+Imet+Ikf+Iks)/Cm;
        Vm=Vm+dV*dt;
    for i in range(len(mtIn[:,0])):
        mt=(1-alphaMet)*mtIn[i,:]+alphaMet*mt
        Imet=(Gmet*mt)*(Vm-EP);
        Ileak=Gleak*(Vm-Eca);
        mk=1/(1+exp(-(Vm-xk)/sk));
        mkf1=(1-alphakf1)*mk+alphakf1*mkf1;
        mks1=(1-alphaks1)*mk+alphaks1*mks1;
        Gkf=Gk*mkf1;
        Gks=Gk*mks1;
        Ikf=Gkf*(Vm-Ekf);
        Iks=Gks*(Vm-Eks);
        dV=-(Ileak+Imet+Ikf+Iks)/Cm;
        Vm=Vm+dV*dt;
        Vsol[i,:]=Vm
    # print(Vsol)
    return Vsol






#class ihc():
#    def __init__(self,cf,mu,fs):
#        size=len(cf)
#        self.cf=cf
#        self.fs=fs
#        dt=1/fs
#        self.alphaMet=exp(-dt/tauMet);
#        self.a=zeros([size,2]);
#        self.mt=zeros(size)
#        self.alphakf1=exp(-dt/tkf1);
#        self.alphaks1=exp(-dt/tks1);
#        self.Vm=zeros(size)+resting_potential
#        self.mt=zeros(size)+1/(1+exp(x0/s0)*(1+exp(x0/s1)));
#        self.Vsol=np.zeros_like(mu)
#        self.mtIn=1/(1+exp((x0-mu)/s0)*(exp((x0-mu)/s1)+1));
#        self.mkf=1/(1+exp(-(self.Vm-xk)/sk));
#        self.mks=1/(1+exp(-(self.Vm-xk)/sk));
#    
#    def solve_step(self):
#        fs=self.fs
#        dt=1/fs
#        zero_time=fs*200e-3;
#        for i in range(tpad):
#            Imet=(Gmet*mt)*(V-EP);
#            Ileak=Gleak*(V-Eca);
#            self.mk=1/(1+exp(-(self.Vm-xk)/sk));
#            self.mkf1=(1-self.alphakf1)*self.mk+self.alphakf1*self.mkf1;
#            self.mks1=(1-self.alphaks1)*self.mk+self.alphaks1*self.mks1;
#            Gkf=Gk*self.mkf1;
#            Gks=Gk*self.mks1;
#            Ikf=Gkf*(self.Vm-Ekf);
#            Iks=Gks*(self.Vm-Eks);
#            dV=-(Ileak+Imet+Ikf+Iks)/Cm;
#            self.Vm=self.Vm+dV*dt;
#        for i in range(len(self.mtIn[:,0])):
#            mt=(1-self.alphaMet)*self.mtIn[:,i]+self.alphaMet*self.mt
#            Imet=(Gmet*mt)*(V-EP);
#            Ileak=Gleak*(V-Eca);
#            self.mk=1/(1+exp(-(self.Vm-xk)/sk));
#            self.mkf1=(1-self.alphakf1)*self.mk+self.alphakf1*self.mkf1;
#            self.mks1=(1-self.alphaks1)*self.mk+self.alphaks1*self.mks1;
#            Gkf=Gk*self.mkf1;
#            Gks=Gk*self.mks1;
#            Ikf=Gkf*(self.Vm-Ekf);
#            Iks=Gks*(self.Vm-Eks);
#            dV=-(Ileak+Imet+Ikf+Iks)/Cm;
#            self.Vm=self.Vm+dV*dt;
#            self.Vsol[:,i]=self.Vm
#
#
#
#
#
#
#
#
