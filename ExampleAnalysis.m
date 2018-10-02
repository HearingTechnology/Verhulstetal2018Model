close all;
clear all;

%this file runs plots the results of a click simulation performed with
%ExampleSimulation.m
load Simulations.mat

%are the different conditions
L=[0 10 20 30 40 50 60 70 80 90 100];
p0=2e-5;
CF=output.cf;
figure,plot(CF),xlabel('Cochlear Channel Number [-]'),ylabel('Characteristic Frequency [Hz]')

fs_c=output.fs_bm; %the sampling frequency of BM, OAE and IHC are higher to avoid numerical errors (see Altoe et al., 2014 JASA)
fs=output.fs_an; %the sampling frequency of AN, CN, IC and waves are 5 times lower
t_c=[0:size(output(1).v,1)-1]/fs_c;
t=[0:size(output(1).anfH,1)-1]/fs;
f_c=0:fs_c/size(output(1).v,1):fs_c-1/size(output(1).v,1);
f=0:fs/size(output(1).anfH,1):fs-1/size(output(1).anfH,1);

% pick a channel number to plot results from. The CF corresponding to the
% channel number depends on whether you chose 'all', 'half' or 'abr'
No=245;

%reorganisation of the data for easier processing
for n=1:numel(L)
   oae(:,n)=output(n).e;
   vrms(:,n)=rms(output(n).v);
   ihcrms(:,n)=rms(output(n).ihc);
   v(:,n)=output(n).v(:,No);
   ihc(:,n)=output(n).ihc(:,No);
end

%% the OAE
%this figure is the ear-canal pressure: 
%to get the reflection-source OAE, do the following simulation:
%OAE_{reflections on}-OAE_{reflections off}
%to get the distortion-source OAE, do a normalisation using a low level
%(linear) simulation that the reflections off
figure,
subplot(2,1,1),plot(1000*t_c,oae)
xlabel('Time [ms]'),ylabel('Ear Canal Pressure [Pa]'),xlim([0 20]),ylim([-0.02 0.02]),legend(num2str(L')),legend('boxoff')
oae(1:200,:)=0; %zeropadding to plot OAE spectrum
subplot(2,1,2),plot(f_c/1000,20*log10(abs(fft(oae/p0))));
xlabel('Frequency [kHz]'),ylabel('EC Magnitude [dB re p0]'),xlim([0 12]),legend(num2str(L')),legend('boxoff')

%% v_bm and V_IHC 
figure,
subplot(2,2,1),plot(1000*t_c,v),
xlabel('Time [ms]'),ylabel('v_{bm} [m/s]'),xlim([0 30]),legend(num2str(L')),legend('boxoff')
title(['CF of ',num2str(round(CF(No))),' Hz'])
subplot(2,2,2),plot(CF/1000,20*log10(vrms)),
xlabel('CF [kHz]'),ylabel('rms of v_{bm} [dB re 1 m/s]'),xlim([0 14]),ylim([max(max(20*log10(vrms)))-100 max(max(20*log10(vrms)))+10])
title(['Excitation Pattern'])
subplot(2,2,3),plot(1000*t_c,ihc),
xlabel('Time [ms]'),ylabel('V_{ihc} [V]'),xlim([0 30]),legend(num2str(L')),legend('boxoff')
title(['CF of ',num2str(round(CF(No))),' Hz'])
subplot(2,2,4),plot(CF/1000,ihcrms),
xlabel('CF [kHz]'),ylabel('rms of V_{ihc} [V]'),xlim([0 14]),%ylim([max(max(20*log10(ihcrms)))-100 max(max(20*log10(ihcrms)))+10])
title(['Excitation Pattern'])

%reorganisation of the data for easier processing
for n=1:numel(L)
   HSR(:,n)=output(n).anfH(:,No);
   MSR(:,n)=output(n).anfM(:,No);
   LSR(:,n)=output(n).anfL(:,No);
   AN(:,n)=output(n).an_summed(:,No);
   CN(:,n)=output(n).cn(:,No);
   IC(:,n)=output(n).ic(:,No);
   W1(:,n)=output(n).w1;
   W3(:,n)=output(n).w3;
   W5(:,n)=output(n).w5;
   EFR(:,n)=output(n).w1+output(n).w3+output(n).w5;
end

%single unit responses
figure,
subplot(3,2,1),plot(1000*t,HSR)
title(['CF of ',num2str(round(CF(No))),' Hz'])
xlim([0 20]),xlabel('Time [ms]'),ylabel('HSR fiber [spikes/s]')
legend(num2str(L')),legend('boxoff')
subplot(3,2,3),plot(1000*t,MSR)
xlim([0 20]),xlabel('Time [ms]'),ylabel('MSR fiber [spikes/s]')
subplot(3,2,5),plot(1000*t,LSR)
xlim([0 20]),xlabel('Time [ms]'),ylabel('LSR fiber [spikes/s]')

subplot(3,2,2),plot(1000*t,AN)
title(['CF of ',num2str(round(CF(No))),' Hz'])
xlim([0 20]),xlabel('Time [ms]'),ylabel('sum AN [spikes/s]')
%Spikes summed across all fibers @ 1 CF
subplot(3,2,4),plot(1000*t,CN)
xlim([0 20]),xlabel('Time [ms]'),ylabel('CN [spikes/s]')
subplot(3,2,6),plot(1000*t,IC)
xlim([0 20]),xlabel('Time [ms]'),ylabel('IC [spikes/s]')

%population responses
figure
subplot(4,1,1),plot(1000*t,1e6*W1)
title(['Population Responses summed across simulated CFs'])
xlim([0 20]),xlabel('Time [ms]'),ylabel('W-1 [\muV]')
legend(num2str(L')),legend('boxoff')
subplot(4,1,2),plot(1000*t,1e6*W3)
xlim([0 20]),xlabel('Time [ms]'),ylabel('W-3 [\muV]')
subplot(4,1,3),plot(1000*t,1e6*W5)
xlim([0 20]),xlabel('Time [ms]'),ylabel('W-5 [\muV]')
subplot(4,1,4),plot(1000*t,1e6*EFR)
xlim([0 20]),xlabel('Time [ms]'),ylabel('EFR [\muV]')

%The model code and interface was written by Alessandro Altoè and Sarah Verhulst (copyright 2012,2014,2015,2016,2018) 
%and is licensed under the UGent acadamic license (see details in license file that is part of this repository). 
%The Verhulstetal2018Model consists of the following files: 
%tridiag.so, cochlea_utils.c, run_model2018.py, model2018.m, cochlear_model2017.py, inner_hair_cell2018.py, auditory_nerve2017.py, ic_cn2017.py, ExampleSimulation.m, ExampleAnalysis.m, the HI profiles in the Poles folder. 
 


