close all;
clear all;

%this file runs an example simulation for a click stimulus of different
%stimulus levels

L=[0 10 20 30 40 50 60 70 80 90 100];
fs=100e3;
p0=2e-5;
dur=50e-3;
t=(0:1/fs:dur);
click_duration=10; % 100 us click
stim=zeros(numel(L),length(t)); %the simulation runs as long as the stimulus

%depending on the number of cores you have, you can run more or less stimuli in one simulation
%the dimensions are: conditions x duration (the code takes care of it in case you do the inverse)
for j=1:numel(L) 
    stim(j,10:10+click_duration)=p0*10^(L(j)/20)*sqrt(2)*2;
end
%to set the amplitude of the stimulus, the digital 1 is multiplied with a
%Pascal value:
%   for rms: make a stimulus with rms of 1: multiply with p0*10^(L_{SPLindB}/20);
%   for pure tone: multiply A with sqrt(2)*p0*10^(L_{SPLindB}/20)
%   for condensation click in peSPL: multiply with
%   2*sqrt(2)*p0*10^(L_{SPLindB}/20);

%% load the starting poles (alpha*,30) for NH or HI models
sheraP=load('StartingPoles.dat');
%for a hearing impaired model, load other startingpoles e.g.:
%sheraP=load('./Poles/Flat00_Slope30/StartingPoles.dat');

%% set synaptopathy and store the results
% decide how many channels you want to store in the simulations ('all','half','abr')
% decide how many fibers you want to include to compute the CN and IC responses
% default: 13 HSR (70 spikes/s), 3 MSR (10) and 3 LSR (1) fibers (across all CF-channels) , you want
% to set CF-dependent synaptopathy profiles, insert a vector of size:   
    %1000x1 for 'all'
    %500x1 for 'half'
    %401x1 for 'abr'
% the indices should correspond to the CFs simulated for that condition (output(1).cf)
% run the model ones using default settings to write out the CFs, and then
% determine your CF dependent synaptopathy profiles.

% decide which responses you want to store: 
    % e= emission 
    % v= velocity 
    % i=ihc 
    % h=hsr 
    % m=msr 
    % l=lsr 
    % b=summed AN responses for each CF as well as CN and IC responses 
    % w=population response waves I, III, V
    
% to store all channels (1000)
 output=model2018(stim,fs,'all',1,'evihmlbw',1,sheraP,0.05,'vel',13,3,3,1,[pwd(),'/']);
% to store half of the channels (500)
% output=model2018(stim,fs,'half',1,'evihmlbw',1,sheraP,0.05,'vel',13*ones(500,1),3,3,1,[pwd(),'/']);
%store only the ABR channels: CFs between 112 Hz and 12 kHz
%output=model2018(stim,fs,'abr',1,'evihmlbw',1,sheraP,0.05,'vel',13,3,3,1,[pwd(),'/']);

save('Simulations.mat','output','-v7.3')




%The model code and interface was written by Alessandro Altoè and Sarah Verhulst (copyright 2012,2014,2015,2016,2018) 
%and is licensed under the UGent acadamic license (see details in license file that is part of this repository). 
%The Verhulstetal2018Model consists of the following files: 
%tridiag.so, cochlea_utils.c, run_model2018.py, model2018.m, cochlear_model2017.py, inner_hair_cell2018.py, auditory_nerve2017.py, ic_cn2017.py, ExampleSimulation.m, ExampleAnalysis.m, the HI profiles in the Poles folder. 
 
