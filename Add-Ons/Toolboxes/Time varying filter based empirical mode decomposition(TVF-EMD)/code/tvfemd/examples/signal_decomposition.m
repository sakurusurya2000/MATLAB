%---------- multi-component signal decomposition using EMD and TVF-EMD
clc;
clear;
close all;
Fs = 1000;   %sampling rate
freq1 = 50; %frequency
w1=2*pi*freq1; %rad frequency
t = 0:1/Fs:4; %time span

%---------- multicomponent signal 
sig1=1*cos(1*w1*t); % linear signal
sig2=1*cos(0.7*w1*t);
sig3=1*cos(0.3*w1*t);

FM=1*cos((0.2+0.04*t).*w1.*t);% chirp signal
FM2=cos((0.04+0.04*t).*w1.*t);

x1=1*sig1+ 1*sig2+1*sig3+0.2*randn(size(sig1)); % noisy linear signal
x2=1*FM+1*FM2+0.1*randn(size(sig1));   % noisy non-stationary signal

 imf1=tvf_emd(x1);% decompose linear signal using TVF-EMD

imf2=tvf_emd(x2); % decompose chirp signal using TVF-EMD