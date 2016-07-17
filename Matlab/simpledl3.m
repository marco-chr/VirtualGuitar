clear all
close all
clc

fs=44100;               % sampling frequency
T=1/fs;                 % sampling time

duration=5;             % duration in sec of synthetized sound
t=zeros(1,fs*duration); % pickup output vector init

%%% H

g=0.9958;
a1=-0.4213;
b0=g*(1+a1);

f1=110.7;                 % string fundamental frequency
L=floor(fs/f1);         % approximated delay line length

omega1=2*pi*f1;
H=b0/(1+a1*exp(-1i*omega1*T));
phase_a=angle(H)/(omega1*T);        % phase for 1 pole LP filter

phase_c=(fs/f1)-L-phase_a;          % phase for AP filter

C=(sin((omega1*T-omega1*T*phase_c)/(2))/(sin(omega1*T+omega1*T*phase_c)/(2)));          % parameter C for AP filter

b_coeff=[1 a1+C a1*C];
a_coeff=[1 a1+C a1*C zeros(1,(L-3)) -(C*(g+g*a1)) -(g+g*a1)];

[Y,FS,NBITS]=wavread('gtrbody3.wav');     % input waveform

x=[Y' zeros(1,length(t)-length(Y))];

h1=filter(b_coeff,a_coeff,x);

%%% V

g=0.99;
a1=-0.4213;
b0=g*(1+a1);

f1=110.5;                % string fundamental frequency
L=floor(fs/f1);          % approximated delay line length

omega1=2*pi*f1;
H=b0/(1+a1*exp(-1i*omega1*T));
phase_a=angle(H)/(omega1*T);        % phase for 1 pole LP filter

phase_c=(fs/f1)-L-phase_a;          % phase for AP filter

C=(sin((omega1*T-omega1*T*phase_c)/(2))/(sin(omega1*T+omega1*T*phase_c)/(2)));          % parameter C for AP filter

b_coeff=[1 a1+C a1*C];
a_coeff=[1 a1+C a1*C zeros(1,(L-3)) -(C*(g+g*a1)) -(g+g*a1)];

[Y,FS,NBITS]=wavread('gtrbody3.wav');     % input waveform

x=[Y' zeros(1,length(t)-length(Y))];

h2=filter(b_coeff,a_coeff,x);


h=0.5*h1+0.5*h2;
figure(1)
plot(h)

sound(h,fs)

