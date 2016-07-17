close all
clear all
clc

% read ir
[Y,FS,NBITS]=wavread('gtrbody.wav');
Y=Y';
dt=1/FS;
t=[0:dt:(length(Y)/FS)-dt];

figure(1)
plot(t,Y)

t1=[0:dt:0.03];
excite1=(0.25/0.03)*t1;
excite=[excite1 zeros(1,length(t)-length(excite1))];

figure(2)
plot(t,excite)

agg=conv(Y,excite);

figure(3)
plot(agg)

sound(agg,22050)
sound(Y,22050)

WAVWRITE(agg,FS,'gtrbody3.wav')