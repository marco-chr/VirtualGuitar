clc
close all
clear all

fid=fopen('LP_filter.dat');
C = textscan(fid,'%f %f %f %f %f');
fclose(fid);

string = C{1};
fret = C{2};
freq = C{3};
g0 = C{4};
a1 = C{5};

for index=1:6

[g0m,g0b]=polyfit(fret(index*8-7:index*8),g0(index*8-7:index*8),2);
[a1m,a1b]=polyfit(fret(index*8-7:index*8),a1(index*8-7:index*8),2);

g0_fit=g0m(1)*fret(index*8-7:index*8).^2+g0m(2)*fret(index*8-7:index*8)+g0m(3);
a1_fit=a1m(1)*fret(index*8-7:index*8).^2+a1m(2)*fret(index*8-7:index*8)+a1m(3);

figure(2*index-1)
plot(fret(index*8-7:index*8),g0(index*8-7:index*8))
hold on
plot(fret(index*8-7:index*8),g0_fit,'--')
hold off
xlabel('Fret number');
ylabel('Loop filter coefficient g');

figure(2*index)

plot(fret(index*8-7:index*8),a1(index*8-7:index*8))
hold on
plot(fret(index*8-7:index*8),a1_fit,'--')
hold off
xlabel('Fret number');
ylabel('Loop filter coefficient a');

fid=fopen('coeff_fit.dat','a');
fprintf(fid,'%f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\r\n',index,g0m(1),g0m(2),g0m(3),a1m(1),a1m(2),a1m(3));
fclose(fid);

end

