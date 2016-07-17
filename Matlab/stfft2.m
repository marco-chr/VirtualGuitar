clear all
close all
clc

% opening file dialog
[FileName,PathName,FilterIndex] = uigetfile('*.wav','Select wav sample to be analyzed');

% note frequencies vectors
E2=[82.51 87.31 92.50 98.00 103.8 110.0 116.5 123.5 130.8 138.6 146.8 155.6 164.8 174.5 185.0 196.0 207.7 220.0 233.1 246.9 261.6];
A2=[110.0 116.5 123.5 130.8 138.6 146.8 155.6 164.8 174.6 185.0 196.0 207.7 220.0 233.1 246.9 261.6 277.2 293.7 311.1 329.6 349.2];
D3=[146.8 155.6 164.8 174.6 185.0 196.0 207.7 220.0 233.1 246.9 261.6 277.2 293.7 311.1 329.6 349.2 370.0 392.0 415.3 440.0 466.2];
G3=[196.0 207.7 220.0 233.1 246.9 261.6 277.2 293.7 311.1 329.6 349.2 370.0 392.0 415.3 440.0 466.2 493.9 523.3 554.4 587.3 622.3];
B3=[246.9 261.6 277.2 293.7 311.1 329.6 349.2 370.0 392.0 415.3 440.0 466.2 493.9 523.3 554.4 587.3 622.3 659.3 698.5 740.0 784.0];
E4=[329.6 349.2 370.0 392.0 415.3 440.0 466.2 493.9 523.3 554.4 587.3 622.3 659.3 698.5 740.0 784.0 830.6 880.0 932.3 987.8 1047.0];

% input string fundamental frequency
prompt1 = 'Please enter string number (1=E2 6=E4):';
string = input(prompt1);

prompt2 = 'Please enter fret number (1=open string 21=last fret):';
fret = input(prompt2);

switch string
    
    case 1
        f0=E2(fret);
    case 2
        f0=A2(fret);   
    case 3
        f0=D3(fret);
    case 4
        f0=G3(fret);
    case 5
        f0=B3(fret);    
    case 6
        f0=E4(fret);
end

fprintf('Working on fundamental frequency: %f Hz\n',f0);

frameSizeMS = 30; % minimum frame length, in ms
overlap = 0.83; % fraction of frame overlapping
windowType = 'blackman'; % type of windowing used for each frame

[signal,fs,bits] = wavread(FileName);
L  = round(fs/f0);
signal = signal(1:floor(2*fs),1);

% calculate STFT frames
minFrameLen = fs*frameSizeMS/1000; 
frameLenPow = nextpow2(minFrameLen);
M = 2^frameLenPow; % frame length = fft size
eval(['frameWindow = ' windowType '(M);']);

R=floor(overlap*M);

[B,F,T] = spectrogram(signal,frameWindow,R,2*M,fs);

% Magnitude in Decibels

B_dB = 20*log10(abs(B));
offset = max(max(B_dB));
B_dBN = B_dB-offset;

% Maximum number of peaks for frame

max_peaks=200;
nframes=size(B,2);

peaks_matrix=zeros(2*nframes,max_peaks); 

% peaks locations matrix (odd rows: peak value , even rows: peak 
% location - column peaks, rows frames)

% the following loop computes for each frame the location of freq. spectrum
% peaks filling peaks_matrix

for frame_i=1:nframes
    [pks,locs] = findpeaks(B_dBN(:,frame_i),'minpeakheight',-50);
    
    npeaks = length(pks);
    
    for i=1:npeaks
        peaks_matrix(2*frame_i-1,i)=pks(i);
        peaks_matrix(2*frame_i,i)=locs(i);
    end        
end

locations=[];

% this section computes how many frequency index occurences are present in
% peaks_matrix ie finding unique entries

for index=1:nframes
    locations=[locations peaks_matrix(2*index,:)];
end

unique_v = unique(locations);

if unique_v(1) == 0
    unique_v = unique_v(2:length(unique_v));
end

unique_count = length(unique_v);
loc_count = [unique_v ; histc(locations,unique_v)];
fprintf('computed %f unique partials\n',size(loc_count,2));

% good_matrix reorder all data per partial - odd row partial magnitude, 
% even row time index

good_matrix=zeros(2*unique_count,max(loc_count(2,(2:unique_count))));
row_pointer=1;

t_offset=M/(2*fs);
t_hop=R/fs;

for i=1:unique_count        
    pointer=1;
    disp(sprintf('working on %.1f frequency index',F(unique_v(i))));
    
    for j=1:nframes
        for k=1:max_peaks
            if peaks_matrix(2*j,k)==unique_v(i)
                
                good_matrix(2*row_pointer-1,pointer)=peaks_matrix(2*j-1,k);
                good_matrix(2*row_pointer,pointer)=(t_offset+(j-1)*t_hop);    
                pointer=pointer+1;
                
                break
            end
        end
    end
    row_pointer=row_pointer+1;
end

fprintf('computed peaks and frequency matrix\n',row_pointer-1);

% fitting data with polyfit and calculating magnitude
fprintf('fitting data with polyfit\n',row_pointer-1);
row_pointer=1;

for i=1:10
    if loc_count(2,i) > 10
    figure(1)
    plot(good_matrix(2*i,1:loc_count(2,i)),good_matrix(2*i-1,1:loc_count(2,i)));
    hold on
    end
end
hold off

for i=1:unique_count
    
    if loc_count(2,i) > 10
        [m,b]=polyfit(good_matrix(2*i,1:loc_count(2,i)),good_matrix(2*i-1,1:loc_count(2,i)),1);
        if m(1) < 0
            G(1,row_pointer)=10^((m(1)*L)/(20*R));
            G(2,row_pointer)=F(loc_count(1,i));
            row_pointer=row_pointer+1;
        end
    end
    
end

fprintf('computed Gk for %f partials\n',row_pointer-1);
