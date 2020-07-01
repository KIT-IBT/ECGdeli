% -------------------------------------------------------
%
%    Notch_Filter  - filter the artifacts seen in the ECG signal with a constant frequency of f0 Hz.
%
%    Ver. 1.0.0
%
%    Created:           Nicolas Pilia (22.06.2020)
%    Last modified:     Claudia Nagel (29.06.2020)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2020 - All rights reserved.
%
% ------------------------------------------------------
%
% Filtered_Signal=Notch_Filter(signal,samplerate,f0,width)
% Filter the artifacts seen in the ECG signal with a constant frequency of f0 Hz.
%
%
% This function is used to filter the artifacts seen in the ECG signal
% with a constant frequency of f0 Hz. They can be removed using the band 
% stop filter implented here. 
%
% Inputs:
%       signal: Multilead or single lead ecg signal. Every channel must be a
%           standing vector.
%       samplerate: Sample frequency of the signal.
%       f0: Frequency in Hz at which the peaks occur
%       width: bandwith in Hz of the notch filter around the f0
%
% Outputs:
%       Filtered_Signal: filtered ecg signal.
%
%
% Example Usage:
%       Filtered_Signal=Notch_Filter(signal,250,50,1)
%
% Revision history:
%
%

function Filtered_Signal=Notch_Filter(signal,samplerate,f0,width)
%The spectrum will have peaks at k*f0Hz. K gives the greatets number n 
%that can be chosen for a harmonic oscillation without going beyond the 
%nyquist frequency 
K=floor(samplerate/2*1/f0);

%Extend signal to avoid boundary effects
extpoints=round(0.5*ceil(samplerate/width));
signal_extended=zeros(size(signal,1)+2*extpoints,size(signal,2));
for i=1:size(signal,2)
    signal_extended(:,i)=wextend('1D','sp0',signal(:,i),extpoints);
end

L=size(signal_extended,1);%length of the signal
f=(0:1:L-1)/L*samplerate;%Frequcuency vector

sigmaf=width;%Standard deviation of gaussian bell used to select frequency
sigma=ceil(L*sigmaf/samplerate);%Sigma discrete
lg=2*round(4*sigma)+1;%Size of gaussian bell
lb=(lg-1)/2;%Position of center of guassian bell
g=fspecial('gaussian',[1,lg],sigma)';%Gaussian bell
g=1/(max(g)-min(g))*(max(g)-g);%Scale gaussian bell to be in interval [0;1]

H=ones(size(signal_extended,1),1);%Filter

%Implementation of periodical gaussian bells at k*f0Hz
for k=1:K
        [~,b]=min(abs(f-k*f0));%Discrete position at which f=k*f0Hz
        H(b-lb:b+lb)=g;%Gaussian bell placed around k*f0Hz
        H(L+2-b-lb:L+2-b+lb)=g;%Gaussian bell placed symmetriclly around samplerate-k*f0Hz   
   
        %Uncomment the following lines if you want to use rectangles instead of gaussian bells           
        %H(b-lb:b+lb)=0;
        %H(L+2-b-lb:L+2-b+lb)=0;
end

H=repmat(H,1,size(signal_extended,2));%Reproduce the filter for all channels
X=fft(signal_extended);%FFT of signal
Y=H.*X;%Filtering process in the Fourier Domain
Filtered_Signal=real(ifft(Y));%Reconstruction of filtered signal
Filtered_Signal=Filtered_Signal(extpoints+1:end-extpoints,:);
