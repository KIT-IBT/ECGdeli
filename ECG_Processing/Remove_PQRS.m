% -------------------------------------------------------
%
%    Remove_PQRS  - Remove P wave and QRS complex from the ecg signal
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
% [replaced_signal,intervalSigmoid]=Remove_PQRS(signal,samplerate,FPT)
% Remove P wave and QRS complex from the ecg signal
%
% Remove P wave and QRS complex from the ecg signal
%
% Inputs:
%       signal: ECG signal
%       samplerate: sampling rate of the ECG 
%       FPT: Fiducial Point Table
%
% Outputs:
%       replaced_signal: ECG with replaced P wave and QRS complex
%       intervalSigmoid: intervals that are replaced by the sigmoid
%       function
%
%
% Example Usage:
%       [replaced_signal,intervalSigmoid]=Remove_PQRS(signal,250,FPT)
%
% Revision history:
%
%

function [replaced_signal,intervalSigmoid]=Remove_PQRS(signal,samplerate,FPT)

%check if signal is of type double
if ~isa(signal,'double')
    signal=double(signal);
end
if ~isvector(signal)
    error('The input ECG signal must be a vector!');
end

%Baseline removal
w=750e-3;
[~,baseline]=ECG_Baseline_Removal(signal,samplerate,w,0.75);
w=2;
[~,baseline]=ECG_Baseline_Removal(baseline,samplerate,w,0.75);
signal=signal-baseline;

%Alocate variables
replaced_signal=signal;
intervalSigmoid=nan(size(FPT,1),2);

%Paramter estimation for every QRS complex and removal using a baseline
%compensation first and an approximation of the QRS complex using the 
%Hermite basis functions
for i=1:size(FPT,1)
    
    %Estimate the with of a normal QRS complex
    if i==1
        Remove_L=round(1/3*(FPT(i,6)));
    else
        Remove_L=round(1/3*(FPT(i,6)-FPT(i-1,6)));
    end
    if ~FPT(i,7)
        Remove_R=round(80e-3*samplerate);
    else
        Remove_R=round(FPT(i,7)-FPT(i,6)+50e-3*samplerate);
    end

    %Parameter for baseline compensation
    n1=ceil(samplerate*10e-3);
    n2=ceil(samplerate*10e-3);
    xc=linspace(-6,6,Remove_L+Remove_R+1)';
    c=1./(1+exp(-xc));
    
    if FPT(i,6)-Remove_L<1
        PQRS=replaced_signal(1:FPT(i,6)+Remove_R);
        y1=mean(PQRS(1:n1));
        y2=mean(PQRS(end-n2+1:end));
        replaced_signal(1:FPT(i,6)+Remove_R)=(y2-y1)*c(end-length(PQRS)+1:end)+y1;
        intervalSigmoid(i,1)=1;
        intervalSigmoid(i,2)=FPT(i,6)+Remove_R;
    elseif FPT(i,6)+Remove_R>length(signal)
        PQRS=replaced_signal(FPT(i,6)-Remove_L:end);
        y1=mean(PQRS(1:n1));
        y2=mean(PQRS(end-n2+1:end));
        replaced_signal(FPT(i,6)-Remove_L:end)=(y2-y1)*c(1:length(PQRS))+y1;
        intervalSigmoid(i,1)=FPT(i,6)-Remove_L;
        intervalSigmoid(i,2)=size(replaced_signal,1);
    else
        PQRS=replaced_signal(FPT(i,6)-Remove_L:FPT(i,6)+Remove_R);
        y1=mean(PQRS(1:n1));
        y2=mean(PQRS(end-n2+1:end));
        replaced_signal(FPT(i,6)-Remove_L:FPT(i,6)+Remove_R)=(y2-y1)*c+y1;
        intervalSigmoid(i,1)=FPT(i,6)-Remove_L;
        intervalSigmoid(i,2)=FPT(i,6)+Remove_R;
    end
end

