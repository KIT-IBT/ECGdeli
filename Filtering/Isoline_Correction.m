% -------------------------------------------------------
%
%    Isoline_Correction  - Function used to estimate and remove the offset of the electrical isoline of an ECG.
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
% [filteredsignal,offset,frequency_matrix,bins_matrix]=Isoline_Correction(signal,varargin)
% Function used to estimate and remove the offset of the electrical isoline of an ECG.
%
% Function used to estimate and remove the offset of the electrical isoline 
% of an ECG. The algorithm estimates the isoline by finding the mode (most 
% frequent amplitude) of the amplitude distribution of the ECG signal and removes it. 
%
% Inputs:
%       signal: multi or single channel ECG. Every channel of the signal must be a
%           column of this matrix
%       varagin: optional input used to give the number of bins used to create the
%           histogram for the estimation of the mode. The default value is 2^10=1024
%
% Outputs:
%       filteredsignal: multi or single channel filtered signal for wich the 
%           offset has been compensated in every channel. The isoline of each ECG 
%           should be now around zero
%       offset: vector or scalar containing the estimated offset (most frequent
%           amplitude) of the each channel in the ECG signal
%       frequency_matrix: the columns of this matrix contain the frequencies used 
%           for the histogram of each ECG channel  
%       bin_matrix: the columns of this matrix contain the amplitudes used for
%           the histogram of each ECG channel
%
%
% Example Usage:
%       [filteredsignal,offset,frequency_matrix,bins_matrix]=Isoline_Correction(signal)
%
% Revision history:
%
%


function [filteredsignal,offset,frequency_matrix,bins_matrix]=Isoline_Correction(signal,varargin)

%Alocate filtered signal
filteredsignal=zeros(size(signal));
%Number of channels in ECG
number_channels=size(signal,2);

%Check of optional input
if isempty(varargin)
    number_bins=min(2^10,size(signal,1)); %default number of bins for histogram
else
    number_bins=varargin{1}; %given number of bins for histogram
end

%Alocate matrix for histogram frequencies
frequency_matrix=zeros(number_bins,number_channels);
%Alocate matrix for bin centers
bins_matrix=frequency_matrix;
%
offset=zeros(number_channels,1);

%Constant offset removal
for i=1:number_channels
    [frequency_matrix(:,i),bins_matrix(:,i)]=hist(signal(:,i),number_bins); %Create histogram
    [~,pos]=max(frequency_matrix(:,i)); %find maximum of histogram
    offset(i)=bins_matrix(pos,i); %find most frequent amplitude in the ECG signal
    filteredsignal(:,i)=signal(:,i)-offset(i); %remove offset
end

