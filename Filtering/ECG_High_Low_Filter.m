% -------------------------------------------------------
%
%    ECG_High_Low_Filter  - This function uses a Fourier-based algorithm to linearly filter a signal
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
% [filteredsignal]=ECG_High_Low_Filter(signal,samplerate,highpass_frequency,lowpass_frequency,varargin)
% This function uses a Fourier-based algorithm to linearly filter a signal 
%
% This function uses a Fourier-based algorithm to linearly filter a signal. 
% It implements a low-pass filter, a high-passfilter and a bandstop filter. 
% The cut-off frequencies can be freely choosen for better performance.
% The user can also choose among three different filters for the low pass 
% filter. Use the third input of the function to define the type of filter. 
% Possible filters are: 'Smooth', 'Gauss', and 'Butterworth'. The default 
% type is Butterworth.
%
% Inputs:
%       signal: ecg signal to be filtered. In the multichannel case every column
%           of the signal matrix should be a channel.
%       samplerate: sample frequency in Hz of the signal 
%       highpass_frequency: cutoff frequency of the highpass filter 
%       lowpass_frequency: cutoff frequency of the lowpass filter 
%       varargin: type of filter 'butterworth', 'gauss', 'smooth'
%
%
% Outputs:
%       filtered_signal: filtered ecg signal.
%
%
% Example Usage:
%       [filteredsignal]=ECG_High_Low_Filter(signal,250,0.5,150)
%
% Revision history:
%
%

function [filteredsignal]=ECG_High_Low_Filter(signal,samplerate,highpass_frequency,lowpass_frequency,varargin)

%In case no filter type is given. Butterworth is used.
if isempty(varargin)
    varargin{1}='Butterworth';
end

%Highpass filtering
filteredsignal=ECG_High_Filter(signal,samplerate,highpass_frequency,varargin{1});

%Lowpass filtering
filteredsignal=ECG_Low_Filter(filteredsignal,samplerate,lowpass_frequency,varargin{1});
