% -------------------------------------------------------
%
%    ECG_Baseline_Removal  - Filter baseline wander using median filters with overlaping windows. 
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
% 
% [filtered_signal,baseline]=ECG_Baseline_Removal(signal,samplerate,window_length,overlap)
% Filter baseline wander using median filters with overlaping windows. 
%
% Function used to filter baseline wander using median filters with overlaping windows. 
% A window is used to calculate the local median. This value is placed in
% the center of the window. The next window is chosen in such a manner that
% a partial overlap exist between the two windows. Afterwards, the next median is calcualted.
% In the end, the obtained values are interpolated to create a baseline estimation with
% the same length as the given signal.	
%
% Inputs:      
%       signal: Single or multichannel ECG signal. The channel must be a column vector.
%       samplerate: Sample frequency of the ECG signal
%       window_length: Window length to obtain the local medians. The length must
%           be given in seconds.
%       overlap: Percentual overlap between neighboring windows. The overlap must
%           be number between 0 and 1. If set to 1 the windows will be centered at
%           every sample value of the given signal.
%
% Outputs:
%       filtered_signal: The filtered signal after baseline removal
%       baseline: Baseline estimation calculated using the mentioned algorithm.
%
% Example Usage:
%       [filtered_signal,baseline]=ECG_Baseline_Removal(signal,250,1,0.5)
%
% Revision history:
%
%
function [filtered_signal,baseline]=ECG_Baseline_Removal(signal,samplerate,window_length,overlap)

%get signal properties
[L,NCH]=size(signal); %number of channels and length of signal
baseline=zeros(size(signal)); %alocate matrix for baseline
filtered_signal=zeros(size(signal)); %alocate filtered signal

window_length=round(window_length*samplerate); %windlow length in sample values
window_length=window_length+1-mod(window_length,2); %window length is set to next odd number
window_half_length=(window_length-1)/2; %half window length

if overlap>=0 && overlap<1
    N=floor((L-window_length*overlap)/(window_length*(1-overlap))); %number of windows to fill signal
    center=round(window_length*(1-overlap)*(0:1:N-1))'+window_half_length+1; %center for window placement
    %center=cat(1,1,center);%First and final values are also center of local smaller windows
    %center=cat(1,center,L);
    %center=[1;center;L]; 
elseif overlap==1
    center=(1:L)'; %number of windows to fill signal
    %center=cat(1,1,center);%First and final values are also center of local smaller windows
    %center=cat(1,center,L);
    %center=[1;center;L]; 
    N=length(center); %number of windows to fill signal
else
    error('overlap must be a number between 0 and 1')
end
    
for j=1:NCH
    baseline_points=zeros(size(center));%alocate memory for baseline points
    %signal_extended=wextend(1,'sp0',signal(:,j),2*window_half_length, 'b');
    for i=1:N
        leftInt=max(center(i)-window_half_length,1);
        rightInt=min(center(i)+window_half_length,length(signal));
%         if center(i)-window_half_length<0
%             baseline_points(i)=median(signal(1:ceil(window_length/4),j));%local median beginng of signal
%         elseif i==N+2
%             baseline_points(i)=median(signal(end-ceil(window_length/4):end,j));%local median end of signal
%         else
            baseline_points(i)=median(signal(leftInt:rightInt,j));%local median
%         end
%         baseline_points(i)=median(signal_extended(center(i):center(i)+2*window_half_length,1));
    end

    baseline(:,j)=pchip(center,baseline_points,1:1:L)'; %interpolate baseline
    filtered_signal(:,j)=signal(:,j)-baseline(:,j); %subtract estimated signal from baseline
    
    %Correction of constant offset
    [filtered_signal(:,j),offset]=Isoline_Correction(filtered_signal(:,j));
    baseline(:,j)=baseline(:,j)+offset;
    
end

