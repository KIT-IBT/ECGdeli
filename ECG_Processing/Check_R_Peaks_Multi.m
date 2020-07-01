% -------------------------------------------------------
%
%    Check_R_Peaks_Multi  - Rechecking the exact postion of the R peaks
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
% [FPT]=Check_R_Peaks_Multi(signal,samplerate,FPT,varargin)
% Rechecking the exact postion of the R peaks 
%
% This function is used to recheck the exact postion of the R peaks. It is
% especially useful after a synchronization of channels or in the case of
% labeled ecg, where the R-peaks were not exactly marked. It should also be
% used when an ecg has been filtered and the FPT should be updated. This
% function is also the final part of the function QRS_Detection.	
%
% Inputs:
%       signal: single channel filtered ECG as a standing 
%           vector. The baseline wander and high frequency noise must be removed first 
%           from the raw data.
%       samplerate: sampling frequency (in Hz) of the ECG signal.
%       varargin: time in seconds of the vecinity of the
%           current location of the R peak, where the exact position should be
%           contained in. Default value is 0.09s.
%       FPT: Fiducial Point Table corresponding to the ECG signal.
%
% Outputs:
%       FPT: corrected Fiducial Point Table corresponding to the ECG signal.
%
%
% Example Usage:
%       [FPT_c]=Check_R_Peaks_Multi(signal,250,FPT)
%
% Revision history:
%
%

function [FPT]=Check_R_Peaks_Multi(signal,samplerate,FPT,varargin)
%This function is used to recheck the exact postion of the R peaks. It is
%specially useful after a synchronization of channels or in the case of
%labeled ecg, where the R-peaks were not exactly marked. It should also be
%used when an ecg has been filtered and the FPT should be updated. This
%function is also the final part of the function QRS_Detection.
%Inputs:
%signal: This variable is a single channel filtered ECG as a standing 
%vector The baseline wander and high frequency noise must be removed first 
%from the raw data.
%samplerate: This variable is the sample frequency of the ECG signal.
%varargin: This varible is the time in seconds of the vecinity of the
%current location of the R peak, where the exact position should be
%contained in. Default value is 0.09s.
%Outputs:
%FPT is the Fidutial Point Table corresponding to the ECG signal.

disp('Checking R Peaks...');

%% Denoise ECG
%Highpass 5Hz to reduce baseline wander, aminorate the T and P waves and
%create a QRS complex with strong Q and S waves
highpass_frequency=0.5;
%Filtering
if samplerate>500
    lowpass_frequency=250;
    %High and low pass filter to reduce noise
    signal=ECG_High_Low_Filter(signal,samplerate,highpass_frequency,lowpass_frequency,'B');
else
    %High pass filter to reduce noise
    signal=ECG_High_Filter(signal,samplerate,highpass_frequency,'B');
end
%% Delineate QRS complexes
%Vector for the positions of the R-peaks
RPOS_vector=zeros(size(FPT,1),1);
QPOS_vector=RPOS_vector;
SPOS_vector=RPOS_vector;
R_Synced=FPT(:,6);

%Check if minimum amount of R peaks is present
if length(R_Synced)<3 
    disp('Warning: Too little QRS complexes were detected. Returning an empty FPT table')
    FPT=[];
    return
end

%Size of the window in order to look for the exact position of the R-peak
%in the proximity of the already found point inside the QRS complex
if isempty(varargin)
    WB=ceil(50e-3*samplerate);
else
    WB=ceil(varargin{1}*samplerate/2);
end
QRS_region(:,1)=FPT(:,6)-WB;
QRS_region(:,2)=FPT(:,6)+WB;

%calculate zero crossings of derivative of ecg for QRS Bondaries
dsignal=diff(signal);
i=1:length(dsignal)-1;
I_ext=find((dsignal(i)>=0 & dsignal(i+1)<0) | (dsignal(i)<0 & dsignal(i+1)>=0))+1;

%Dominant type of normal R peaks. A template is built to detect if dominant
%R peak is positive or negative
RR=diff(R_Synced);
X=[RR(1:end-1),RR(2:end)];
index=1:1:size(X,1);
SCORE=(X-[mean(X(index,1))*ones(size(X,1),1),mean(X(index,2))*ones(size(X,1),1)])*1/sqrt(2)*[1,-1;1,1];
D1=abs(SCORE(:,1));
Thl1=2.5*std(D1);
index=(SCORE(:,1)>=-Thl1 & SCORE(:,2)<=0);
Ind_QRS_normal=find(index)+1;
Ind_QRS_normal=Ind_QRS_normal(2:end-1);
QRS_Matrix=zeros(2*WB+1,length(Ind_QRS_normal));
MP=zeros(length(Ind_QRS_normal),2);
for k=1:length(Ind_QRS_normal)
    QRS_Matrix(:,k)=signal(R_Synced(Ind_QRS_normal(k))-WB:R_Synced(Ind_QRS_normal(k))+WB);
    MP(k,:)=[max(QRS_Matrix(:,k)),min(QRS_Matrix(:,k))];
end
Th11=quantile(MP(:,1),0.25);
Th12=quantile(MP(:,1),0.75);
Th21=quantile(MP(:,2),0.25);
Th22=quantile(MP(:,2),0.75);
QRS_Matrix_selected=QRS_Matrix(:,MP(:,1)>=Th11 & MP(:,1)<=Th12 & MP(:,2)>=Th21 & MP(:,2)<=Th22);
if isempty(QRS_Matrix_selected)
    Template=mean(QRS_Matrix,2);
else
    Template=mean(QRS_Matrix_selected,2);
end
R_type=sign(max(Template)+min(Template));

%In every detected QRS complex the position of the R-peak is chosen by 
%finding the most significant peak inside the QRS complex and taking into
%consideration the type of the complex
biph_crit=2/5; %Criterium for biphasic QRS complexes
w_crit=9/10; %Criterium for w-shape QRS complexes
for i=1:length(RPOS_vector)
    
    %find extrema in interval number i
    tmp_ZC=find(I_ext>=QRS_region(i,1)-WB & I_ext<=QRS_region(i,2)+WB);
    
    if isempty(tmp_ZC)
        RPOS_vector(i)=round((QRS_region(i,1)+QRS_region(i,2))/2);
        QPOS_vector(i)=QRS_region(i,1); %position Q peak in case no Q peak was found
        SPOS_vector(i)=QRS_region(i,2); %position S peak in case no S peak was found
    elseif length(tmp_ZC)==1 %in case only 1 extremum is found in the given interval
        RPOS_vector(i)=I_ext(tmp_ZC);
        QPOS_vector(i)=QRS_region(i,1); %position Q peak in case no Q peak was found
        SPOS_vector(i)=QRS_region(i,2); %position S peak in case no S peak was found
    else
        [amplitude,index]=sort(signal(I_ext(tmp_ZC)));%sort amplitude of posible R peaks
        
        if min(abs(amplitude(1)/amplitude(end)),abs(amplitude(end)/amplitude(1)))>biph_crit %biphasic QRS complex
            if R_type>=0
                if abs(amplitude(end-1)/amplitude(end))<w_crit
                    RPOS_vector(i)=I_ext(tmp_ZC(index(end)));%position of the maximum (R peak)
                    Qpeak=index(end)-1;
                    Speak=index(end)+1;
                else %"M" like QRS complex
                    RPOS_vector(i)=min(I_ext(tmp_ZC(index(end))),I_ext(tmp_ZC(index(end-1))));%position of the maximum (R peak)
                    Qpeak=min(index(end-1),index(end))-1;
                    Speak=max(index(end-1),index(end))+1;
                end      
            else
                if abs(amplitude(2)/amplitude(1))<w_crit 
                    RPOS_vector(i)=I_ext(tmp_ZC(index(1)));%position of the maximum (R peak)
                    Qpeak=index(1)-1;
                    Speak=index(1)+1;
                else %"W" like QRS complex
                    RPOS_vector(i)=min(I_ext(tmp_ZC(index(1))),I_ext(tmp_ZC(index(2))));%position of the minimum (R peak)
                    Qpeak=min(index(2),index(1))-1;
                    Speak=max(index(2),index(1))+1;
                end
            end
            if Qpeak>0
                QPOS_vector(i)=I_ext(tmp_ZC(Qpeak)); %position of Q peak
            else
                QPOS_vector(i)=RPOS_vector(i)-WB; %position Q peak in case no Q was found
            end
            if Speak<=length(tmp_ZC)
                SPOS_vector(i)=I_ext(tmp_ZC(Speak)); %position of S peak
            else
                SPOS_vector(i)=RPOS_vector(i)+WB; %position S peak in case no S peak was found
            end
        elseif abs(amplitude(end))>abs(amplitude(1)) %positiv QRS complex
            if abs(amplitude(end-1)/amplitude(end))<w_crit
                RPOS_vector(i)=I_ext(tmp_ZC(index(end)));%position of the maximum (R peak)
                Qpeak=index(end)-1;
                Speak=index(end)+1;
            else %"M" like QRS complex
                RPOS_vector(i)=min(I_ext(tmp_ZC(index(end))),I_ext(tmp_ZC(index(end-1))));%position of the maximum (R peak)
                Qpeak=min(index(end),index(end-1))-1;
                Speak=max(index(end),index(end-1))+1;
            end
            if Qpeak>0
                QPOS_vector(i)=I_ext(tmp_ZC(Qpeak)); %position of Q peak
            else
                QPOS_vector(i)=RPOS_vector(i)-WB; %position Q peak in case no Q was found
            end
            if Speak<=length(tmp_ZC)
                SPOS_vector(i)=I_ext(tmp_ZC(Speak)); %position of S peak
            else
                SPOS_vector(i)=RPOS_vector(i)+WB; %position S peak in case no S peak was found
            end
        else %negativ QRS complex
            if abs(amplitude(2)/amplitude(1))<w_crit
                RPOS_vector(i)=I_ext(tmp_ZC(index(1)));%position of the maximum (R peak)
                Qpeak=index(1)-1;
                Speak=index(1)+1;
            else %"W" like QRS complex
                RPOS_vector(i)=min(I_ext(tmp_ZC(index(1))),I_ext(tmp_ZC(index(2))));%position of the maximum (R peak)
                Qpeak=min(index(2),index(1))-1;
                Speak=max(index(2),index(1))+1;
            end
            if Qpeak>0
                QPOS_vector(i)=I_ext(tmp_ZC(Qpeak)); %position of Q peak
            else
                QPOS_vector(i)=RPOS_vector(i)-WB; %position Q peak in case no Q was found
            end
            if Speak<=length(tmp_ZC)
                SPOS_vector(i)=I_ext(tmp_ZC(Speak)); %position of S peak
            else
                SPOS_vector(i)=RPOS_vector(i)+WB; %position S peak in case no S peak was found
            end
        end
    end
end

%Place QRSon 20ms prior to Q wave and QRSoff 20 ms after S wave
donoff=round(20e-3*samplerate);
QRSonPOS_vector=QPOS_vector-donoff;
QRSoffPOS_vector=SPOS_vector+donoff;

%Recheck QRSon and QRSoff of first and last peak
if ~isempty(R_Synced)
    if QRSonPOS_vector(1)<1
        QRSonPOS_vector(1)=1;
    end
    if QRSoffPOS_vector(end)>length(signal)
        QRSoffPOS_vector(end)=length(signal);
    end

    %Recheck Q and S waves of first and last peak
    if QPOS_vector(1)<2 
        QPOS_vector(1)=2;
    end
    if SPOS_vector(end)>length(signal)-1
        SPOS_vector(end)=length(signal)-1;
    end
end

%Place vectors with QRS peaks in the FPT table
FPT(:,4)=QRSonPOS_vector;
FPT(:,5)=QPOS_vector;
FPT(:,6)=RPOS_vector;
FPT(:,7)=SPOS_vector;
FPT(:,8)=QRSoffPOS_vector;

disp('Done');