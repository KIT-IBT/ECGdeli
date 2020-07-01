% -------------------------------------------------------
%
%    Check_T_Wave - Recheck the exact postion of the T waves
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
% FPT=Check_T_Wave(signal,samplerate,FPT,varargin)
% Recheck the exact postion of the T waves
%
% This function is used to recheck the exact postion of the T waves. It is
% especially useful after a synchronization of channels or in the case of
% labeled ecgs, where the T waves were not exactly marked. It should also be
% used after an ecg has been filtered and the FPT should be updated. This
% function is very similar to the function T_Detection.
% 
% Inputs:
%       signal: filtered ECG. The baseline wander
%           and high frequency noise must be removed first from the raw data
%   samplerate: sample frequency of the ECG signal
%   varargin: time in seconds of the vecinity of the
%       current location of the T peak, where the exact position should be
%       contained in. Default value is 0.1s
%
% Outputs:
%       FPT: Fiducial Point Table corresponding to the ECG signal
%
%
% Example Usage:
%       FPT_c=Check_T_Wave(signal,250,FPT)
%
% Revision history:
%
%

function FPT=Check_T_Wave(signal,samplerate,FPT,varargin)

display('Checking T Waves...')

%check if signal is of type double
if ~isa(signal,'double')
    signal=double(signal);
end
if ~isvector(signal)
    error('The input ECG signal must be a vector!');
end

%Check if minimum amount of R peaks is present
if size(FPT,1)<3
    disp('Warning: Too little QRS complexes were detected. Returning an empty FPT table')
    FPT=[];
    return
end

%% Filter signal and remove P wave and QRS complexes
%Low and highpass filter together with baseline removal
signal=ECG_High_Low_Filter(signal,samplerate,0.3,20,'B');
w=750e-3;
[~,baseline]=ECG_Baseline_Removal(signal,samplerate,w,0.75);
w=2;
[~,baseline]=ECG_Baseline_Removal(baseline,samplerate,w,0.75);
signal=signal-baseline;

%QRS complexes are removed from the signal. After removal only T waves are
%present. They are then easier to detect.
signal=Remove_PQRS(signal,samplerate,FPT);

%% Wavelet tranformation
% x: level in which we find frequency contents around 7 Hz. Here we find
% the most important components for QRS detection.
x=round(log2(samplerate/2/7));
if x<1
    x=1;
end

%% ECG signal is extended to have a length equal to the next power of 2
l=2^ceil(log2(length(signal)));
l1=floor((l-length(signal))/2);
l2=l-length(signal)-l1;
ecg_w = wextend(1,'sp0',signal,l1,'r');
ecg_w = wextend(1,'sp0',ecg_w,l2,'l');

%% Discrete Wavelet Transformation
% using db1 = Haar wavelet from level 1 up to to level=x
[~,swd]=swt(ecg_w,x,'rbio3.3');
Dx1=swd(x,:)';
Dx1=Dx1((l2+1:end-l1));

t=length(ecg_w):-1:1;
[~,swd]=swt(ecg_w(t),x,'rbio3.3');
Dx2=swd(x,:)';
Dx2=Dx2(t);
Dx2=Dx2((l2+1:end-l1));
Dx=-(Dx1+Dx2);

[~,locs]=findpeaks(abs(Dx));
TPOS_vector=FPT(:,11);


%Dominant type of normal T wave. A template is built to detect if dominant
%T wave is positive or negative
L_width=round(mean(abs(TPOS_vector(1:end-1)-FPT(1:end-1,6)))+25e-3*samplerate);
R_width=round(0.5*mean(abs(FPT(2:end-1,6)-TPOS_vector(1:end-2))));

TPEAK_vector=zeros(size(TPOS_vector)); %position T peak
TON_vector=TPEAK_vector; %position T onset 
TOFF_vector=TPEAK_vector; %position T offset

%Dominant type of normal T wave peaks. T wave peak is positive or negative
RR=diff(FPT(:,6));
X=[RR(1:end-1),RR(2:end)];
index=1:1:size(X,1);
SCORE=(X-[mean(X(index,1))*ones(size(X,1),1),mean(X(index,2))*ones(size(X,1),1)])*1/sqrt(2)*[1,-1;1,1];
D1=abs(SCORE(:,1));
Thl1=2.5*std(D1);
index=(SCORE(:,1)>=-Thl1 & SCORE(:,2)<=0);
Ind_QRS_normal=find(index)+1;
Ind_QRS_normal=Ind_QRS_normal(2:end-1);
MP=zeros(length(Ind_QRS_normal),2);
for i=1:length(Ind_QRS_normal)
    if TPOS_vector(i)-L_width<1 || TPOS_vector(i)+R_width>length(signal)
        continue
    else
        MP(i,:)=[max(Dx(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width)),min(Dx(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width))];
    end
end
Th11=quantile(MP(:,1),0.25);
Th12=quantile(MP(:,1),0.75);
Th21=quantile(MP(:,2),0.25);
Th22=quantile(MP(:,2),0.75);
MP_plus=MP(MP(:,1)>=Th11 & MP(:,1)<=Th12,1);
MP_minus=MP(MP(:,2)>=Th21 & MP(:,2)<=Th22,2);
if isempty(MP_plus) || isempty(MP_minus)
    T_type=sign(median(MP(:,1))+median(MP(:,2)));
else
    T_type=sign(median(MP_plus)+median(MP_minus));
end

%Find real position of peak of T wave peak 
for i=1:size(FPT,1)-1
    
    L_width=round(TPOS_vector(i)-FPT(i,6)-25e-3*samplerate);    
    R_width=round(0.5*(FPT(i+1,6)-TPOS_vector(i)));
        
    if TPOS_vector(i)-L_width<1 || TPOS_vector(i)+R_width>length(signal)
        TPEAK_vector(i)=TPOS_vector(i);
        TON_vector(i)=FPT(i,6)+2;
        TOFF_vector(i)=TPOS_vector(i)+1;
    else 
        Dy_wave=Dx(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width,1);
        [~,I_ext]=findpeaks(abs(Dy_wave));

        if isempty(I_ext)
            TPEAK_vector(i)=TPOS_vector(i);
        else
            [amplitude,index]=sort(Dy_wave(I_ext));%sort amplitude of posible T peaks

            if T_type>0 
                %closest significant maximum to greatest wavelet coeff.
                pos2=I_ext(index)>=L_width-round(140e-3*samplerate) & I_ext(index)<=L_width+round(140e-3*samplerate);
                I_ext2=I_ext(index(pos2));
                if isempty(I_ext2) && ~all(Dy_wave(I_ext2)<0)
                    pos_correction=I_ext(index(amplitude>=Dy_wave(L_width)))-L_width;
                    [~,a]=min(abs(pos_correction));
                    pos_correction=pos_correction(a);
                    TPEAK_vector(i)=pos_correction+TPOS_vector(i);
                elseif isempty(I_ext2) && all(Dy_wave(I_ext2)<0)
                    pos_correction=I_ext(index)-L_width;
                    [~,a]=min(abs(pos_correction));
                    pos_correction=pos_correction(a);
                    TPEAK_vector(i)=pos_correction+TPOS_vector(i);
                else
                    [~,a]=max(Dy_wave(I_ext2));
                    TPEAK_vector(i)=I_ext2(a)-L_width+TPOS_vector(i);
                end              
            else
                %closest significant maximum to greatest wavelet coeff.
                pos2=I_ext(index)>=L_width-round(140e-3*samplerate) & I_ext(index)<=L_width+round(140e-3*samplerate);
                I_ext2=I_ext(index(pos2));
                if isempty(I_ext2) && ~all(Dy_wave(I_ext2)>0)
                    pos_correction=I_ext(index(amplitude<=Dy_wave(L_width)))-L_width;
                    [~,a]=min(abs(pos_correction));
                    pos_correction=pos_correction(a);
                    TPEAK_vector(i)=pos_correction+TPOS_vector(i);
                elseif isempty(I_ext2) && all(Dy_wave(I_ext2)>0)
                    pos_correction=I_ext(index)-L_width;
                    [~,a]=min(abs(pos_correction));
                    pos_correction=pos_correction(a);
                    TPEAK_vector(i)=pos_correction+TPOS_vector(i);
                else
                    [~,a]=min(Dy_wave(I_ext2));
                    TPEAK_vector(i)=I_ext2(a)-L_width+TPOS_vector(i);
                end
            end
        end
    end
end

%Recalcualte width of T wave using the obtained positions of T wave peaks 
 
threshold_distribution_on=0.3;
threshold_distribution_off=0.3;

%% Delination of T wave
%In every detected T wave the position of the peak is chosen by 
%finding the most significant peak inside the T wave and taking into
%consideration the type of the wave
for i=1:size(FPT,1)-1
    L_width=max(1,round(TPEAK_vector(i)-FPT(i,6)-25e-3*samplerate));    
    R_width=max(1,round(0.5*(FPT(i+1,6)-TPEAK_vector(i))));
    
    Dy_density=abs(Dx(TPEAK_vector(i)-L_width:TPEAK_vector(i)+R_width));
    %Dy_density=Dy_density/;
    Dy_distribution_L=cumsum(Dy_density(1:L_width)/sum(Dy_density(1:L_width)));
    Dy_distribution_R=cumsum(Dy_density(L_width+1:end)/sum(Dy_density(L_width+1:end)));
    
    TON_vector(i)=TPEAK_vector(i)-L_width-1+find(Dy_distribution_L>=threshold_distribution_on,1,'first');
    [~,pos]=min(abs(TON_vector(i)-locs));
    TON_vector(i)=locs(pos);

    if TON_vector(i)<=FPT(i,6)+round(75e-3*samplerate)
        if TPEAK_vector(i)+R_width>length(signal) || TPEAK_vector(i)-L_width<1
            TON_vector(i)=FPT(i,6)+round(75e-3*samplerate)+1;
        else
            TON_vector(i)=TPEAK_vector(i)-L_width-1+find(Dy_distribution_L>=threshold_distribution_on,1,'first');
        end
    end
       
    TOFF_vector(i)=TPEAK_vector(i)+find(Dy_distribution_R>=1-threshold_distribution_off,1,'first');
    [~,pos]=min(abs(TOFF_vector(i)-locs));
    TOFF_vector(i)=locs(pos);

    if TOFF_vector(i)>=FPT(i+1,6)-round(75e-3*samplerate)
        if TPEAK_vector(i)+R_width>length(signal) || TPEAK_vector(i)-L_width<1
            TOFF_vector(i)=FPT(i+1,6)-round(75e-3*samplerate)-1;
        else
            TOFF_vector(i)=TPEAK_vector(i)+find(Dy_distribution_R>=threshold_distribution_off,1,'first');               
        end
    end
end

%Place local T wave delineation in global FPT
FPT(1:end-1,10)=TON_vector(1:end-1);
FPT(1:end-1,11)=TPEAK_vector(1:end-1);
FPT(1:end-1,12)=TOFF_vector(1:end-1);
FPT(end,10)=length(signal)-2;
FPT(end,11)=length(signal)-1;
FPT(end,12)=length(signal);

%Place the L point. The L point is used to diagnose an ST elevation or 
%depresion. It is located around the middle between QRS offset and T onset
if all(FPT(:,10)~=0) && all(FPT(:,8)~=0) %QRS offset and T onset are given
    k=0.5;
    FPT(:,9)=round(k*FPT(:,10)+(1-k)*FPT(:,8)); %L point is located around the middle between QRS offset and T onset 
elseif all(FPT(:,10)~=0) && all(FPT(:,7)~=0) %S wave and T onset are given
    k=0.55;
    FPT(:,9)=round(k*(FPT(:,10)+(1-k)*FPT(:,7))); %L point is located around the middle between QRS offset and T onset 
else %R peak and T onset are given
    k=0.6;
    FPT(:,9)=round(k*FPT(:,10)+(1-k)*FPT(:,6)); %L point is located around the middle between QRS offset and T onset
end

display('Done');
