% -------------------------------------------------------
%
%    QRS_Detection  - Detection of the QRS complex in an ECG signal
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
% FPT=QRS_Detection(signal,samplerate,varargin)
% Detection of the QRS complex in an ECG signal
%
% Function used to detect and delineate the QRS complexes peaks in an ECG signal.	
%
% Inputs:
%       signal: signal containing an ECG
%       samplerate: sample frequency used to measure the signal
%       definition of the QRS peaks: set varargin{1} to 'peaksQRS' to use
%           the maximum value of the QRS complex as R peak, and the minima for
%           Q and S. Otherwise the dominant type and the shape of the R peaks are
%           evaluated. 
%       muting option: type 'm' or 'mute' for supressing information (optional).
%           Errors/Warnings are shown.
%
% Outputs:
%       FPT: Fidutial Point Table containing the position of the R peaks in the
%           signal. The sixth column of FPT is reserved for this purpose.
%
%
% Example Usage:
%       QRS_Detection(signal,250,varargin)
%
% Revision history:
%
%

function FPT=QRS_Detection(signal,samplerate,varargin)
flag_posR = false; 
if ~isempty(varargin)
    if strcmp(varargin{1},'peaksQRS')
        flag_posR = true; 
    end
    
    if size(varargin,1) == 2 && strcmp(varargin{2},'mute') || strcmp(varargin{1},'m')
        mute=1;
    else
        %warning('Illegal function parameter. Muting is not activated.')
        mute=0;
    end
else
    mute=0;
end

if mute==0
    disp('Detecting R Peaks...')
end

%check if signal is of type double and different than cero
if ~isa(signal,'double')
    signal=double(signal);
end
if ~isvector(signal)
    error('The input ECG signal must be a vector!');
elseif size(signal,2)>size(signal,1)
    signal=signal';
end
if all(abs(signal)<eps)
    warning('The signal values are too small to process. Returning empty FPT table')
    FPT=[];
    disp('Done')
    return
end

%% Denoise ECG
%Highpass 5Hz to reduce baseline wander, aminorate the T and P waves and
%create a QRS complex with strong Q and S waves
highpass_frequency=0.5;
%Lowpass 30Hz to reduce high frequency noise
lowpass_frequency=30;
%Filtering
signal=ECG_High_Low_Filter(signal,samplerate,highpass_frequency,lowpass_frequency,'B');

%% Downsampling
%If samplerate>400Hz the signal is downsampled to 400Hz
%Downsampling accelerates QRS detection
%Use Check_R_Peaks after QRS detection if you need more precision
fdownsample=400;
flagdownsample=false;
if samplerate>fdownsample
    oldsamplerate=samplerate;
    oldsignal=signal;
    
    r=floor(samplerate/fdownsample);
    signal=decimate(signal,r); %downsampling
    samplerate=samplerate/r; %new sampling rate
    flagdownsample=true;
end

%% Wavelet tranformation
% x: level in which we find frequency contents around 75 Hz. Here we find
% the most important components for QRS detection.
x = ceil(log2(samplerate/2/30));
if x<1
    x=1;
end

%% ECG signal is extended to have a length equal to the next power of 2
if log2(length(signal))==nextpow2(length(signal))
    l=2^(nextpow2(length(signal))+1);
else
    l=2^nextpow2(length(signal));
end
l1=floor((l-length(signal))/2);
l2=l-length(signal)-l1;
ecg_w = wextend(1,'sp0',signal,l1,'r');
ecg_w = wextend(1,'sp0',ecg_w,l2,'l');

%% Discrete Wavelet Transformation
% using db1 = Haar wavelet from level 1 up to to level=x
[~,swd]=swt(ecg_w,x,'haar');
Dx=swd(x,:)';
Dx=Dx((l2+1:end-l1));

[~,swd]=swt(flipud(ecg_w),x,'haar');
Dx2=swd(x,:)';
Dx2=flipud(Dx2);
Dx2=Dx2((l2+1:end-l1));
Dx=abs(Dx+Dx2);
Dx=Dx/std(Dx);
saturation=quantile(Dx,0.99);
Dx(Dx>saturation)=saturation;


%% Thresholds and initial values
Thbegin=1;
Thend=quantile(Dx,0.95)/saturation;
threshold=linspace(Thbegin,Thend,20);

Tl=4;
nrep=3;
R_Cell=cell(length(Tl),1);

%% Generate basic rms. It is increased by threshold 0.25 during for loop (H)

for j=1:nrep
    
    NR_vec=zeros(size(threshold));%Vector with number of detected QRS complexes
    
    n1=floor(samplerate*Tl);
    n2=floor(length(signal)/n1)-1;
    for i=1:n2+1
        if ~n2
            rms_Dx_base(:,1)=quantile(Dx(round(0.1*samplerate):length(Dx)-round(0.1*samplerate)),0.95)*ones(size(Dx));
        else
            if i==1
                rms_Dx_base(1:n1,1)=quantile(Dx(round(0.1*samplerate):n1),0.95)*ones(size(Dx((i-1)*n1+1:i*n1)));
            elseif i==n2+1
                rms_Dx_base((i-1)*n1+1:length(Dx),1)=quantile(Dx((i-1)*n1+1:length(Dx)-round(0.1*samplerate)),0.95)*ones(size(Dx((i-1)*n1+1:length(Dx))));
            else
                rms_Dx_base((i-1)*n1+1:i*n1,1)=quantile(Dx((i-1)*n1+1:i*n1),0.95)*ones(size(Dx((i-1)*n1+1:i*n1)));
            end
        end
    end
    
    for H=1:length(threshold)
        %% Generate thresholds for assumed QRS regions
        % the value of "threshold" might be different for diferent leads and
        % patients
        
        %If no optimal threshold value is found, an intermediate value will be
        %chosen. This value should lead to a fairly good detection of the QRS
        %complexes in normal ECGs
        if H==length(threshold)
            [~,mt]=min(diff(NR_vec(1:H-1)));
            rms_Dx=threshold(mt)*rms_Dx_base;
        else
            rms_Dx=threshold(H)*rms_Dx_base;
        end
        
        
        %% Find regions assumed to involve a QRS complex
        %find regions by using all samples higher than the rms_DX threshold
        candidates_Dx=Dx>rms_Dx;
        
        % Signal of the swd around 15Hz with all samples set to zero if not in a
        % QRS region
        Can_Sig_Dx=zeros(size(Dx));
        Can_Sig_Dx(candidates_Dx)=1;
        Can_Sig_Dx(1)=0;
        Can_Sig_Dx(end)=0;
        
        %Get boundaries of regions (Dx)
        i = 1:length(Can_Sig_Dx)-1;
        Bound_A=find(Can_Sig_Dx(i)==0 & Can_Sig_Dx(i+1)>0)+1;
        Bound_B=find(Can_Sig_Dx(i)>0 & Can_Sig_Dx(i+1)==0);
        
        %Combine close regions to one region
        while min(Bound_A(2:end)-Bound_B(1:end-1))/samplerate<0.1
            ind=find((Bound_A(2:end)-Bound_B(1:end-1))/samplerate<0.1);
            Bound_B(ind)=[];
            Bound_A(ind+1)=[];
        end
        
        %Create QRS region and eliminate too long ones
        ind=(Bound_B-Bound_A)/samplerate<5e-3 | (Bound_B-Bound_A)/samplerate>0.25;
        Bound_B(ind)=[];
        Bound_A(ind)=[];
        
        %When the change in the number of detected QRS regions stops
        %decelerating and starts accelerating the adaptive threshold is no
        %longer changed.
        NR_vec(H)=length(Bound_A);
        if H>1
            dNR=NR_vec(H)-NR_vec(H-1);
            if dNR<=0 || H==length(threshold)
                if isempty(Bound_A) || isempty(Bound_B) || size(Bound_A,1) == 1 || size(Bound_B,1) == 1 
                    %For the case that Bound_A or Bound B are empty and no
                    %QRS complex was detected, empty, the next iteration continues
                    continue
                else
                    Tl=quantile(diff(Bound_A')/samplerate,0.98)*4;
                    break
                end
            end
        end
    end
    %Center of region is assumed to be the R peak
    QRS_pos=1/2*(Bound_A+Bound_B);
    R_Cell{j}=round(QRS_pos);
    
end

%% Find Q, R, ans S Peaks

%Sync the QRS detections using different methods
FPT_Cell=cell(length(R_Cell));
for i=1:length(R_Cell)
    FPT_Cell{i}=zeros(length(R_Cell{i}),13);
    FPT_Cell{i}(:,6)=R_Cell{i};
end

FPT_Synced=Sync_R_Peaks(FPT_Cell,samplerate);
R_Synced=FPT_Synced(:,6);
%R_Synced=Sync_R_Peaks(R_Cell,samplerate);

% Remove downsampling
if flagdownsample
    samplerate=oldsamplerate;
    signal=oldsignal;
    R_Synced=R_Synced*r;
end

%Fidutial Point Table is created. If R_Sync does not contain at least 3
%beats it is returned empty
if isempty(R_Synced)
    warning('No QRS complexes were found. Returning an empty FPT table')
    FPT=[];
    return
else
    %Create possible QRS regions to look for extrema
    WB=round(0.05*samplerate);
    QRS_region=[R_Synced-WB,R_Synced+WB];
    
    %Chek first and final R peaks
    if R_Synced(1)-WB<1
        ind=find(R_Synced-WB>=1,1,'first');
        R_Synced=R_Synced(ind:end);
    end
    if R_Synced(1)+WB>length(signal)
        ind=find(R_Synced+WB<=length(signal),1,'last');
        R_Synced=R_Synced(ind:end);
    end
    FPT=zeros(size(R_Synced,1),13);
end
if length(R_Synced)<3
    warning('Too little QRS complexes were detected. Returning an empty FPT table')
    FPT=[];
    return
end

%Vector for the positions of the R-peaks
RPOS_vector=zeros(size(FPT,1),1);
QPOS_vector=RPOS_vector;
SPOS_vector=RPOS_vector;

if flag_posR == false
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
    donoff=round(25e-3*samplerate);
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
else
    QRS_region(QRS_region<=0) = 1; 
    QRS_region(QRS_region>size(signal,1)) = size(signal,1); 
    for i = 1:size(RPOS_vector,1)
        [~, rpeak] = max(signal(QRS_region(i,1):QRS_region(i,2),1));
        RPOS_vector(i,1) = rpeak + QRS_region(i,1) - 1; 
        [~, speak] = min(signal(RPOS_vector(i,1):QRS_region(i,2),1));
        SPOS_vector(i,1) = RPOS_vector(i,1) + speak - 1; 
        [~, qpeak] = min(signal(QRS_region(i,1):RPOS_vector(i,1),1));
        QPOS_vector(i,1) = QRS_region(i,1) + qpeak - 1; 

    end 
    
    donoff=round(25e-3*samplerate);
    QRSonPOS_vector=max(1,QPOS_vector-donoff);
    QRSoffPOS_vector=min(SPOS_vector+donoff,size(signal,1));
end

%Place vectors with QRS peaks in the FPT table
FPT(:,4)=QRSonPOS_vector;
FPT(:,5)=QPOS_vector;
FPT(:,6)=RPOS_vector;
FPT(:,7)=SPOS_vector;
FPT(:,8)=QRSoffPOS_vector;

%Remove unphysiological beats with RR<250ms
remove=find(diff(FPT(:,6))/samplerate<0.25);
FPT(remove+1,:)=[];

if mute==0
    disp('Done')
end

