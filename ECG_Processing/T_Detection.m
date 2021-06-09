% -------------------------------------------------------
%
%    T_Detection  - Function used to detect the T wave in an ECG signal
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
% [FPT,signal,intervalSigmoid,Dx,isoPosition,T_type,TON,TOFF,Dx1,Dx2] = T_Detection(signal,samplerate,FPT)
% Function used to detect the T wave in an ECG signal
%
% Inputs:
%       signal: ECG signal
%       samplerate: sampling rate in Hz
%       FPT: Fiducial Point Table
%
% Outputs:
%       FPT: Fidutial Point Table containing the position of the T boundaries in the
%           signal. The 10th, 11th and 12th columns contain the T wave onset, peak and offset points respectively 
%       signal: filtered signal
%       intervalSigmoid: PQRS intervals that are replaced by the sigmoidal
%           function
%       Dx: wavelet transform
%       isoPosition: isoline position
%       T_type: polarity of the T wave
%       TON: T wave onset
%       TOFF: T wave offset
%
%
% Example Usage:
%       [FPT,signal,intervalSigmoid,Dx,isoPosition,T_type,TON,TOFF,Dx1,Dx2] = T_Detection(signal,250,FPT)
%
% Revision history:
%
%


function [FPT] = T_Detection(signal,samplerate,FPT)

disp('Detecting T Waves...')

%check if signal is of type double
if ~isa(signal,'double')
    signal=double(signal);
end
if ~isvector(signal)
    error('The input ECG signal must be a vector!');
elseif size(signal,1)==1
    signal=signal';   
end
%% Filter signal and remove P wave and QRS complexes
%QRS complexes are removed from the signal. After removal only T waves are
%present. They are then easier to detect.
signal_prefiltered=signal;
[signal,intervalSigmoid]=Remove_PQRS(signal,samplerate,FPT);
%Low and highpass filter together with baseline removal
signal=ECG_High_Low_Filter(signal,samplerate,0.3,20,'B');
w=750e-3;
[~,baseline]=ECG_Baseline_Removal(signal,samplerate,w,0.75);
w=2;
[~,baseline]=ECG_Baseline_Removal(baseline,samplerate,w,0.75);
signal=signal-baseline;

%% Wavelet tranformation
% x: level in which we find frequency contents around 7 Hz. Here we find
% the most important components for QRS detection.
x=max(floor(log2(samplerate/2/7))+1,1);
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
[~,swd1_all]=swt(ecg_w,x,'rbio3.3');
t=length(ecg_w):-1:1;
[~,swd2_all]=swt(ecg_w(t),x,'rbio3.3');
Dx1=nan(length(signal),x);
Dx2=nan(length(signal),x);

for lvl=1:1:x
    Dx1tmp=swd1_all(lvl,:)';
    Dx1(:,lvl)=Dx1tmp((l2+1:end-l1));
    Dx2tmp=swd2_all(lvl,:)';
    Dx2tmp=Dx2tmp(t);
    Dx2(:,lvl)=Dx2tmp((l2+1:end-l1));
end
Dx=-(Dx1+Dx2);

[pks,locs]=findpeaks(abs(Dx(:,lvl)));
TPOS_vector=zeros(size(FPT,1),1);
%find possible TPEAK positions
RR = diff(FPT(:,6));
next = min(RR)*1.9;

for i=1:size(FPT,1)
    if i==size(FPT,1) || FPT(i+1,6)-FPT(i,6)>=next
        if i==1 || FPT(i,6)-FPT(i-1,6)>=next
            I_ext=find(locs>FPT(i,8)+round(75e-3*samplerate) & locs<(FPT(i,8)+round(0.4*samplerate)));
        else
            I_ext=find(locs>FPT(i,8)+round(75e-3*samplerate) & locs<round(0.5*(FPT(i-1,6)+FPT(i,6))));
        end
    else
        I_ext=find(locs>FPT(i,8)+round(75e-3*samplerate) & locs<round(0.5*(FPT(i,6)+FPT(i+1,6))));
    end
    if isempty(I_ext)
        if i==size(FPT,1) || FPT(i+1,6)-FPT(i,6)>=next
            if i==1 || FPT(i,6)-FPT(i-1,6)>=next
                TPOS_vector(i)= FPT(i,8)+round(75e-3*samplerate);
            elseif i==size(FPT,1) && FPT(i,6)+abs(FPT(i-1,6)-TPOS_vector(i-1))>length(signal)
                TPOS_vector(i) = length(signal)-3;
            else
                TPOS_vector(i) = FPT(i,6)+abs(FPT(i-1,6)-TPOS_vector(i-1));
            end
        else
            TPOS_vector(i)=round(2/3*FPT(i,6)+1/3*FPT(i+1,6));
        end
    else
        [~,indmax]=max(pks(I_ext));
        TPOS_vector(i)=locs(I_ext(indmax));
    end
end

%% Find isoline
x=floor(log2(samplerate/2/7))-1;
isoPosition = nan(size(FPT,1)-1,1);
signum = cell(size(FPT,1),1);
%find zerocrossings of Dx of each beat
for i=2:size(FPT,1)
    if FPT(i,6)-FPT(i-1,6)>=next
         if FPT(i,4)-1*samplerate<1
             signum{i}=sign(Dx(FPT(i,4)-round(0.5*samplerate):FPT(i,4),x));
             sig=1;
         else
            signum{i}=sign(Dx(FPT(i,4)-round(1*samplerate):FPT(i,4),x));
            sig=2;
         end
         for k=1:length(signum{i})-1
             if abs(signum{i}(k)+signum{i}(k+1))==0
                 zerocrossings(k) = 1;
             else
                zerocrossings(k) = 0;
             end
         end  
         if sig==1
            zerocrossings = find(zerocrossings) + FPT(i,4)-round(0.5*samplerate);
         else
            zerocrossings = find(zerocrossings) + FPT(i,4)-round(1*samplerate);
         end
    else
        signum{i} = sign(Dx(FPT(i-1,4):FPT(i,4),x)); 
        for k=1:abs(FPT(i-1,4)-FPT(i,4))
            if abs(signum{i}(k)+signum{i}(k+1))==0
                zerocrossings(k) = 1;
            else
                zerocrossings(k) = 0;
            end
        end
        zerocrossings = find(zerocrossings) + FPT(i-1,4);
    end
    if ~isempty(zerocrossings)
        %find isoline
        diffiso=nan(1,size(zerocrossings,2));            
        for k=1:size(zerocrossings,2)
            diffiso(k) = abs(zerocrossings(k)-FPT(i,4));
        end
        [~,isoPosition(i-1)] = min(diffiso);
        if isnan(isoPosition(i-1)) || isempty(diffiso(k)) || isoPosition(i-1)==0
            isoPosition(i-1)=NaN;
        else
        isoPosition(i-1) = zerocrossings(isoPosition(i-1)); 
        end
    end
end

% Go through all isoPositions that are NaN and estimate those
correctionIsoPos=find(isnan(isoPosition(:,1)));
for isoP=1:1:sum(isnan(isoPosition(:,1)))
    if isnan(isoPosition(1,1)) && isoP == 1
        isoPosition(1)=[];
    else
    dIso=nanmean(diff(isoPosition(max(correctionIsoPos(isoP)-6,1):min(correctionIsoPos(isoP)-1,length(isoPosition)),1)));
    isoPosition(correctionIsoPos(isoP))=isoPosition(correctionIsoPos(isoP)-1)+round(dIso);
    end
end

%% Build a template to find T type
TPEAK_vector=zeros(size(TPOS_vector,1),1); %position T peak

RR = diff(FPT(:,6));
X = [RR(1:end-1),RR(2:end)];
index=1:1:size(X,1);
SCORE=(X-[mean(X(index,1))*ones(size(X,1),1),mean(X(index,2))*ones(size(X,1),1)])*1/sqrt(2)*[1,-1;1,1];
D1=abs(SCORE(:,1));
Thl1=2.5*std(D1);
index=(SCORE(:,1)>=-Thl1 & SCORE(:,2)<=0);
Ind_QRS_normal=find(index)+1;
if length(Ind_QRS_normal)<2
    Ind_QRS_normal=2:1:size(X,1);
else
    Ind_QRS_normal=Ind_QRS_normal(2:end-1);
end
MP=zeros(length(Ind_QRS_normal),2);
SP=MP;

for i=1:length(Ind_QRS_normal)
    L_width=max(1,round(TPOS_vector(i)-FPT(i,6)-25e-3*samplerate));
    if i==length(Ind_QRS_normal) || FPT(i+1,6)-FPT(i,6)>=next
         if i==1 || FPT(i,6)-FPT(i-1,6)>=next
             R_width= round(0.2*samplerate);
         else
             R_width=max(1,round(0.5*(FPT(i,6)-TPOS_vector(i-1))));    
         end  
    else
        R_width=max(1,round(0.5*(FPT(i+1,6)-TPOS_vector(i))));
    end
    if TPOS_vector(i)-L_width<1 || TPOS_vector(i)+R_width>length(signal)
        continue
    else
        MP(i,:)=[max(Dx(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width,lvl)),min(Dx(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width,lvl))];
        SP(i,:)=[max(signal(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width)),min(signal(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width))];
    end
end
Th11=quantile(MP(:,1),0.25);
Th12=quantile(MP(:,1),0.75);
Th21=quantile(MP(:,2),0.25);
Th22=quantile(MP(:,2),0.75);

ThS11=quantile(SP(:,1),0.25);
ThS12=quantile(SP(:,1),0.75);
ThS21=quantile(SP(:,2),0.25);
ThS22=quantile(SP(:,2),0.75);

T_type=sign(median(MP(MP(:,1)>=Th11 & MP(:,1)<=Th12,1))+median(MP(MP(:,2)>=Th21 & MP(:,2)<=Th22,2))+median(SP(SP(:,1)>=ThS11 & SP(:,1)<=ThS12,1))+median(SP(SP(:,2)>=ThS21 & SP(:,2)<=ThS22,2)));

%correct T type
if T_type==-1 && abs(median(MP(MP(:,1)>=Th11 & MP(:,1)<=Th12,1))+median(MP(MP(:,2)>=Th21 & MP(:,2)<=Th22,2))+median(SP(SP(:,1)>=ThS11 & SP(:,1)<=ThS12,1))+median(SP(SP(:,2)>=ThS21 & SP(:,2)<=ThS22,2)))< 0.3*Th11
    if median(MP(MP(:,1)>=Th11 & MP(:,1)<=Th12,1))>=abs(median(MP(MP(:,2)>=Th21 & MP(:,2)<=Th22,2)))
        T_type=1;
    elseif median(MP(MP(:,1)>=Th11 & MP(:,1)<=Th12,1))>=0.9*abs(median(MP(MP(:,2)>=Th21 & MP(:,2)<=Th22,2))) && median(SP(SP(:,1)>=ThS11 & SP(:,1)<=ThS12,1))>=0.6*abs(median(SP(SP(:,2)>=ThS21 & SP(:,2)<=ThS22,2)))
        T_type=1;
    elseif median(MP(MP(:,1)>=Th11 & MP(:,1)<=Th12,1))>=0.8*abs(median(MP(MP(:,2)>=Th21 & MP(:,2)<=Th22,2)))&& median(SP(SP(:,1)>=ThS11 & SP(:,1)<=ThS12,1))>=0.7*abs(median(SP(SP(:,2)>=ThS21 & SP(:,2)<=ThS22,2)))
        T_type=1;
    end
end

%% Find real position of peak of T wave peak 
for i=1:size(FPT,1)
    L_width=max(1,round(TPOS_vector(i)-FPT(i,8)-75e-3*samplerate));
    if i==size(FPT,1) || FPT(i+1,6)-FPT(i,6)>=next
        if i==1 || FPT(i,6)-FPT(i-1,6)>=next
            R_width= round(0.2*samplerate);
        else
            R_width=max(1,round(0.5*(FPT(i,6)-TPOS_vector(i-1))));
        end
    else
        R_width=max(1,round(0.5*(FPT(i+1,6)-TPOS_vector(i))));
    end
    if TPOS_vector(i)-L_width<1 || TPOS_vector(i)+R_width>length(signal)
        TPEAK_vector(i)=TPOS_vector(i);
    else 
        %find all peaks of WT from lvl to lvl-2
        wt_lvl=Dx(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width,lvl);
        [~,I_ext_pos] = findpeaks(wt_lvl);
        [~,I_ext_neg] = findpeaks(-1*wt_lvl);
        wt_lvl1=Dx(TPOS_vector(i)-L_width:TPOS_vector(i)+R_width,lvl-1);
        
        %T wave must be located before next Sigmoid interval
        if i<size(FPT,1)
            if FPT(i+1,6)-FPT(i,6)>next
                I_ext_pos(I_ext_pos+TPOS_vector(i)-L_width>FPT(i,6)+round(0.8*min(RR)))=[];
                I_ext_neg(I_ext_neg+TPOS_vector(i)-L_width>FPT(i,6)+round(0.8*min(RR)))=[];
            else
                I_ext_pos(I_ext_pos+TPOS_vector(i)-L_width>intervalSigmoid(i+1,1)) = [];
                I_ext_neg(I_ext_neg+TPOS_vector(i)-L_width>intervalSigmoid(i+1,1)) = [];
            end
        end
        %if WT has no peak in search area take TPOS
        if isempty(I_ext_pos) && isempty(I_ext_neg)
            TPEAK_vector(i)=TPOS_vector(i);
        else        
            %sort positive and negative peaks
            [~,pos_pos] = sort(wt_lvl(I_ext_pos));
            [~,pos_neg] = sort(abs(wt_lvl(I_ext_neg)));
            %isoline
            if i==1
                iso = signal_prefiltered(isoPosition(i));
            else
                iso = signal_prefiltered(isoPosition(i-1));
            end
            if T_type > 0   
                %obvious peaks
                if length(pos_pos)==1 && isempty(pos_neg)
                    TPEAK_vector(i) = I_ext_pos(pos_pos)+TPOS_vector(i)-L_width;
                    normal=1;
                elseif isempty(pos_pos) && ~isempty(pos_neg)
                    TPEAK_vector(i) = I_ext_neg(pos_neg(end))+TPOS_vector(i)-L_width;
                    normal=0;
                elseif length(pos_pos)>1 && wt_lvl(I_ext_pos(pos_pos(end))) > 5* wt_lvl(I_ext_pos(pos_pos(end-1)))
                    TPEAK_vector(i) = I_ext_pos(pos_pos(end))+TPOS_vector(i)-L_width;
                    normal=1;
                %second biggest positive peak   
                elseif i>1 && length(pos_pos)>1 && I_ext_pos(pos_pos(end))+TPOS_vector(i)-L_width-FPT(i,6) > 1.4* mean(TPEAK_vector(1:i-1)-FPT(1:i-1,6)) && I_ext_pos(pos_pos(end-1))+TPOS_vector(i)-L_width-FPT(i,6) < 1.4* mean(TPEAK_vector(1:i-1)-FPT(1:i-1,6))
                    TPEAK_vector(i) = I_ext_pos(pos_pos(end-1))+TPOS_vector(i)-L_width;
                    normal=1;    
                %negative T-wave
                elseif abs(wt_lvl(I_ext_neg(pos_neg(end)))) > 5* wt_lvl(I_ext_pos(pos_pos(end)))
                    TPEAK_vector(i) = I_ext_neg(pos_neg(end))+TPOS_vector(i)-L_width;
                    normal=0;   
                elseif abs(wt_lvl(I_ext_neg(pos_neg(end)))-iso) > 1.5*abs(wt_lvl(I_ext_pos(pos_pos(end)))-iso) && abs(wt_lvl1(I_ext_neg(pos_neg(end)))-iso) > 5*abs(wt_lvl1(I_ext_pos(pos_pos(end)))-iso)
                    TPEAK_vector(i) = I_ext_neg(pos_neg(end))+TPOS_vector(i)-L_width;
                    normal=0;
                %normal positive T-wave
                else
                    TPEAK_vector(i) = I_ext_pos(pos_pos(end))+TPOS_vector(i)-L_width;
                    normal=1;
                end
                %search peak in time domain
                if TPEAK_vector(i)- round(0.04*samplerate)<1 || TPEAK_vector(i)+ round(0.04*samplerate)>length(signal)
                    continue
                end
                if normal==1
                    [~,tpeak] = max(signal(TPEAK_vector(i)- round(0.04*samplerate):TPEAK_vector(i)+ round(0.04*samplerate)));
                    TPEAK_vector(i) = TPEAK_vector(i) + tpeak - round(0.04*samplerate); 
                else
                    [~,tpeak] = min(signal(TPEAK_vector(i)- round(0.04*samplerate):TPEAK_vector(i)+ round(0.04*samplerate)));
                    TPEAK_vector(i) = TPEAK_vector(i) + tpeak - round(0.04*samplerate); 
                end
            %T type=-1    
            else    
                % obvious peak
                if length(pos_neg)==1 && isempty(pos_pos)
                    TPEAK_vector(i) = I_ext_neg(pos_neg)+TPOS_vector(i)-L_width;
                    normal=1;
                elseif isempty(pos_neg) && ~isempty(pos_pos)
                    TPEAK_vector(i) = I_ext_pos(pos_pos(end))+TPOS_vector(i)-L_width;
                    normal=0;
                elseif length(pos_neg)>1 && abs(wt_lvl(I_ext_neg(pos_neg(end)))) > 5* abs(wt_lvl(I_ext_neg(pos_neg(end-1))))
                    TPEAK_vector(i) = I_ext_neg(pos_neg(end))+TPOS_vector(i)-L_width;
                    normal=1;
                %second biggest negative peak
                elseif i>1 && length(pos_neg)>1 && I_ext_neg(pos_neg(end))+TPOS_vector(i)-L_width-FPT(i,6) > 1.4* mean(TPEAK_vector(1:i-1)-FPT(1:i-1,6)) && I_ext_neg(pos_neg(end-1))+TPOS_vector(i)-L_width-FPT(i,6) < 1.4* mean(TPEAK_vector(1:i-1)-FPT(1:i-1,6))
                    TPEAK_vector(i) = I_ext_neg(pos_neg(end-1))+TPOS_vector(i)-L_width;
                    normal=1;
                %positive T-wave
                elseif wt_lvl(I_ext_pos(pos_pos(end))) > 5*abs(wt_lvl(I_ext_neg(pos_neg(end))))
                    TPEAK_vector(i) = I_ext_pos(pos_pos(end))+TPOS_vector(i)-L_width;
                    normal=0;  
                elseif abs(wt_lvl(I_ext_pos(pos_pos(end)))-iso) > 1.5*abs(wt_lvl(I_ext_neg(pos_neg(end)))-iso) && abs(wt_lvl1(I_ext_pos(pos_pos(end)))-iso) > 5*abs(wt_lvl1(I_ext_neg(pos_neg(end)))-iso)
                    TPEAK_vector(i) = I_ext_pos(pos_pos(end))+TPOS_vector(i)-L_width;
                    normal=0;
                else
                    %normal negative T-wave 
                    TPEAK_vector(i) = I_ext_neg(pos_neg(end))+TPOS_vector(i)-L_width; 
                    normal=1;
                end 
                %search peak in time domain
                if TPEAK_vector(i)- round(0.04*samplerate)<1 || TPEAK_vector(i)+ round(0.04*samplerate)>length(signal)
                    continue
                end
                if normal==1
                    [~,tpeak] = min(signal(TPEAK_vector(i)- round(0.04*samplerate):TPEAK_vector(i)+ round(0.04*samplerate)));
                    TPEAK_vector(i) = TPEAK_vector(i) + tpeak - round(0.04*samplerate); 
                else
                    [~,tpeak] = max(signal(TPEAK_vector(i)- round(0.04*samplerate):TPEAK_vector(i)+ round(0.04*samplerate))); 
                    TPEAK_vector(i) = TPEAK_vector(i) + tpeak - round(0.04*samplerate); 
                end
            end
        end
    end
end

%% Delination of T wave
TON_vector=nan(size(TPEAK_vector,1),1);
TOFF_vector=nan(size(TPEAK_vector,1),1);
threshold_distribution=0.3;

for i=1:size(FPT,1)
    L_width=max(1,round(TPEAK_vector(i)-FPT(i,6)-25e-3*samplerate));
    if i==size(FPT,1) || FPT(i+1,6)-FPT(i,6)>=next
        if i==1 || FPT(i,6)-FPT(i-1,6)>=next
            R_width= round(0.2*samplerate);
        else
            R_width=max(1,round(0.5*(FPT(i,6)-TPEAK_vector(i-1))));
        end
    else
        R_width=max(1,round(0.5*(FPT(i+1,6)-TPEAK_vector(i))));
    end
    if TPEAK_vector(i)-L_width<FPT(i,8)
        L_width = TPEAK_vector(i)-FPT(i,8);
    end
    if i<size(FPT,1) && TPEAK_vector(i)+R_width>intervalSigmoid(i+1,1)
        R_width = intervalSigmoid(i+1,1)-TPEAK_vector(i);
    end
    if i==1
        iso=signal_prefiltered(isoPosition(1));
    else
        iso=signal_prefiltered(isoPosition(i-1));
    end
    if L_width<2 || R_width<2 || TPEAK_vector(i)+R_width>length(Dx)
        TON_vector(i)= FPT(i,8)+round(75e-3*samplerate);
        TOFF_vector(i)=TPEAK_vector(i)+round(75e-3*samplerate);
    else    
        %find all possible peaks for Ton at lvl and lvl-1
        [~,loco_on1]= findpeaks(Dx(TPEAK_vector(i)-L_width:TPEAK_vector(i),lvl));
        [~,loco_on2]= findpeaks(-1*(Dx(TPEAK_vector(i)-L_width:TPEAK_vector(i),lvl)));
        [~,loco_on1a]= findpeaks(Dx(TPEAK_vector(i)-L_width:TPEAK_vector(i),lvl-1));
        [~,loco_on2a]= findpeaks(-1*(Dx(TPEAK_vector(i)-L_width:TPEAK_vector(i),lvl-1)));
        loco_on_lvl = [loco_on1;loco_on2] +TPEAK_vector(i)-L_width;
        loco_on_lvl1 = [loco_on1a;loco_on2a] +TPEAK_vector(i)-L_width;
        %delete WT peak of TPEAK
        loco_on_lvl(find(TPEAK_vector(i)-loco_on_lvl<=round(0.06*samplerate)))=[];
        loco_on_lvl1(find(TPEAK_vector(i)-loco_on_lvl1<=round(0.04*samplerate)))=[];
        %sort WT peaks
        loco_on_lvl=sort(loco_on_lvl);
        loco_on_lvl1=sort(loco_on_lvl1);
        
        %find possible peaks for Toff at lvl and lvl-1
        [~,loco_off1]= findpeaks(Dx(TPEAK_vector(i):TPEAK_vector(i)+R_width,lvl));
        [~,loco_off2]= findpeaks(-1*(Dx(TPEAK_vector(i):TPEAK_vector(i)+R_width,lvl)));
        [~,loco_off1a]= findpeaks(Dx(TPEAK_vector(i):TPEAK_vector(i)+R_width,lvl-1));
        [~,loco_off2a]= findpeaks(-1*(Dx(TPEAK_vector(i):TPEAK_vector(i)+R_width,lvl-1)));
        loco_off_lvl = [loco_off1;loco_off2]+TPEAK_vector(i);
        loco_off_lvl1 = [loco_off1a;loco_off2a]+TPEAK_vector(i);
        %delete WT peak of TPEAK
        loco_off_lvl(find(loco_off_lvl-TPEAK_vector(i)<=round(0.08*samplerate)))=[];
        loco_off_lvl1(find(loco_off_lvl1-TPEAK_vector(i)<=round(0.04*samplerate)))=[];
        %sort WT peaks
        loco_off_lvl=sort(loco_off_lvl);
        loco_off_lvl1=sort(loco_off_lvl1);
        
        %calculate density of WT at lvl-1
        Dy_density=abs(Dx(TPEAK_vector(i)-L_width:TPEAK_vector(i)+R_width,lvl-1));
        Dy_distribution_L=cumsum(Dy_density(1:L_width)/sum(Dy_density(1:L_width)));
        Dy_distribution_R=cumsum(Dy_density(L_width+1:end)/sum(Dy_density(L_width+1:end)));
    
        %Ton
        if ~isempty(find(Dy_distribution_L>=threshold_distribution,1,'first'))
            TON(i) = TPEAK_vector(i)-L_width+find(Dy_distribution_L>=threshold_distribution,1,'first');
        else
            TON(i) = TPEAK_vector(i)-80*(samplerate/1000); % assuming fixed length of 80ms. 
        end
        %no peaks at lvl found
        if isempty(loco_on_lvl)
            if ~isempty(loco_on_lvl1)
                if length(loco_on_lvl1)>1 && (abs(loco_on_lvl1(end-1)-TON(i))<abs(loco_on_lvl1(end)-TON(i)) || abs(signal(loco_on_lvl1(end)))>0.5*abs(signal(TPEAK_vector(i))))
                    TON_vector(i)=loco_on_lvl1(end-1);
                else
                    TON_vector(i)=loco_on_lvl1(end);
                end
            else
                if FPT(i,8)+round(75e-3*samplerate)-TPEAK_vector(i)<50e-3*samplerate
                    TON_vector(i) = TPEAK_vector(i)-L_width;
                else
                    TON_vector(i)=FPT(i,8)+round(75e-3*samplerate);
                end
            end
        %second last peak
        elseif length(loco_on_lvl)>1 && signal(TPEAK_vector(i))<iso && (Dx(loco_on_lvl(end-1),5)>iso || abs(Dx(loco_on_lvl(end-1),5))<abs(10*iso))
            TON_vector(i)=loco_on_lvl(end-1);
        %peak at lvl-1   
        elseif length(loco_on_lvl1)>1 && length(loco_on_lvl1)<4 && abs(loco_on_lvl1(end)-TON(i))<abs(loco_on_lvl(end)-TON(i)) && length(loco_on_lvl)==1 && loco_on_lvl(end)-FPT(i,8)<round(0.06*samplerate)
            TON_vector(i)=loco_on_lvl1(end);
        %last peak at lvl
        else
            TON_vector(i)=loco_on_lvl(end);
        end
    
        %Toff
        if ~isempty(find(Dy_distribution_R>=1-threshold_distribution,1,'first'))
            TOFF(i)=TPEAK_vector(i)-1+find(Dy_distribution_R>=1-threshold_distribution,1,'first');
        else
            TOFF(i) = TPEAK_vector(i)-1+80*(samplerate/1000); % assuming fixed length of 80ms. 
        end
        %no peaks at lvl found
        if isempty(loco_off_lvl)
            if ~isempty(loco_off_lvl1)
                if length(loco_off_lvl1)>1 && (abs(loco_off_lvl1(2)-TOFF(i))<abs(loco_off_lvl1(1)-TOFF(i)) || abs(signal(loco_off_lvl1(1)))>0.5*abs(signal(TPEAK_vector(i))))
                    TOFF_vector(i)=loco_off_lvl1(2);
                else
                    TOFF_vector(i)=loco_off_lvl1(1);
                end
            else
                TOFF_vector(i)=TPEAK_vector(i)+round(75e-3*samplerate);
            end
        %first peak of WT at lvl
        elseif isempty(loco_off_lvl1)
            TOFF_vector(i) = loco_off_lvl(1);
        elseif loco_off_lvl(1)+round(0.04*samplerate)<length(signal) && loco_off_lvl1(1)+round(0.04*samplerate)<length(signal) && abs(loco_off_lvl(1)-TOFF(i))<round(0.02*samplerate) && abs(sum(diff(signal(loco_off_lvl(1):loco_off_lvl(1)+round(0.04*samplerate)))))<3*abs(sum(diff(signal(loco_off_lvl1(1):loco_off_lvl1(1)+round(0.04*samplerate))))) || abs(TPEAK_vector(i)-loco_off_lvl1(1))<0.6*nanmean(TOFF_vector-TPEAK_vector) || abs(TPEAK_vector(i)-loco_off_lvl1(1))>1.4*nanmean(TOFF_vector-TPEAK_vector)
            TOFF_vector(i)=loco_off_lvl(1);
        %first peak of WT at lvl-1
        else
            TOFF_vector(i)=loco_off_lvl1(1);
        end
    end
end
 
%Place local T wave delineation in global FPT
FPT(:,10)=TON_vector;
FPT(:,11)=TPEAK_vector;
FPT(:,12)=TOFF_vector;

%Place the L point. The L point is used to diagnose an ST elevation or 
%depresion. It is located around the middle between QRS offset and T onset
if all(FPT(:,10)~=0) && all(FPT(:,8)~=0) %QRS offset and T onset are given
    k=0.5;
    FPT(:,9)=round(k*FPT(:,10)+(1-k)*FPT(:,8)); %L point is located around the middle between QRS offset and T onset 
elseif all(FPT(:,10)~=0) && all(FPT(:,7)~=0) %S wave and T onset are given
    k=0.55;
    FPT(:,9)=round(k*FPT(:,10)+(1-k)*FPT(:,7)); %L point is located around the middle between QRS offset and T onset 
else %R peak and T onset are given
    k=0.6;
    FPT(:,9)=round(k*FPT(:,10)+(1-k)*FPT(:,6)); %L point is located around the middle between QRS offset and T onset
end

disp('Done');