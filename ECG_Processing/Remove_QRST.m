% -------------------------------------------------------
%
%    Remove_QRST  - Remove T wave and QRS complex from the ecg signal
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
% [ sig_killed, perc ] = Remove_QRST(signal,samplerate,FPT)
% Remove T wave and QRS complex from the ecg signal
%
% Remove T wave and QRS complex from the ecg signal
%
% Inputs:
%       signal: ECG signal
%       samplerate: sampling rate of the ECG 
%       FPT: Fiducial Point Table
%
% Outputs:
%       sig_killed: ECG with replaced T wave and QRS complex
%       perc: percentage of an RR interval that is overwritten by the
%           sigmoidal function
%
%
% Example Usage:
%       [ sig_killed, perc ] = Remove_QRST(signal,250,FPT)
%
% Revision history:
%
%


function [ sig_killed, perc ] = Remove_QRST(signal,samplerate,FPT)
% Function replaces known QRS-T-parts with with sigmoidal function.
% perc (=percentage of an RR interval) says how much of the RR-interval is overwritten.

%% calculate the percentage of the cutted RR-Intveral with a template
%The size of the template is based on the size of the median RR interval
L1=round(min(1/3*median(diff(FPT(:,6))),500/1000*samplerate));
L2=round(min(2/3*median(diff(FPT(:,6))),1000/1000*samplerate));
L=L1+L2+1;

%The beats used to create the template are selected first. Beats having 
%strong deviations in their rhythmical properties are not considered. The
%deviations are measured using the poincare plot.
RR=diff(FPT(:,6));
X=[RR(1:end-1),RR(2:end)];
index=1:1:size(X,1);
for i=1:2
    SCORE=(X-[mean(X(index,1))*ones(size(X,1),1),mean(X(index,2))*ones(size(X,1),1)])*1/sqrt(2)*[1,-1;1,1];
    D1=abs(SCORE(:,1));
    D2=abs(SCORE(:,2));
    Thl1=2.5*std(D1);
    Thl2=0.7*std(D2);
    index=D1<Thl1 & D2<Thl2;
    if all(~index) %all index not equal to 0
        break
    end
end
Preselected_Beats_Index=find(index)+1;

%In case too little beats are selected to build a template a new
%approach is taken
if length(Preselected_Beats_Index)<1/3*size(FPT,1)
    if all(~index)
        index=(SCORE(:,1)>=-Thl1 & SCORE(:,2)<=Thl2);
    else
        index=(SCORE(:,1)>=min(SCORE(index,1)) & SCORE(:,2)<=max(SCORE(index,2)));
    end
    Preselected_Beats_Index=find(index)+1;
    L1=round(min(1/3*median(diff(FPT(:,6))),500/1000*samplerate));
    L2=round(min(2/3*median(X(index,2)),1000/1000*samplerate));
    L=L1+L2+1;
end

TemplateLength=L;
posRpeak=L1+1;
perc=(TemplateLength-posRpeak-0.05*samplerate)/TemplateLength;

%find q-on and t-off in FPT
sig_killed=signal;
rpeak=FPT(:,6);
r_before=circshift(rpeak,1);
RRint=rpeak-r_before;
RRint=RRint(2:end);
%n=round(samplerate*10e-3);
%start 0.0150s before detected values for Q
%replace values with linear progression
for i=1:1:length(FPT(:,6))-1
    Remove_L=round(FPT(i,6)-FPT(i,5)+30*10^(-3)*samplerate);
    %Remove_L=round(FPT(i,6)-FPT(i,5));
    Remove_R=round((RRint(i)*perc));
    %except entries with 0
    if FPT(i,6)-Remove_L<=1
        xc=linspace(-6,6,FPT(i,6)+Remove_R);
        c=1./(1+exp(-xc));
        y1=signal(1);
        y2=signal(FPT(i,6)+Remove_R+1);
        sig_killed(1:FPT(i,6)+Remove_R,1)=(y2-y1)*c+y1;
    elseif FPT(i+1,6)+Remove_R>length(signal)
        xc=linspace(-6,6,(length(signal)-FPT(i,6)+Remove_L+1));
        c=1./(1+exp(-xc));
        y1=signal(FPT(i,6)-Remove_L-1);
        y2=signal(end);
        sig_killed(FPT(i,6)-Remove_L:end,1)=(y2-y1)*c+y1;
    else
        xc=linspace(-6,6,(1+Remove_R+Remove_L));
        c=1./(1+exp(-xc));
        y1=signal(FPT(i,6)-Remove_L-1);
        y2=signal(FPT(i,6)+Remove_R+1);
        sig_killed(FPT(i,6)-Remove_L:FPT(i,6)+Remove_R,1)=(y2-y1)*c+y1;
    end
end