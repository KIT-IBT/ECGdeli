% -------------------------------------------------------
%
%   PT_Detection_single_v2  - Detect P and T waves in a signal
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
% [FPT,isoPosition,templateCorr]=PT_Detection_single_v2(signal,samplerate,FPT)
% Detect P and T waves in a signal
%
% Function used to detect the T and P waves in an ECG signal
%
% Inputs:
%       signal: signal containing an ECG
%       samplerate: sample frequency used to measure the signal
%       FPT: Fidutial Point Table containing the position of the R peaks in the
%           signal. The 6th column of FPT is reserved for this purpose.
%
% Outputs:
%       FPT: Fidutial Point Table containing the position of the P and T peaks in the
%           signal. The first, second and third column of FPT is reserved for the P wave (onset, peak, offset). The 10th,
%           11th and 12th columns contain the T wave onset, peak and offset points
%           respectively.
%       isoPosition: isoline positions in the templates
%       templateCorr: correlation of the built templates
%
%
% Example Usage:
%       [FPT,isoPosition,templateCorr]=PT_Detection_single_v2(signal,250,FPT)
%
% Revision history:
%
%

function [FPT,isoPosition,templateCorr]=PT_Detection_single_v2(signal,samplerate,FPT)

disp('Detecting P and T Waves...')

%check if signal is of type double
if ~isa(signal,'double')
    signal=double(signal);
end
if ~isvector(signal)
    error('The input ECG signal must be a vector!');
elseif size(signal,1)==1
    signal=signal';
end

%% Calculate needed signals
org_sig=signal;
signal=ECG_Low_Filter(org_sig,samplerate,20);
[signalPTNF,intSigm,intNotSigm]=Remove_QRS_single(org_sig,samplerate,FPT);
signalPT=ECG_Low_Filter(signalPTNF,samplerate,20);

%% Variables to set
RRintP=0.25;
nBeatsTemplate=30;
critMWshapeT=0.3;
critBiT=0.3;
critBiWT=0.5;

%specimen
PlengthMin=0.15;%in s


%% Calculate Wavelet transform for all needed levels in advance
% x: level in which we find frequency contents around 14 Hz. Here we find
% the most important components for QRS detection.
x=max(floor(log2(samplerate/2/7))+1,1);

% ECG signal is extended to have a length equal to the next power of 2
if log2(length(signal))==nextpow2(length(signal))
    l=2^(nextpow2(length(signal))+1);
else
    l=2^nextpow2(length(signal));
end
l1=floor((l-length(signal))/2);
l2=l-length(signal)-l1;
ecg_w = wextend(1,'sp0',signal,l1,'r');
ecg_w = wextend(1,'sp0',ecg_w,l2,'l');
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


%% Find all extrema and zero crossings in each RR interval and detect isoline points
x=floor(log2(samplerate/2/7))-1;
dxExtrZerXRR=cell(size(FPT,1)-1,7);
sigExtrRR=cell(size(FPT,1)-1,2);
isoPosition=nan(size(FPT,1)-1,1);
for rri=1:1:size(FPT,1)-1
    [~,dxExtrZerXRR{rri,1}]=findpeaks(Dx(FPT(rri,6):FPT(rri+1,6),x)); %Max
    dxExtrZerXRR{rri,1}=dxExtrZerXRR{rri,1}+FPT(rri,6)+1;
    dxExtrZerXRR{rri,1}(dxExtrZerXRR{rri,1}<=intNotSigm(rri,1)|dxExtrZerXRR{rri,1}>=intNotSigm(rri,2))=[];
    [~,dxExtrZerXRR{rri,2}]=findpeaks(-Dx(FPT(rri,6):FPT(rri+1,6),x));
    dxExtrZerXRR{rri,2}=dxExtrZerXRR{rri,2}+FPT(rri,6)+1; %Min
    dxExtrZerXRR{rri,2}(dxExtrZerXRR{rri,2}<=intNotSigm(rri,1)|dxExtrZerXRR{rri,2}>=intNotSigm(rri,2))=[];
    dxExtrZerXRR{rri,3}=find((Dx1(FPT(rri,6)+1:FPT(rri+1,6),x)<0 & Dx1(FPT(rri,6):FPT(rri+1,6)-1,x)>=0) | (Dx1(FPT(rri,6)+1:FPT(rri+1,6),x)>0 & Dx1(FPT(rri,6):FPT(rri+1,6)-1,x)<=0))+FPT(rri,6)+1;
    dxExtrZerXRR{rri,3}(dxExtrZerXRR{rri,3}<=intNotSigm(rri,1)|dxExtrZerXRR{rri,3}>=intNotSigm(rri,2))=[]; %Zero-crossings Dx1
    dxExtrZerXRR{rri,4}=find((Dx2(FPT(rri,6)+1:FPT(rri+1,6),x)<0 & Dx2(FPT(rri,6):FPT(rri+1,6)-1,x)>=0) | (Dx2(FPT(rri,6)+1:FPT(rri+1,6),x)>0 & Dx2(FPT(rri,6):FPT(rri+1,6)-1,x)<=0))+FPT(rri,6)+1;
    dxExtrZerXRR{rri,4}(dxExtrZerXRR{rri,4}<=intNotSigm(rri,1)|dxExtrZerXRR{rri,4}>=intNotSigm(rri,2))=[]; %Zero-crossings Dx2
    dxExtrZerXRR{rri,5}=find((Dx(FPT(rri,6)+1:FPT(rri+1,6),x)<0 & Dx(FPT(rri,6):FPT(rri+1,6)-1,x)>=0) | (Dx(FPT(rri,6)+1:FPT(rri+1,6),x)>0 & Dx(FPT(rri,6):FPT(rri+1,6)-1,x)<=0))+FPT(rri,6)+1;
    dxExtrZerXRR{rri,5}(dxExtrZerXRR{rri,5}<=intNotSigm(rri,1)|dxExtrZerXRR{rri,5}>=intNotSigm(rri,2))=[]; %Zero-crossings Dx
    
    [~,sigExtrRR{rri,1}]=findpeaks(signal(FPT(rri,6):FPT(rri+1,6))); %Max signal
    sigExtrRR{rri,1}=sigExtrRR{rri,1}+FPT(rri,6)+1;
    sigExtrRR{rri,1}(sigExtrRR{rri,1}<=intNotSigm(rri,1) | sigExtrRR{rri,1}>=intNotSigm(rri,2))=[];
    [~,sigExtrRR{rri,2}]=findpeaks(-signal(FPT(rri,6):FPT(rri+1,6))); %Min signal
    sigExtrRR{rri,2}=sigExtrRR{rri,2}+FPT(rri,6)+1;
    sigExtrRR{rri,2}(sigExtrRR{rri,2}<=intNotSigm(rri,1)| sigExtrRR{rri,2}>=intNotSigm(rri,2))=[];
    
    [~,dxExtrZerXRR{rri,6}]=findpeaks(Dx1(FPT(rri,6):FPT(rri+1,6),x)); %Max Dx1
    dxExtrZerXRR{rri,6}=dxExtrZerXRR{rri,6}+FPT(rri,6)+1;
    dxExtrZerXRR{rri,6}(dxExtrZerXRR{rri,6}<=intNotSigm(rri,1)|dxExtrZerXRR{rri,6}>=intNotSigm(rri,2))=[];
    
    [~,dxExtrZerXRR{rri,7}]=findpeaks(-Dx1(FPT(rri,6):FPT(rri+1,6),x));
    dxExtrZerXRR{rri,7}=dxExtrZerXRR{rri,7}+FPT(rri,6)+1; %Min Dx1
    dxExtrZerXRR{rri,7}(dxExtrZerXRR{rri,7}<=intNotSigm(rri,1)|dxExtrZerXRR{rri,7}>=intNotSigm(rri,2))=[];
    if ~isempty(dxExtrZerXRR{rri,5})
        diffIso=dxExtrZerXRR{rri,5}-intNotSigm(rri,2);
        [~,isoPosition(rri,1)]=min(abs(diffIso));
        isoPosition(rri,1)=dxExtrZerXRR{rri,5}(isoPosition(rri,1));
    end
end

% Go through all isoPositions that are NaN and estimate those
correctionIsoPos=find(isnan(isoPosition(:,1)));
for isoP=1:1:sum(isnan(isoPosition(:,1)))
    dIso=nanmean(diff(isoPosition(max(correctionIsoPos(isoP)-6,1):min(correctionIsoPos(isoP)-1,length(isoPosition)),1)));
    isoPosition(correctionIsoPos(isoP))=isoPosition(correctionIsoPos(isoP)-1)+round(dIso);
end


%% Build a template for every nBeatsTemplate beats to detect the polarity of P and T wave; additionally have a look at the shape
nTemplates=floor(length(FPT(:,6))/nBeatsTemplate);
btStart=(1:nBeatsTemplate:nBeatsTemplate*nTemplates)';
RRmean=mean(diff(FPT(:,6)));
TemplatePT=nan(round(RRmean),nTemplates);
Template=nan(round(RRmean),nTemplates);
templateParameters=cell(nTemplates,2); %1: T-Wave 2: P-Wave
templateCorr=zeros(length(FPT(:,6)),1);
for bt=1:1:nTemplates
    TemplatePT(:,bt)=Create_Template(signalPT,samplerate,FPT(btStart(bt):btStart(bt)+29,:),'ECG',[1,round(RRmean)-2]);
    Template(:,bt)=Create_Template(signal,samplerate,FPT(btStart(bt):btStart(bt)+29,:),'ECG',[1,round(RRmean)-2]);
    
    for i=1:1:nBeatsTemplate
        if size(signalPT,1) >= FPT(btStart(bt)+i-1,6)+round(RRmean)-1
            templateCorr((bt-1)*nBeatsTemplate+i,1)=corr(signalPT(FPT(btStart(bt)+i-1,6):FPT(btStart(bt)+i-1,6)+round(RRmean)-1,1),TemplatePT(:,bt));
        end
    end
    
    % Correct isoline of templates
    isoPositionTemplate=round(mean(isoPosition(btStart(bt):min(btStart(bt)+29,size(isoPosition,1)))-FPT(btStart(bt):min(btStart(bt)+29,size(isoPosition,1)),6)));
    %TemplatePT(:,bt)=Isoline_Correction(TemplatePT(:,bt));
    TemplatePT(:,bt)=TemplatePT(:,bt)-TemplatePT(min(isoPositionTemplate,size(TemplatePT(:,bt),1)),bt);
    
    [pksplus,locsplus]=findpeaks(TemplatePT(round(0.01*samplerate):end-round(0.01*samplerate),bt));
    [pksminus,locsminus]=findpeaks(-TemplatePT(round(0.01*samplerate):end-round(0.01*samplerate),bt));
    locsplus=locsplus+round(0.01*samplerate)-1;
    locsminus=locsminus+round(0.01*samplerate)-1;
    
    % P-wave can be found in the last  of the RR interval
    valplus=[];
    valminus=[];
    d=0;
    while isempty(valplus) || isempty(valminus)
        pcandplusTemplate=locsplus>=(8/12-d/12)*RRmean;
        pcandminusTemplate=locsminus>=(8/12-d/12)*RRmean;
        [valplus,posplus]=max(abs(pksplus(pcandplusTemplate)));
        [valminus,posminus]=max(abs(pksminus(pcandminusTemplate)));
        if ~isempty(valminus) && ~isempty(valplus)
            if valplus >= valminus
                tmp=locsplus(pcandplusTemplate);
                templateParameters{bt,2}=[tmp(posplus),+1,d];
            else
                tmp=locsminus(pcandminusTemplate);
                templateParameters{bt,2}=[tmp(posminus),-1,d];
            end
        elseif isempty(valplus) && ~isempty(valminus)
            tmp=locsminus(pcandminusTemplate);
            templateParameters{bt,2}=[tmp(posminus),-1,d];
        elseif ~isempty(valplus) && isempty(valminus)
            tmp=locsplus(pcandplusTemplate);
            templateParameters{bt,2}=[tmp(posplus),+1,d];
        else
            templateParameters{bt,2}=[0,0,0];
        end
        d=d+1;
        if d==8
            break;
        end
    end
    
    % T-wave can be found in the first 60% of the RR interval
    tcandplusTemplate=locsplus<2/3*RRmean & locsplus>1/10*RRmean;
    tcandminusTemplate=locsminus<2/3*RRmean & locsminus>1/10*RRmean;
    [valplus,posplus]=sort(abs(pksplus(tcandplusTemplate)),'descend');
    [valminus,posminus]=sort(abs(pksminus(tcandminusTemplate)),'descend');
    tmpp=locsplus(tcandplusTemplate);
    tmpm=locsminus(tcandminusTemplate);
    if isempty(valplus) && ~isempty(valminus)
        templateParameters{bt,1}=[tmpm(posminus(1)),-1,1,nan];
        if size(valminus,1)>1 && critMWshapeT*valminus(1)<=valminus(2)
            templateParameters{bt,1}(1,3)=5;
        end
    elseif ~isempty(valplus) && isempty(valminus)
        
        templateParameters{bt,1}=[tmpp(posplus(1)),+1,1,nan];
        if size(valplus,1)>1 && critMWshapeT*valplus(1)<=valplus(2)
            templateParameters{bt,1}(1,3)=5;
        end
    elseif isempty(valplus) && isempty(valminus)
        templateParameters{bt,1}=[0,0,0,nan];
    else
        if valplus(1) >= valminus(1)
            templateParameters{bt,1}=[tmpp(posplus(1)),+1,1,nan];
            if length(valplus)>1
                % check for M-shaped waves
                if critMWshapeT*valplus(1)<=valplus(2)
                    templateParameters{bt,1}(1,3)=5;
                else
                    % Check for biphasic waves
                    if abs(valplus(1))*critBiT < abs(valminus(1))
                        if posplus(1) > posminus(1)
                            templateParameters{bt,1}(1,3)=4;
                        elseif posplus(1) < posminus(1)
                            templateParameters{bt,1}(1,3)=3;
                        end
                    end
                end
            end
        else
            templateParameters{bt,1}=[tmpm(posminus(1)),-1,2,nan];
            if length(valminus)>1
                % check for W-shaped waves
                if critMWshapeT*valminus(1)<=valminus(2)
                    templateParameters{bt,1}(1,3)=6;
                else
                    % Check for biphasic waves
                    if abs(valplus(1)) > valminus(1)*critBiT
                        if posplus(1) > posminus(1)
                            templateParameters{bt,1}(1,3)=4;
                        elseif posplus(1) < posminus(1)
                            templateParameters{bt,1}(1,3)=3;
                        end
                    end
                end
            end
        end
    end
    
    % check if T and P wave are on the same positions
    if templateParameters{bt,1}(1,1)>templateParameters{bt,2}(1,1)
        templateParameters{bt,2}(1,1)=nan;
    end
end

% Calculate correlation for last templates
for i=1:1:(length(FPT(:,6))-1-(nBeatsTemplate*nTemplates))
    if FPT(btStart(end)+nBeatsTemplate+i-1,6)+round(RRmean)-1<=length(signalPT)
        templateCorr(nTemplates*nBeatsTemplate+i,1)=corr(signalPT(FPT(btStart(end)+nBeatsTemplate+i-1,6):FPT(btStart(end)+nBeatsTemplate+i-1,6)+round(RRmean)-1,1),TemplatePT(:,end));
    end
end



%% Detect P/T wave peaks
tt=1;
% We allow maximum 2 peaks with T and P waves
Tmax=nan(size(FPT,1)-1,3);
Ton=nan(size(FPT,1)-1,2);
Toff=nan(size(FPT,1)-1,2);
Tpeak=nan(size(FPT,1)-1,1);

ToffCand=cell(size(FPT,1)-1,1);
TonCand=cell(size(FPT,1)-1,1);

Pmax=nan(size(FPT,1)-1,3);
Pon=nan(size(FPT,1)-1,2);
Poff=nan(size(FPT,1)-1,2);
Ppeak=nan(size(FPT,1)-1,1);

pPosGuess=cell2mat(templateParameters(:,2));
pPosQuantile=[];
if sum(pPosGuess(pPosGuess(:,3)==0 & pPosGuess(:,1)~=0))>5
    pPosQuantile=quantile(pPosGuess(pPosGuess(:,3)==0 & pPosGuess(:,1)~=0,1),0.9);
end

pPolarity=cell2mat(templateParameters(:,2));
pPolarity=sign(sum(pPolarity(:,2)));

tPolarity=cell2mat(templateParameters(:,1));
tPolarity=sign(sum(tPolarity(:,2)));

suspiciousPeaks=nan(size(FPT,1)-1,1);

for rri=1:1:size(FPT,1)-1
    templateNo=min(floor(rri/nBeatsTemplate)+1,nTemplates);
    
    % All peaks in the RR interval based on the Wavelet transform
    AllpeakCandidatesPos=([dxExtrZerXRR{rri,1}, signal(dxExtrZerXRR{rri,1})-signal(isoPosition(rri,1))]);
    AllpeakCandidatesNeg=([dxExtrZerXRR{rri,2}, signal(dxExtrZerXRR{rri,2})-signal(isoPosition(rri,1))]);
    
    % T Wave peak candidates
    if ~isempty(AllpeakCandidatesPos)
        TpeakCandidatesPos=sortrows(AllpeakCandidatesPos,2,'descend');
    else
        TpeakCandidatesPos=nan;
    end
    if ~isempty(AllpeakCandidatesNeg)
        TpeakCandidatesNeg=sortrows(AllpeakCandidatesNeg,2,'ascend');
    else
        TpeakCandidatesNeg=nan;
    end
    
    % If we are very sure about P wave position in the template, we can use
    % it!
    if pPosGuess(templateNo,3)==0 && pPosGuess(templateNo,1)~=0
        TpeakCandidatesPos=TpeakCandidatesPos(TpeakCandidatesPos(:,1)<FPT(rri,6)+pPosGuess(templateNo,1)-PlengthMin*samplerate);
        TpeakCandidatesNeg=TpeakCandidatesNeg(TpeakCandidatesNeg(:,1)<FPT(rri,6)+pPosGuess(templateNo,1)-PlengthMin*samplerate);
    elseif (pPosGuess(templateNo,3)>0 || pPosGuess(templateNo,1)==0) && ~isempty(pPosQuantile)
        TpeakCandidatesPos=TpeakCandidatesPos(TpeakCandidatesPos(:,1)<FPT(rri,6)+pPosQuantile-PlengthMin*samplerate);
        TpeakCandidatesNeg=TpeakCandidatesNeg(TpeakCandidatesNeg(:,1)<FPT(rri,6)+pPosQuantile-PlengthMin*samplerate);
    end
    
    % find corresponding minimum peaks in the signal and delete possible double
    % peaks
    if ~isempty(TpeakCandidatesNeg)
        TminCand=[TpeakCandidatesNeg(:,1),nan(size(TpeakCandidatesNeg(:,1),1),4)];
        if ~isempty(sigExtrRR{rri,2}) && ~isempty(TpeakCandidatesNeg)
            [TminCand(:,2),TminCand(:,3)]=dsearchn(sigExtrRR{rri,2},TpeakCandidatesNeg(:,1));
            TminCand(:,4)=sigExtrRR{rri,2}(TminCand(:,2),1);
            uniPeaks=unique(TminCand(:,2));
            % Double peaks are deleted
            for dd=1:1:size(unique(TminCand(:,2)),1)
                idx=find(TminCand(:,2)==uniPeaks(dd));
                [~,pos]=min(TminCand(idx,3));
                delp=TminCand(:,3)>TminCand(idx(pos),3) & TminCand(:,2)==uniPeaks(dd);
                TminCand(delp,:)=[];
            end
            % Wavelet peaks that are far away from another signal peak are deleted!
            TminCand(TminCand(:,3)>0.05*samplerate,:)=[];
            TminCand(:,5)=signal(TminCand(:,4))-signal(isoPosition(rri,1));
        end
    end
    % find corresponding maximum peaks in the signal and delete possible double
    % peaks
    if ~isempty(TpeakCandidatesPos)
        TmaxCand=[TpeakCandidatesPos(:,1),nan(size(TpeakCandidatesPos(:,1),1),4)];
        if ~isempty(sigExtrRR{rri,1})
            [TmaxCand(:,2),TmaxCand(:,3)]=dsearchn(sigExtrRR{rri,1},TpeakCandidatesPos(:,1));
            TmaxCand(:,4)=sigExtrRR{rri,1}(TmaxCand(:,2),1);
            uniPeaks=unique(TmaxCand(:,2));
            % Double peaks are deleted
            for dd=1:1:size(unique(TmaxCand(:,2)),1)
                idx=find(TmaxCand(:,2)==uniPeaks(dd));
                [~,pos]=min(TmaxCand(idx,3));
                delp=TmaxCand(:,3)>TmaxCand(idx(pos),3) & TmaxCand(:,2)==uniPeaks(dd);
                TmaxCand(delp,:)=[];
            end
            % Wavelet peaks that are far away from another signal peak are deleted!
            TmaxCand(TmaxCand(:,3)>0.05*samplerate,:)=[];
            TmaxCand(:,5)=signal(TmaxCand(:,4))-signal(isoPosition(rri,1));
        end
    end
    
        % We now decide with the templates on the polarity and distinguish
        % between negative and positive peaks
        %if (templateParameters{templateNo,1}(1,2)==1 || (templateParameters{templateNo,1}(1,2)==0 && tPolarity==1)) && ~isempty(TpeakCandidatesPos)
        if tPolarity==1 && ~isempty(TpeakCandidatesPos)
            % If we have a noisy signal, it could be that there is nothing left
            if ~isempty(TmaxCand)
                Tmax(rri,1)=TmaxCand(1,4);
                Tmax(rri,3)=2;
                
                % Now check for M shaped waves
                if size(TmaxCand,1)>1 && critMWshapeT*abs(TmaxCand(1,5))<=abs(TmaxCand(2,5))
                    % Go through all peaks hitting this threshold and mark the
                    % last one
                    idxPeaks=critMWshapeT*abs(TmaxCand(1,5))<=abs(TmaxCand(1:end,5));
                    idxPeaks=find([idxPeaks;0]==0,1,'first')-1; % add a zero so we have at least one there!
                    for pk=2:idxPeaks
                        % it is not allowed that the signal in between gets near to
                        % the isoline
                        if all(abs(signal(min(TmaxCand(1:pk,4)):max(TmaxCand(1:pk,4)))-signal(isoPosition(rri,1)))>0.25*abs(TmaxCand(1,5)))
                            Tmax(rri,1)=TmaxCand(1,4);
                            Tmax(rri,2)=nanmax(Tmax(rri,2),TmaxCand(pk,4));
                            Tmax(rri,3)=6;
                        end
                    end
                elseif ~isempty(TminCand)
                    % Check for biphasic waves: there must be a crossing of the
                    % isolione which is checked by the different sign
                    if abs(TmaxCand(1,5))>abs(TminCand(1,5))
                        if critBiT*abs(TmaxCand(1,5)) < abs(TminCand(1,5)) && sign(TminCand(1,5))~=sign(TmaxCand(1,5)) && abs(critBiWT*Dx(TmaxCand(1,1)))<abs(Dx(TminCand(1,1)))
                            Tmax(rri,1)=TmaxCand(1,4);
                            Tmax(rri,2)=TminCand(1,4);
                            Tmax(rri,3)=7; %3+4
                        end
                    end
                end
            else % if the array is empty, we just take the highest Wavelet candidate
                Tmax(rri,1)=TpeakCandidatesPos(1,1);
                Tmax(rri,3)=2;
            end
            
        %lseif (templateParameters{templateNo,1}(1,2)==-1  || (templateParameters{templateNo,1}(1,2)==0 && tPolarity==-1)) && ~isempty(TpeakCandidatesNeg)
        elseif tPolarity==-1 && ~isempty(TpeakCandidatesNeg)

            % If we have a noisy signal, it could be that there is nothing left
            if ~isempty(TminCand)
                Tmax(rri,1)=TminCand(1,4);
                Tmax(rri,3)=2;
                
                % Now check for W shaped waves
                
                if size(TminCand,1)>1 && critMWshapeT*abs(TminCand(1,5))<=abs(TminCand(2,5))
                    % last one
                    idxPeaks=critMWshapeT*abs(TminCand(1,5))<=abs(TminCand(1:end,5));
                    idxPeaks=find([idxPeaks;0]==0,1,'first')-1; % add a zero so we have at least one there!
                    for pk=2:idxPeaks
                        % it is not allowed that the signal in between gets near to
                        % the isoline
                        if all(abs(signal(min(TminCand(1:pk,4)):max(TminCand(1:pk,4)))-signal(isoPosition(rri,1)))>0.25*abs(TminCand(1,5)))
                            Tmax(rri,1)=TminCand(1,4);
                            Tmax(rri,2)=nanmax(Tmax(rri,2),TminCand(pk,4));
                            Tmax(rri,3)=6;
                        end
                    end
                elseif ~isempty(TmaxCand)
                    % Check for biphasic waves
                    if abs(TminCand(1,5))>abs(TmaxCand(1,5))
                        if critBiT*abs(TminCand(1,5)) < abs(TmaxCand(1,5)) && sign(TminCand(1,5))~=sign(TmaxCand(1,5)) && abs(critBiWT*Dx(TminCand(1,1)))<abs(Dx(TmaxCand(1,1)))
                            Tmax(rri,1)=TminCand(1,4);
                            Tmax(rri,2)=TmaxCand(1,4);
                            Tmax(rri,3)=7; %3+4
                        end
                    end
                end
            else
                Tmax(rri,1)=TpeakCandidatesNeg(1,1);
                Tmax(rri,3)=2;
            end
            
            if max(Tmax(rri,1:2))>FPT(rri,6)+0.85*(FPT(rri+1,6)-FPT(rri,6))
                suspiciousPeaks(rri,1)=1;
            end
        else % just take the maximum/minimum in the interval
            if tPolarity==-1
                [~,Tmax(rri,1)]=min(signalPT(round(FPT(rri,6)+0.1*(FPT(rri+1,6)-FPT(rri,6))):round(FPT(rri,6))+0.85*(FPT(rri+1,6)-FPT(rri,6))));
            else
                [~,Tmax(rri,1)]=max(signalPT(round(FPT(rri,6)+0.1*(FPT(rri+1,6)-FPT(rri,6))):round(FPT(rri,6))+0.85*(FPT(rri+1,6)-FPT(rri,6))));
            end
        end
    
    % Set the peak to the last maximum
    Tpeak(rri,1)=nanmax(Tmax(rri,1:2));

    
    % Delineation -> next extremum
    minmaxidx=[2,1];
    if templateParameters{templateNo,1}(1,2)==-1
        minmaxidx=[1,2];
    end
    
    %_______________________%
    % Finding Toff
    if ~isempty(find(dxExtrZerXRR{rri,minmaxidx(1)}>nanmax(Tmax(rri,1:2))+round(0.01*samplerate)))% & dxExtrZerXRR{rri,minmaxidx(1)}<Ppeak(rri,1),1,'first'))
        % Toff has to be after the maximum and before Ppeak - Save all the  possibilities for a check
        ToffCand{rri,1}=find(dxExtrZerXRR{rri,minmaxidx(1)}>nanmax(Tmax(rri,1:2))+0.01*samplerate);% & dxExtrZerXRR{rri,minmaxidx(1)}<Ppeak(rri,1));
        ToffCand{rri,1}=[ToffCand{rri,1},nan(size(ToffCand{rri,1},1),2),zeros(size(ToffCand{rri,1},1),1)];
        % Toff has to be very near at the isoline
        ToffCand{rri,1}(:,2)=dxExtrZerXRR{rri,minmaxidx(1)}(ToffCand{rri,1}(:,1),1);
        % Toff has to be very near at the isoline
        ToffCand{rri,1}(:,3)=signal(dxExtrZerXRR{rri,minmaxidx(1)}(ToffCand{rri,1}(:,1),1))-signal(isoPosition(rri,1));
        % It is not allowed that there is a maximum between Toff candidate
        % unless the wave is M shaped or biphasic!
        for cd=1:1:size(ToffCand{rri,1},1)
            if ~isempty(find(dxExtrZerXRR{rri,minmaxidx(2)}>(nanmax(Tmax(rri,1:2))+0.01*samplerate) & dxExtrZerXRR{rri,minmaxidx(2)}<ToffCand{rri,1}(cd,2))) && ~isempty(find(dxExtrZerXRR{rri,minmaxidx(1)}>nanmax(Tmax(rri,1:2))+0.01*samplerate & dxExtrZerXRR{rri,minmaxidx(1)}<ToffCand{rri,1}(cd,2)))
                ToffCand{rri,1}(cd,4)=1;
            end
        end
        
        % Now we take the first Tend Candidate that has a significant
        % lower amplitude than the peak
        chosenCand=find(abs(ToffCand{rri,1}(:,3))<0.2*abs(signal(Tpeak(rri,1))-signal(isoPosition(rri,1))),1,'first');
        if ~isempty(chosenCand)
            Toff(rri,2)=ToffCand{rri,1}(chosenCand,2);
        else
            Toff(rri,2)=ToffCand{rri,1}(1,2);
        end
    else
        Toff(rri,2)=1;
    end
    
    %_______________________%
    % Finding Ton
    if ~isempty(find(dxExtrZerXRR{rri,minmaxidx(1)}<nanmin(Tmax(rri,1:2))-round(0.01*samplerate)))% & dxExtrZerXRR{rri,minmaxidx(1)}<Ppeak(rri,1),1,'first'))
        % Toff has to be after the maximum and before Ppeak - Save all the  possibilities for a check
        TonCand{rri,1}=find(dxExtrZerXRR{rri,minmaxidx(1)}<nanmin(Tmax(rri,1:2))-0.01*samplerate);% & dxExtrZerXRR{rri,minmaxidx(1)}<Ppeak(rri,1));
        TonCand{rri,1}=[TonCand{rri,1},nan(size(TonCand{rri,1},1),2),zeros(size(TonCand{rri,1},1),1)];
        % Toff has to be very near at the isoline
        TonCand{rri,1}(:,2)=dxExtrZerXRR{rri,minmaxidx(1)}(TonCand{rri,1}(:,1),1);
        % Toff has to be very near at the isoline
        TonCand{rri,1}(:,3)=signal(dxExtrZerXRR{rri,minmaxidx(1)}(TonCand{rri,1}(:,1),1))-signal(isoPosition(rri,1));
        % It is not allowed that there is a maximum between Toff candidate
        % unless the wave is M shaped or biphasic!
        for cd=1:1:size(TonCand{rri,1},1)
            if ~isempty(find(dxExtrZerXRR{rri,minmaxidx(2)}<(nanmin(Tmax(rri,1:2))-0.01*samplerate) & dxExtrZerXRR{rri,minmaxidx(2)}>TonCand{rri,1}(cd,2))) && ~isempty(find(dxExtrZerXRR{rri,minmaxidx(1)}<nanmin(Tmax(rri,1:2))-0.01*samplerate & dxExtrZerXRR{rri,minmaxidx(1)}>TonCand{rri,1}(cd,2)))
                TonCand{rri,1}(cd,4)=1;
            end
        end
        
        % Now we take the first Tend Candidate that has a significant
        % lower amplitude than the peak
        chosenCand=find(abs(TonCand{rri,1}(:,3))<0.2*abs(signal(Tpeak(rri,1))-signal(isoPosition(rri,1))),1,'first');
        if ~isempty(chosenCand)
            Ton(rri,2)=TonCand{rri,1}(chosenCand,2);
        else
            Ton(rri,2)=TonCand{rri,1}(1,2);
        end
    else
        Ton(rri,2)=1;
    end
    %_______________________%
    
    % Delineation -> zero-crossing after last Tmax is Toff
    if ~isempty(find(dxExtrZerXRR{rri,3}<Toff(rri,2) & dxExtrZerXRR{rri,3}>nanmax(Tmax(rri,1:2)),1,'last'))
        Toff(rri,1)=find(dxExtrZerXRR{rri,3}<Toff(rri,2) & dxExtrZerXRR{rri,3}>nanmax(Tmax(rri,1:2)),1,'last');
        Toff(rri,1)=dxExtrZerXRR{rri,3}(Toff(rri,1),1);
    else
        Toff(rri,1)=Toff(rri,2);
    end
    
    if ~isempty(find(dxExtrZerXRR{rri,3}<Ton(rri,2) & dxExtrZerXRR{rri,3}<nanmin(Tmax(rri,1:2)),1,'last'))
        Ton(rri,1)=find(dxExtrZerXRR{rri,3}<Ton(rri,2) & dxExtrZerXRR{rri,3}<nanmin(Tmax(rri,1:2)),1,'last');
        Ton(rri,1)=dxExtrZerXRR{rri,3}(Ton(rri,1),1);
    else
        Ton(rri,1)=Ton(rri,2);
    end
    
    
    
end

%% We have now collected the fiducial points for P and T waves. We now check if the points are valid points
TR=Tpeak-FPT(1:end-1,6);
TRfiltered=medfilt1(TR,5,'truncate');

for rri=1:1:size(FPT,1)-1
    % if there were more candidates, it is worth to check
    if size(ToffCand{rri,1},1)>1
        if abs(TRfiltered(rri,1)-TR(rri,1))>0.005*samplerate % RTend is not fitting into neighbourhood
            idxChosen=find(ToffCand{rri,1}(:,2)==Toff(rri,2));
            [~,idxNearestNeigRTend]=min(abs((ToffCand{rri,1}(:,2)-FPT(rri,6))-TRfiltered(rri,1)));
            if abs(ToffCand{rri,1}(idxChosen,3)) > abs(ToffCand{rri,1}(idxNearestNeigRTend,3))
                Toff(rri,2)=ToffCand{rri,1}(idxNearestNeigRTend,2);
                disp('Hit');
            end
        end
    elseif isempty(ToffCand{rri,1})
        %disp('Hit');
    end
end






%%
%figure; plot(signal); hold on; plot(Tmax(:,1),org_sig(Tmax(:,1)),'ok','MarkerFaceColor','k'); hold on; plot(Tmax(~isnan(Tmax(:,2)),2),org_sig(Tmax(~isnan(Tmax(:,2)),2)),'or','MarkerFaceColor','r'); hold on; plot(isoPosition,org_sig(isoPosition),'xr','MarkerFaceColor','r'); hold on; plot(Dx(:,x)/5);

% figure; plot(signal); hold on; plot(Tmax(:,1),org_sig(Tmax(:,1)),'ok','MarkerFaceColor','k');...
%     hold on; plot(Tmax(~isnan(Tmax(:,2)),2),signal(Tmax(~isnan(Tmax(:,2)),2)),'or','MarkerFaceColor','r'); hold on;...
%     plot(Pmax(:,1),signal(Pmax(:,1)),'og','MarkerFaceColor','g'); hold on; plot(isoPosition,signal(isoPosition),'xr','MarkerFaceColor','r');...
%     hold on; plot(Dx(:,x)/5);
% figure; plot(signal); hold on; plot(Tmax(:,1),org_sig(Tmax(:,1)),'ok','MarkerFaceColor','k');...
%     hold on; plot(Tmax(~isnan(Tmax(:,2)),2),signal(Tmax(~isnan(Tmax(:,2)),2)),'or','MarkerFaceColor','r'); hold on;...
%     plot(isoPosition,signal(isoPosition),'xr','MarkerFaceColor','r');...
%     hold on; plot(Dx(:,x)/5);

nmeth=1;
% figure; plot(signal); hold on; plot(Tpeak(:,1),signal(Tpeak(:,1)),'ok','MarkerFaceColor','k'); hold on; ...
%     plot(Ton(:,nmeth),signal(Ton(:,nmeth)),'ob','MarkerFaceColor','b'); hold on; plot(Toff(:,nmeth),signal(Toff(:,nmeth)),'or','MarkerFaceColor','r'); ...
%     hold on; plot(Dx(:,x-1)*4);
%
% hold on; plot(Ppeak(:,1),signal(Ppeak(:,1)),'ok','MarkerFaceColor','k'); hold on; ...
%     plot(Pon(:,1),signal(Pon(:,1)),'ob','MarkerFaceColor','b'); hold on; plot(Poff(:,1),signal(Poff(:,1)),'or','MarkerFaceColor','r');

%% Check Tpeak and Ppeak for consistency
%figure; plot(Ppeak-FPT(1:end-1,6))


%% Post-processing T Wave
FPT(1:end-1,10)=Ton(:,nmeth);
FPT(1:end-1,11)=Tpeak(:,1);
FPT(1:end-1,12)=Toff(:,nmeth);
FPT(end,10)=length(signal)-2;
FPT(end,11)=length(signal)-1;
FPT(end,12)=length(signal);

%% Post-processing P Wave
% FPT(1:end-1,1)=Pon(:,1);
% FPT(1:end-1,2)=Ppeak(:,1);
% FPT(1:end-1,3)=Poff(:,1);
% FPT(end,1)=1;
% FPT(end,2)=2;
% FPT(end,3)=3;



