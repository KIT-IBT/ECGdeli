function [PMorph, Peaks] = Get_P_Morphology(signal, samplerate, FPT)
% -------------------------------------------------------
%
%    Get_P_Morphology.m  -  Evaluation of P-wave morphology
%    Ver. 1.0.0
%
%    Created:           Fabian Anzlinger (29.06.2020)
%    Last modified:     Fabian Anzlinger (29.06.2020)
%
%    Institute of Biomedical Engineering
%    Karlsruhe Institute of Technology
%
%    http://www.ibt.kit.edu
%
%    Copyright 2000-2022 - All rights reserved.
%
% ------------------------------------------------------
%
% [PMorph, Peaks] = Get_P_Morphology(signal, samplerate, FPT)
%
% This script evaluates P-wave-morphology for every detected P-wave within
% an ECG-signal, whether they are mono- biphasic or M-shaped.
%
% WARNING: This code is not tested for robustness on a wide range of signals.
%          Expect NaN or non-reliable outputs.
%
% Input:
%        signal: ECG-recording (unfiltered signal)
%
%        samplerate: Samplerate of the recording in Hz.
%
%        FPT: FPT-table (cell) of signal calculated by Annotate_ECG_Multi()
%
% Output: PMorph: Table including P-wave-morphology for each detected P-wave
%                Possible contents: 1 = 'monophasic', 2 = 'biphasic' , 3 = 'M-shaped' and a '+' or '-' for positive or negative
%                (beats x lead - table)
%        Peaks:  Table with timestamps/location of peaks
%        depending on detected P-wave shape. 


%% Remove Baseline and correct Isoline
original_signal = signal;
window_length = 0.3;
overlap = 0.3;
[signal,~]=ECG_Baseline_Removal(signal,samplerate,window_length,overlap);
[signal,~,~,~]=Isoline_Correction(signal);
window_length =0.3; overlap = 0.3;
sig_killed = zeros(length(signal(:,1)),length(signal(1,:)));
%Removing QRST-Complex of ECG
for i = 1: length(FPT)
    [ sig_killed(:,i), ~ ] = Remove_QRST(signal(:,i),samplerate,FPT{i,1});
end

%% Lowpass Filter equiripple
lpfrequency = 15;
lpfilt = designfilt('lowpassfir', 'PassbandFrequency', lpfrequency, ...
                     'StopbandFrequency', lpfrequency + 5, 'PassbandRipple', 0.5, ...
                     'StopbandAttenuation', 60, 'SampleRate', 1000, ...
                     'DesignMethod', 'equiripple');

%creating filtered signal:
signalLP = filtfilt(lpfilt, sig_killed); %zero-phase filtering of the sigmoid replaced signal

%--------------------------------------------------------------------------
%% Actual morphology determination
%--------------------------------------------------------------------------
% preallocate arrays
pks = zeros(size(FPT{1,1},1),length(signalLP(1,:)));
invpks = zeros(size(FPT{1,1},1),length(signalLP(1,:)));
negpks = zeros(size(FPT{1,1},1),length(signalLP(1,:)));
PMorph = zeros(size(FPT{1,1},1),length(signal(1,:)));
PadArray = zeros(length(signal(:,1))+25,length(signal(1,:))); 
Peaks = cell(length(FPT),1);

%loop through all leads
for lead = 1:length(signalLP(1,:))
    %loop through all detected P-waves
    for number = 1: size(FPT{lead,1},1)
        % Adjusting window and padding, to bypass possibly wrong detected p-onsets
        if FPT{lead,1}(number,1) > 25   % corresponds to 0.025ms
            xRange = (FPT{lead,1}(number,1)) - 25 : FPT{lead,1}(number,3)+25;
            PadArray(1,lead) = NaN;
        %if first P-wave is wrongly detected, padding    
        else 
            xRange = (1: FPT{lead,1}(number,3) + 25);
            PadArray(1:FPT{lead,1}(number,4) + 25,lead) = padarray(signalLP(1: FPT{lead,1}(number,4),lead),25,signalLP(1,lead), 'pre');
            
        end
        
        %
        if ~isnan(PadArray(1,lead))
            [peaks,locs,width,prom] =  findpeaks(PadArray(1:xRange(end),lead),1:xRange(end),'MinPeakProminence', 0.005,'SortStr','descend');
            negmonosig = - PadArray;
            [~,posidx]= max(prom);
            posloc = locs(posidx);
            [negmonopks,negmonoloc,~,negmonoprom1] =  findpeaks(negmonosig(1:xRange(end),lead),1:xRange(end),'MinPeakProminence', 0.005,'SortStr','descend');
            
            [negmonoprom,negativeindex] = max(negmonoprom1);
            negativeloc = negmonoloc(negativeindex);
           
        else
            [peaks,locs,width,prom] =  findpeaks(signalLP(xRange,lead),xRange,'MinPeakProminence', 0.005,'SortStr','descend');
            negmonosig = - signalLP;
            %searching negative peak
            [negmonopks,negmonoloc,~,negmonoprom1] =  findpeaks(negmonosig(xRange,lead),xRange,'MinPeakProminence', 0.005,'SortStr','descend');
            [negmonoprom, negmonoidx] = max(negmonoprom1);
            negativeloc = negmonoloc(negmonoidx);
        end
        %if one or more peaks are found 
        if ~isempty(peaks) 
            [prom, index] = max(prom); 
            wd = width(index);              
            pks(number,lead)= peaks(index);
            posloc = locs(index);            
        else
            pks(number,lead) = NaN;
            
        end
        if ~isempty(locs)
               invsignalLP = - signalLP;
               %Checking for neg. Peak behind first detected peak to distinguish positive biphasic
               [invpeaks,invloc,invwidth,invprom] = findpeaks(invsignalLP(locs(index):xRange(end),lead),locs(index):xRange(end),'MinPeakProminence', 0.005);
               %Checking for neg Peak behind first detected peak to distinguish negative biphasic
               [negpeaks,negloc,negwidth,negprom] = findpeaks(invsignalLP(xRange(1):locs(index),lead),xRange(1):locs(index),'MinPeakProminence', 0.005); % in front of first positive peak
            if ~isempty(invpeaks)                   %if peak is found in inverted signal go on as done before
                [invprom, invidx] = max(invprom);   
                invpks(number,lead)= invpeaks(invidx);
                invwd = invwidth(invidx);
                invloc = invloc(invidx);
            else
                invpks(number,lead) = NaN;
                invprom = 0;
            end %looking for inverted peak in front of first positive detected peaks
            if ~isempty(negpeaks)                   %if peak is found in inverted signal go on as done before
                [negprom, negidx] = max(negprom);   
                negpks(number,lead)= negpeaks(negidx);
                negwd = negwidth(negidx);
                negloc = negloc(negidx);
            else
                negpks(number,lead) = NaN;
                negprom = 0;
                
            end
        end
        %------------------------------------------------------------------
        % adding PromRel as criteria for biphasic P-wave 
        %------------------------------------------------------------------
        PromRel = NaN;
        WidthRel = 0;
        negPromRel = NaN;
        negWidthRel = 0;
        if ~isnan(pks(number,lead)) && ~isnan(invpks(number,lead)) % looking for filled table or unfilled Table
            PromRel = prom / invprom;
            WidthRel = wd/invwd;            
        end
        
        if ~isnan(pks(number,lead)) && ~isnan(negpks(number,lead)) % looking for filled table or unfilled Table
            negPromRel = prom / negprom;
            negWidthRel = wd/negwd;
        end
        % checking for positive biphasic    
        if  ~isnan(invpks(number,lead)) &&  PromRel > 7/20 && PromRel < 20/7 && WidthRel> 2/3 && WidthRel < 3/2 && invprom > negprom 
            PMorph(number,lead) = 2;
            Peaks{lead,1}(number,1) = posloc;
            Peaks{lead,1}(number,2) = invloc;
        % checking for neg. biphasic P-wave    
        elseif ~isnan(negpks(number,lead)) &&  negPromRel > 7/20 && negPromRel < 20/7 && negWidthRel> 2/3 && negWidthRel < 3/2 && invprom < negprom 
            PMorph(number,lead) = -2; 
            Peaks{lead,1}(number,1) = negloc;
            Peaks{lead,1}(number,2) = posloc;
        %else it is considered as monophasic
        else %distinguish between neg. and pos. monophasic P-wave by comparing prominences of detected pos. and neg. peak (if existent)
            if (isnan(pks(number,lead)) && ~isempty(negmonopks)) | (~isempty(negmonoprom) && negmonoprom/prom > 1 )
                PMorph(number,lead) = -1;
                if exist('negativeloc','var') 
                    Peaks{lead,1}(number,1) = negativeloc;
                    Peaks{lead,1}(number,2) = NaN;
                else
                    Peaks{lead,1}(number,1) = NaN;
                    Peaks{lead,1}(number,2) = NaN;
                end
            else   
                PMorph(number,lead) = 1;
                if ~isempty(posloc) && exist('posloc','var') && ~isnan(pks(number,lead))
                    Peaks{lead,1}(number,1) = posloc;
                    Peaks{lead,1}(number,2) = NaN;
                else
                    Peaks{lead,1}(number,1) = NaN;
                    Peaks{lead,1}(number,2) = NaN;
                end
            end
        end
    end
    
end
%% Display Done and return number of P-waves where no shape detection was possible (most of the times last P-wave due to removal of QRST)
if sum(isnan(pks(number,:))) == length(signal(number,:))
    display(['No shape detection possible for P-wave no. ', num2str(number)]);
    for i = 1: length(signal(number,:))
        PMorph(number,i) = 0;
    end
end
disp('Checking for biphasic P-waves done, continuing with checking for M-shaped P-waves')        
        
   
% -------------------------------------------------------------------------
% checking for M-shaped P-wave 
% -------------------------------------------------------------------------

clear PadArray;
clear invpks;
%replace QRST-Komplex with sigmodial-function
%removing Baseline and correctiong isoline
[Msignal,~]=ECG_Baseline_Removal(original_signal,samplerate,window_length,overlap);
[Msignal,~,~,~]=Isoline_Correction(Msignal);
Msig_killed = zeros(length(Msignal(:,1)),length(Msignal(1,:)));
for i = 1: length(FPT)
    [ Msig_killed(:,i), ~ ] = Remove_QRST(Msignal(:,i),samplerate,FPT{i,1});
end
%--------------------------------------------------------------------------
% Lowpassfilter

lpfrequency = 45;
 Mlpfilt = designfilt('lowpassfir', 'PassbandFrequency', lpfrequency, ...
                     'StopbandFrequency', lpfrequency + 5, 'PassbandRipple', 0.5, ...
                     'StopbandAttenuation', 60, 'SampleRate', 1000, ...
                     'DesignMethod', 'equiripple');
MsignalLP = filtfilt(Mlpfilt, Msig_killed); %zero-phase filtering
negMsig = -MsignalLP;
%--------------------------------------------------------------------------
PadArray = zeros(length(Msignal(:,1))+25,length(Msignal(1,:))); 
%loop through all leads
for lead = 1:length(MsignalLP(1,:)) 
    %loop through all detected P-waves
    for number = 1: size(FPT{lead,1},1)
        % Adjusting window and padding, to bypass possibly wrong detected
        % p-onsets
        if FPT{lead,1}(number,1) > 25   % corresponds to 0.02ms
            xRange = (FPT{lead,1}(number,1)) - 25 : FPT{lead,1}(number,3)+25;
            PadArray(1,lead) = NaN;
            
        else %padding
            xRange = (1: FPT{lead,1}(number,4) + 25);
            PadArray(1:FPT{lead,1}(number,4) + 25,lead) = padarray(MsignalLP(1: FPT{lead,1}(number,4),lead),25,MsignalLP(1,lead), 'pre');
            
        end
%--------------------------------------------------------------------------        
% Checking for postive M-Shape 
%--------------------------------------------------------------------------
        %In current PMorph M-shaped P-waves are recognized as monophasic
        %due to filtering frequency
        if PMorph (number,lead ) == 1
            if ~isnan(PadArray(1,lead))
                [peaks,locs,~,~] =  findpeaks(PadArray(1:xRange(end),lead),1:xRange(end),'MinPeakProminence', 'MinPeakDistance', 0.005,'SortStr', 'descend');
            else
                [peaks,locs,~,~] =  findpeaks(MsignalLP(xRange,lead),xRange,'MinPeakProminence', 0.005,'MinPeakDistance', 15,'SortStr', 'descend');
            end
            
            if ~isempty(peaks) && length(peaks) > 1
                
                pks1 = peaks(1);
                index1 = locs(1);
               
                pks2 = peaks(2);
                index2 = locs(2);
            else
                pks1 = NaN;
                pks2 = NaN;
            end

        if exist('index1','var') && exist('index2','var')
            if ~isnan(pks1)  && ~isnan(pks2) && abs(index1-index2) < 45   && abs(index1-index2) > 15 && abs(pks1)/(abs(pks1)+abs(pks2)) < 0.6 && abs(pks1)/(abs(pks1)+abs(pks2)) > 0.45 
                PMorph(number,lead) = 3;
                if index1 < index2
                    Peaks{lead,1}(number,1) = index1;
                    Peaks{lead,1}(number,2) = index2;  
                else
                    Peaks{lead,1}(number,1) = index2;
                    Peaks{lead,1}(number,2) = index1;
                end
            end
        end
       
        
        end
%--------------------------------------------------------------------------
%Checking for neg. M-Shape:
%--------------------------------------------------------------------------
    if PMorph (number,lead ) == -1
            if ~isnan(PadArray(1,lead))
                [peaks,locs,~,~] =  findpeaks(PadArray(1:xRange(end),lead),1:xRange(end),'MinPeakProminence','MinPeakDistance',15, 0.005,'SortStr', 'descend');
            else
                [peaks,locs,~,~] =  findpeaks(negMsig(xRange,lead),xRange,'MinPeakProminence',0.005,'MinPeakDistance',15,'SortStr', 'descend');
            end
            
            if ~isempty(peaks) && length(peaks) > 1
                pks1 = peaks(1);
                index1 = locs(1);

                pks2 = peaks(2);
                index2 = locs(2);

            else
                pks1 = NaN;
                pks2 = NaN;
            end
      
        if exist('index1','var') && exist('index2','var')
            if ~isnan(pks1) && ~isnan(pks2) && abs(index1-index2) < 45   && abs(index1-index2) > 15 && abs(pks1)/(abs(pks1)+abs(pks2)) < 0.6 && abs(pks1)/(abs(pks1)+abs(pks2)) > 0.45 % looking for filled table or unfilled Table
                PMorph(number,lead) = -3;
                if index1 < index2
                    Peaks{lead,1}(number,1) = index1;
                    Peaks{lead,1}(number,2) = index2;  
                else
                    Peaks{lead,1}(number,1) = index2;
                    Peaks{lead,1}(number,2) = index1;
                end
            end
        end
       
%--------------------------------------------------------------------------        
    end

    end
end
%--------------------------------------------------------------------------
% return an estimation of how reliable the waveforms are detected
%--------------------------------------------------------------------------
Confidencelvl = zeros(1,lead);
for lead = 1 : length(FPT(:,1))
    Shape(1,1) = length(find(PMorph(:,lead)==1));
    Shape(1,2) = length(find(PMorph(:,lead)==-1));
    Shape(1,3) = length(find(PMorph(:,lead)==2));
    Shape(1,4) = length(find(PMorph(:,lead)==-2));
    Shape(1,5) = length(find(PMorph(:,lead)==3));
    Shape(1,6) = length(find(PMorph(:,lead)==-3));
    [~, domShape] = max(Shape);
    
    if domShape == 1
        Confidencelvl(1,lead) = Shape(1,1) / length(PMorph(:,lead));
        display(['Lead no. ', num2str(lead), ' is pos. monophasic with a probability of ', num2str(Confidencelvl(1,lead))]);
        
    elseif domShape == 2    
        Confidencelvl(1,lead) = Shape(1,2) / length(PMorph(:,lead));
        display(['Lead no. ', num2str(lead), ' is neg. monophasic with a probability of ', num2str(Confidencelvl(1,lead))]);
        
    elseif domShape == 3
        Confidencelvl(1,lead) = Shape(1,3) / length(PMorph(:,lead));
        display(['Lead no. ', num2str(lead), ' is pos. biphasic with a probability of ', num2str(Confidencelvl(1,lead))]);
    elseif domShape == 4
        Confidencelvl(1,lead) = Shape(1,4) / length(PMorph(:,lead));
        display(['Lead no. ', num2str(lead), ' is neg. biphasic with a probability of ', num2str(Confidencelvl(1,lead))]);  
    elseif domShape == 5
        Confidencelvl(1,lead) = Shape(1,5) / length(PMorph(:,lead));
        display(['Lead no. ', num2str(lead), ' is pos. m-shaped with a probability of ', num2str(Confidencelvl(1,lead))]);
       
    elseif domShape == 6
        Confidencelvl(1,lead) = Shape(1,6) / length(PMorph(:,lead));
        display(['Lead no. ', num2str(lead), ' is neg. m-shaped/notched with a probability of ', num2str(Confidencelvl(1,lead))]);
       
    else
        Confidencelvl(1,lead) = NaN;
    end
   
end

disp('Done');

