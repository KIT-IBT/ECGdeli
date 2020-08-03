% -------------------------------------------------------
%
%    Annotate_ECG_Multi  - Function used to automatically annotate all ECG waves in a multilead ECG signal. The QRS classification is not performed.
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
% [FPT_MultiChannel,FPT_Cell]=Annotate_ECG_Multi(signal,samplerate,varargin)
% Function used to automatically annotate all ECG waveforms in a multilead ECG signal.
%
% This function can be used to perform the ECG delination in a multilead
% ECG signal. The on- and offset as well as the peak of the P wave, the QRS
% complex and the T wave are entered in the columns of FPT_Cell for each
% lead separately. The different rows represent the detected beats in the
% signal.
%
% Inputs:
%        signal: multilead ecg signal. Every lead should be a standing vector.
%        samplerate: Rate in Hz at which the signal has been sampled. (Sampling Frequency)
%        varargin: This optional variable can have two inputs:
%                   1. The first input must be a string denoting the ECG waves that should be detected. The posible inputs are: 
%                       'all' or 'PQRST', 'PQRS', 'QRST', 'P', 'T' and 'PT'
%                   2. A multichannel FPT table that has been previously given or compute
%
%
% Outputs:
%       FPT_Multichannel: Fidutial Point Table. Table containing the position of the detected QRS complexes and their classification. Voting algorithm is used to synchronize information among channels.
%       FPT_Cell: A cell array containing the FPT table of every channel.
%
%
% Example Usage:
%       [FPT_MultiChannel,FPT_Cell]=Annotate_ECG_Multi(signal,250,'PQRST')
%
% Revision history:
%
%

function [FPT_MultiChannel,FPT_Cell]=Annotate_ECG_Multi(signal,samplerate,varargin)

addpath(genpath('../Filtering'))

if isvector(signal)
    NOC=1;%Number of channels
    if size(signal,1)==1
        signal=signal';
    end
    FPT_Cell=cell(1,1);%Cell array containing the Fidutial Point Table of every lead 
else
    NOC=size(signal,2);%Number of channels
    FPT_Cell=cell(NOC,1);%Cell array containing the Fidutial Point Table of every lead   
end

%Check the optional inputs
if isempty(varargin)
    process_flag='all'; %Default value
    fpt_flag=false;
elseif length(varargin)==1
    process_flag=varargin{1};
    fpt_flag=false;
elseif length(varargin)==2
    process_flag=varargin{1};
    FPT_MultiChannel=varargin{2};
    fpt_flag=true;
    if size(FPT_MultiChannel,2)<13 && size(FPT_MultiChannel,1)>=1
        error('FPT table does not have the correct size')
    elseif ~all(FPT_MultiChannel(:,6)>0)
        error('R peak positions in the FPT are not correct')
    end
else
    error('Number of optional must not be more than 2')
end

if ~(strcmp(process_flag,'all') || strcmp(process_flag,'PQRST') || strcmp(process_flag,'QRS') || strcmp(process_flag,'PQRS') || strcmp(process_flag,'QRST')...
        || strcmp(process_flag,'P') || strcmp(process_flag,'T') || strcmp(process_flag,'PT'))
    error('ECG waves to detect not recognized');
end

flagemptyFPT=0; %Flag for labeling channels with no signal

%QRS detection in the ecg signal in every channel
for i=1:NOC
    if fpt_flag
        %Check QRS complexes in the ecg signal using the given FPT table
        FPT_Cell{i}=Check_R_Peaks_Multi(signal(:,i),samplerate,FPT_MultiChannel);
    else
        %Detect QRS complexes in the ecg signal
        FPT_Cell{i}=QRS_Detection(signal(:,i),samplerate);%Detection of R peaks in every ecg lead
    end    
    if isempty(FPT_Cell{i}) %Check if FPT table is empty
        flagemptyFPT=flagemptyFPT+1;
        continue
    end
end

%Synchronize detected beats
if flagemptyFPT==NOC %Check if all FPT tables are empty
    FPT_MultiChannel=[];%Return empty FPT table
    disp('Warning: Too little QRS komplexes detected. Returning an empty FPT table');
    return
elseif ~fpt_flag
    [FPT_MultiChannel,FPT_Cell]=Sync_Beats(FPT_Cell,samplerate); 
    if isempty (FPT_MultiChannel)
        disp('Warning: Too little QRS komplexes were detected. Returning an empty FPT table');
        return
    end
end

for i=1:NOC
    
    %Detect T Wave
    if (strcmp(process_flag,'all') || strcmp(process_flag,'PQRST') || strcmp(process_flag,'QRST') || strcmp(process_flag,'T') || strcmp(process_flag,'PT'))
        FPT_Cell{i}=T_Detection(signal(:,i),samplerate,FPT_Cell{i});
    end
    
    %Detect P Wave
    if (strcmp(process_flag,'all') || strcmp(process_flag,'PQRST') || strcmp(process_flag,'PQRS')|| strcmp(process_flag,'P') || strcmp(process_flag,'PT'))
        FPT_Cell{i}=P_Detection(signal(:,i),samplerate,FPT_Cell{i});
    end
end

%Resync channels and classes
[FPT_MultiChannel,FPT_Cell]=Sync_Beats(FPT_Cell,samplerate);
