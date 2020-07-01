% -------------------------------------------------------
%
%    Sync_Beats  - Synchronize FPTs comming from different channels of the same signal
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
% [FPT_Synced,FPT_Cell_Synced]=Sync_Beats(FPT_Cell,samplerate)
% Synchronize FPTs comming from different channels of the same signal
%
% This function is used to synchronize FPTs comming from different channels 
% of the same signal. The classes of the QRS complexes are not synchronized.
% Every FPT must be placed in the cell array FPT_Cell.
%
% Inputs:
%       FPT_Cell: FPT Cell Array
%       samplerate: sampling frequency in Hz
%
% Outputs:
%       FPT_Synced: synchronized FPT
%       FPT_Cell_Synced: FPTs of all available channels with deleted beats
%           that were removed by the synchronization
%
% Example Usage:
%       [FPT_Synced,FPT_Cell_Synced]=Sync_Beats(FPT_Cell,250)
%
% Revision history:
%
%

function [FPT_Synced,FPT_Cell_Synced]=Sync_Beats(FPT_Cell,samplerate)

disp('Syncing beats...');

%Number of channels
N_Channels=length(FPT_Cell);

%Check for the size of FPT_Cell
if N_Channels==1
    disp('Only one channel present in structure. No synchronization of channels needed.');
    FPT_Synced=FPT_Cell{1};
    FPT_Cell_Synced=FPT_Cell;
    return
end

%Matrix for the futured Synced FPT
M=[];

%FPT_Cell for the synchronized cells
FPT_Cell_Synced=cell(N_Channels,1);

%Maximal time difference in ms between R peaks from different channels to 
%be accepted as one R peak
Time_Limit=100;

%Threshold for the voting decision
Voting_Threshhold=1/2;

%Every FPT is compared with every other FPT. After each comparison the
%beats that are accepted as QRS complexes are removed form all FPT's
%and placed in the synced FPT. After every iteration of the first loop the 
%amount of beats in every table gets smaller.
for n=1:N_Channels
    
    %Components of comp_mat show if a beat located in the ith row of the
    %nth FPT is also located on other channels. Every colum represents a 
    %channel and the number 1 means the beats is also present in the jth
    %table
    comp_mat=zeros(size(FPT_Cell{n},1),N_Channels); 
    pos_mat=comp_mat; %gives position of the QRS complex in each FPT table. Every colum represents a channel.
    samp_point_Pon_mat=comp_mat; %gives sample time of the ith P wave onset of the nth FPT in jth FPT table (every colum represents a channel)
    samp_point_Ppeak_mat=comp_mat; %gives sample time of the P wave peak in each FPT table (every colum represents a channel)
    samp_point_Poff_mat=comp_mat; %gives sample time of the P wave offset in each FPT table (every colum represents a channel)
    samp_point_QRSon_mat=comp_mat; %gives sample time of the ith QRS onset of the nth FPT in jth FPT table (every colum represents a channel)
    samp_point_R_mat=comp_mat; %gives sample time of the ith R peak of the nth FPT in jth FPT table (every colum represents a channel)
    samp_point_Q_mat=comp_mat; %gives sample time of the ith Q peak of the nth FPT in jth FPT table (every colum represents a channel)
    samp_point_S_mat=comp_mat; %gives sample time of the ith S peak of the nth FPT in jth FPT table (every colum represents a channel)
    samp_point_QRSoff_mat=comp_mat; %gives sample time of the ith QRS offset of the nth FPT in jth FPT table (every colum represents a channel)
    samp_point_J_mat=comp_mat; %gives sample time of the J point in each FPT table (every colum represents a channel)
    samp_point_Ton_mat=comp_mat; %gives sample time of the ith T wave onset of the nth FPT in jth FPT table (every colum represents a channel)
    samp_point_Tpeak_mat=comp_mat; %gives sample time of the T wave in each FPT table (every colum represents a channel)
    samp_point_Toff_mat=comp_mat; %gives sample time of the ith T wave onset of the nth FPT in jth FPT table (every colum represents a channel)
    
    %Check all QRS complexes in nth FPT
    for i=1:size(FPT_Cell{n},1)
        
        %Compare the ith QRS complex from the nth FPT with the other FPT's 
        for j=1:N_Channels
            if isempty(FPT_Cell{j})
                continue
            else
                %Timedif: minimal time difference of a QRS complex in the 
                %nth FPT compared to jth FPT 
                %pos: position in the jth FPT
                [Timedif,pos]=min(1000/samplerate*abs(FPT_Cell{n}(i,6)-FPT_Cell{j}(:,6))); 
                if Timedif<=Time_Limit %check if minimal time difference is smaller than Time_Limit
                    comp_mat(i,j)=1; %a QRS complex that was located in the ith row in the nth FPT is also located in the jth FPT 
                    pos_mat(i,j)=pos; %position of the QRS complex in the jth FPT
                    samp_point_Pon_mat(i,j)=FPT_Cell{j}(pos,1); %sample point of the ith P wave onset of the nth FPT in the jth FPT
                    samp_point_Ppeak_mat(i,j)=FPT_Cell{j}(pos,2); %sample point of the ith P wave peak of the nth FPT in the jth FPT
                    samp_point_Poff_mat(i,j)=FPT_Cell{j}(pos,3); %sample point of the ith P wave offset of the nth FPT in the jth FPT
                    samp_point_QRSon_mat(i,j)=FPT_Cell{j}(pos,4); %sample point of the ith ORS onset of the nth FPT in the jth FPT
                    samp_point_Q_mat(i,j)=FPT_Cell{j}(pos,5); %sample point of the ith Q peak of the nth FPT in the jth FPT
                    samp_point_R_mat(i,j)=FPT_Cell{j}(pos,6); %sample point of the ith R peak of the nth FPT in the jth FPT 
                    samp_point_S_mat(i,j)=FPT_Cell{j}(pos,7); %sample point of the ith S peak of the nth FPT in the jth FPT
                    samp_point_QRSoff_mat(i,j)=FPT_Cell{j}(pos,8); %sample point of the ith ORS offset of the nth FPT in the jth FPT
                    samp_point_J_mat(i,j)=FPT_Cell{j}(pos,9); %sample point of the ith J point of the nth FPT in the jth FPT
                    samp_point_Ton_mat(i,j)=FPT_Cell{j}(pos,10); %sample point of the ith T wave onset of the nth FPT in the jth FPT
                    samp_point_Tpeak_mat(i,j)=FPT_Cell{j}(pos,11); %sample point of the ith T wave peak of the nth FPT in the jth FPT
                    samp_point_Toff_mat(i,j)=FPT_Cell{j}(pos,12); %sample point of the ith T wave offset of the nth FPT in the jth FPT
                    
                end
            end
        end
        
        V=sum(comp_mat(i,:))/N_Channels; %Voting for the acceptance of a QRS complex 
        if V>=Voting_Threshhold %QRS complex is accepted if the relative number of votes is bigger than Voting_Threshold
            PonPOS=round(median(samp_point_Pon_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted P wave onset
            PpeakPOS=round(median(samp_point_Ppeak_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted P wave peak
            PoffPOS=round(median(samp_point_Poff_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted P wave offset
            QRSonPOS=round(median(samp_point_QRSon_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted QRS onset
            RPOS=round(median(samp_point_R_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted R peak
            QPOS=round(median(samp_point_Q_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted Q peak
            SPOS=round(median(samp_point_S_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted S peak
            QRSoffPOS=round(median(samp_point_QRSoff_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted QRS offset
            JPOS=round(median(samp_point_J_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted R peak
            TonPOS=round(median(samp_point_Ton_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted T wave onset
            TpeakPOS=round(median(samp_point_Tpeak_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted T wave peak
            ToffPOS=round(median(samp_point_Toff_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted T wave peak
            
            M=cat(1,M,[PonPOS,PpeakPOS,PoffPOS,QRSonPOS,QPOS,RPOS,SPOS,QRSoffPOS,JPOS,TonPOS,TpeakPOS,ToffPOS,0]); %aggregate new fiducial points to the existing ones
            for l=1:N_Channels
                if comp_mat(i,l)
                    FPT_Cell_Synced{l}=cat(1,FPT_Cell_Synced{l},[FPT_Cell{l}(pos_mat(i,l),1:12),0]);
                else
                    FPT_Cell_Synced{l}=cat(1,FPT_Cell_Synced{l},M(end,:));
                end
            end
        end
    end
    for j=1:N_Channels
        FPT_Cell{j}(pos_mat(pos_mat(:,j)>0,j),:)=[]; %Remove beats that were accepted
    end 
end
    
%Check for an empty synchronized table
if isempty(M)
    disp('Warning: Too little QRS complexes were detected. Returning an empty FPT table')
    FPT_Synced=[];
    return
end

%Sort the position of the QRS complexes in the signal and proof uniqueness
%of every QRS complex
FPT_Synced=M;
[~,ind]=sort(M(:,6)); %Sort the position of the QRS complexes in the signal
FPT_Synced(:,1:12)=M(ind,1:12); %Relocate positions of all ECG waves
for i=1:N_Channels
    FPT_Cell_Synced{i}=FPT_Cell_Synced{i}(ind,:); %Relocate positions in each FPT
end

[~,ind]=unique(FPT_Synced(:,6)); %Proof uniqueness
FPT_Synced=FPT_Synced(ind,:); %Relocate positions
for i=1:N_Channels
    FPT_Cell_Synced{i}=FPT_Cell_Synced{i}(ind,:); %Relocate positions
end

%Check for small (non physiological) RR intervals. One of the two R peaks 
%shorter building an RR interval shorter than 250ms is removed
[FPT_Synced,removed_index]=Check_Small_RR(FPT_Synced,samplerate);
for i=1:N_Channels
    FPT_Cell_Synced{i}(removed_index,:)=[];
end

%Check that all synchronized ECG waves do not appear earlier or later than 
%they should within a beat
[FPT_Synced,remove]=Check_Position_ECG_Waves(FPT_Synced);
for i=1:N_Channels
    FPT_Cell_Synced{i}(remove,:)=[];
end

%Check if the synchornized FPT tables is too small to process the ECG 
%signal. In that case return empty FPT tables
if size(FPT_Synced,1)<3 
    disp('Warning: Too little QRS complexes were detected. Returning an empty FPT table')
    FPT_Synced=[];
    for i=1:N_Channels
        FPT_Cell_Synced{i}=[];
    end
end

disp('Done');