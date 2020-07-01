% -------------------------------------------------------
%
%    Sync_R_Peaks  - Synchronize detected R peaks 
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
% [FPT_Synced,FPT_Cell_Synced]=Sync_R_Peaks(FPT_Cell,samplerate)
% Synchronize detected R peaks 
%
% Inputs:
%       FPT_Cell: FPT Cell Array
%       samplerate: sampling frequency in Hz
%
% Outputs:
%       FPT_Synced: syncronized FPT
%       FPT_Cell_Synced: FPTs of all available channels with deleted beats
%           that were removed by the synchronization
%
%
% Example Usage:
%       [FPT_Synced,FPT_Cell_Synced]=Sync_R_Peaks(FPT_Cell,250)
%
% Revision history:
%
%


function [FPT_Synced,FPT_Cell_Synced]=Sync_R_Peaks(FPT_Cell,samplerate)

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

%Number of channels
N_Channels=length(FPT_Cell);

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
    samp_point_R_mat=comp_mat; %gives sample time of the ith R peak of the nth FPT in jth FPT table (every colum represents a channel)
    
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
                    samp_point_R_mat(i,j)=FPT_Cell{j}(pos,6); %sample point of the ith R peak of the nth FPT in the jth FPT                    
                end
            end
        end
        
        V=sum(comp_mat(i,:))/N_Channels; %Voting for the acceptance of a QRS complex 
        if V>=Voting_Threshhold %QRS complex is accepted if the relative number of votes is bigger than Voting_Threshold
            RPOS=round(median(samp_point_R_mat(i,comp_mat(i,:)==1))); %recalculate sample point for the accepted R peak
            
            M=cat(1,M,[zeros(1,5),RPOS,zeros(1,7)]); %aggregate new fiducial points to the existing ones
            for l=1:N_Channels
                if comp_mat(i,l)
                    FPT_Cell_Synced{l}=cat(1,FPT_Cell_Synced{l},[zeros(1,5),FPT_Cell{l}(pos_mat(i,l),6),zeros(1,7)]);
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
    
FPT_Synced=zeros(size(M,1),size(FPT_Cell{1},2));
[~,ind]=sort(M(:,6)); %Sort the position of the QRS complexes in the signal
FPT_Synced(:,6)=M(ind,6); %Relocate positions of R peaks
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