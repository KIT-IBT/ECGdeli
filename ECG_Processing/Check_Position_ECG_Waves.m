% -------------------------------------------------------
%
%    Check_Position_ECG_Waves  - Function for checking if the ECG waves
%    detected in a signal appear in the right order within one beat.
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
% [FPT_Checked,removed_index]=Check_Position_ECG_Waves(FPT)
% Function for checking if the ECG waves detected in a signal appear in the
% right order within one beat. 
%
% This functions is used to check if the ECG waves detected in a signal
% appear in the right order within one beat and among all beats in the
% signal. The following order is expected: ...-P-QRS-T-P-QRS-T-P...
% All beats detected that do not fulfil this criterion are removed from the
% FPT table.
%
% Inputs:
%       FPT: FPT table containing the positions of all the detected waves in the ECG.  
%
% Outputs:
%       FPT_Checked: FPT table containing the positions of all the detected waves 
%           in the ECG that appear in the right order within one beat and among all 
%           beats in the signal
%       removed_index: Index of the removed beats in the original table
%
%
% Example Usage:
%       [FPT_Checked,removed_index]=Check_Position_ECG_Waves(FPT)
%
% Revision history:
%
%

function [FPT_Checked,removed_index]=Check_Position_ECG_Waves(FPT)

%Check that all ECG waves do not appear earlier or later than 
%they should within a beat
ind_waves=sum(FPT,1)>0 & [false,true,false(1,3),true,false(1,4),true,false(1,2)]; %Check if the P, R an T peaks are present in the table
FPT_Checked=FPT;
removed_index=[];
if sum(ind_waves)>1
    M=FPT(:,ind_waves); %Check if the P, R an T peaks are in the right position
    A=diff(M')';
    remove1=~all(A>0,2);

    %Check that all synchronized ECG waves do not appear earlier or later than 
    %they should among two subsequent beat
    i=1:1:size(M,1)-1;
    remove2=~(M(i,end)<M(i+1,1)); % Check that waves from a given beat do not appear later than the next beat
    remove2=[remove2;false];

    remove=remove1 | remove2;
    FPT_Checked(remove,:)=[]; %Create the checked FPT table with the non removed beats
    removed_index=find(remove); %Inde of the removed beats in the original table
end