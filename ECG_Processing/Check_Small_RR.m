% -------------------------------------------------------
%
%    Check_Small_RR - Remove very short unphysiological RR intervals
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
% [FPT_Checked,removed_index]=Check_Small_RR(FPT,samplerate)
% Remove very short unphysiological RR intervals
%
% This is a function used to look in an FPT for very short unphysiological 
% RR intervals. If an RR interval is below 250ms it is removed.	
%
% Inputs:
%       FPT: FPT table that should be checked
%       samplerate: Sampling frequency of the signal for which the given FPT table
%           was created
%
% Outputs:
%       FPT_Checked: corrected FPT without too short RR intervals
%       removed_index: Index of the removed beats in the original table

%
%
% Example Usage:
%       [FPT_Checked,removed_index]=Check_Small_RR(FPT,250)
%
% Revision history:
%
%

function [FPT_Checked,removed_index]=Check_Small_RR(FPT,samplerate)

%Check for small (non physiological) RR intervals. One of the two R peaks 
%shorter building an RR interval shorter than 250ms is removed 
RR=abs(diff(FPT(:,6)))/samplerate*1000;
RR250pos=find(RR<=250);
remove=[];
FPT_Checked=FPT;
while ~isempty(RR250pos)
    for i=1:length(RR250pos)
        if RR250pos(i)>=3 && RR250pos(i)<=size(FPT_Checked,1)-3
            Rloc=FPT_Checked([RR250pos(i)-2;RR250pos(i)-1;RR250pos(i);RR250pos(i)+1;RR250pos(i)+2;RR250pos(i)+3],6);
            RR1=diff([Rloc(1:2);Rloc(4:6)]);
            RR2=diff([Rloc(1:3);Rloc(5:6)]);
            d1=sum(abs(diff(RR1)));
            d2=sum(abs(diff(RR2)));
            if d1>d2
                remove=cat(1,remove,RR250pos(i));
            else
                remove=cat(1,remove,RR250pos(i)+1);
            end
        elseif RR250pos(i)<3
            remove=cat(1,remove,RR250pos(i)+1);
        elseif RR250pos(i)>size(FPT_Checked,1)-3
            remove=cat(1,remove,RR250pos(i));
        end
    end
    
    FPT_Checked(remove,:)=[];
    remove=[];
    RR=abs(diff(FPT_Checked(:,6)))/samplerate*1000;
    RR250pos=find(RR<=250);
end

removed_index=find(~ismember(FPT(:,6),FPT_Checked(:,6)));
