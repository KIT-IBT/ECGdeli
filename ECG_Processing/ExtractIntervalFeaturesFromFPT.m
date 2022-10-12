% -------------------------------------------------------
%
%    ExtractIntvervalFeaturesFromFPT  - Calculate ECG timing features
%    based on FPT tables
%
%    Ver. 1.0.0
%
%    Created:           Claudia Nagel (12.10.2022)
%    Last modified:     Claudia Nagel (12.10.2022)
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
% [features_leadwise, features_sync] = ExtractIntervalFeaturesFromFPT(FPT, FPT_MultiChannel)
% Calculate timing features in ECG based on FPT table
%
%
% Inputs:
%       FPT_MultiChannel: Fiducial Point Table (synchronized)
%       FPT: Fiducial Point Table (lead-wise)
%
% Outputs:
%       features_leadwise: extracted timing features based on FPT
%       annotations in ECG signal (LxBx7 matrix, L: nbr of leads, B:
%       detected beats in FPT_Cell table, 7 timing features - P wave
%       duration, QRS duration, T wave duration, PQ interval, PR interval,
%       QT interval, RR interval)
%
%       features_sync: synchronized timing features across all leads (Bx7
%       matrix, B: detected beats in FPT_MultiChannel table, 8 synchronized
%       timing features - P wave duration, QRS duration, T wave duration, PQ interval, PR interval,
%       QT interval, QTc interval, RR interval)
%
%
% Example Usage:
%       [Timing_feature_12leads, Timing_feature_sync] = ExtractIntervalFeaturesFromFPT(FPT_Cell, FPT_MultiChannel);
%
% Revision history:
%
%

function [features_leadwise, features_sync] = ExtractIntervalFeaturesFromFPT(FPT, FPT_MultiChannel)

FPT_mat = cat(3,FPT{:}); % reshape FPT table
%% extract features lead by lead
features_leadwise = zeros(size(FPT,1),size(FPT_mat,1), 7);
features_leadwise(:,:,1) = 2.*squeeze(FPT_mat(:,3,:)-FPT_mat(:,1,:))'; %p wave duration
features_leadwise(:,:,2) = 2.*squeeze(FPT_mat(:,8,:)-FPT_mat(:,4,:))'; % qrs duration
features_leadwise(:,:,3) = 2.*squeeze(FPT_mat(:,12,:)-FPT_mat(:,10,:))'; % t wave duration
features_leadwise(:,:,4) = 2.*squeeze(FPT_mat(:,4,:)-FPT_mat(:,1,:))'; % pq interval
features_leadwise(:,:,5) = 2.*squeeze(FPT_mat(:,6,:)-FPT_mat(:,1,:))'; % pr interval
features_leadwise(:,:,6) = 2.*squeeze(FPT_mat(:,12,:)-FPT_mat(:,4,:))'; % qt interval
rr = 2.*squeeze(diff(FPT_mat(:,6,:)))';
rr(:,end+1) = rr(:,end);
features_leadwise(:,:,7) = rr; % rr interval
%% find features synchronized over 12 leads
FPT_mat(FPT_mat<5) = NaN;

features_sync = zeros(size(FPT_mat,1),8);
% synchronized timing features, e.g. for P wave duration, are calculated as the 2nd latest detectable P wave ending
% across all leads minus the 2nd earliest detectable P wave beginning of
% all leads. In this way, possibly wrong annotated waveform on- and offsets
% in single leads won't impair the feature calculation

pwd = 2.*(maxk(squeeze(FPT_mat(:,3,:)),2,2)-mink(squeeze(FPT_mat(:,1,:)),2,2));
features_sync(:,1) = pwd(:,2); % p wave duration
qrsd =  2.*(maxk(squeeze(FPT_mat(:,8,:)),2,2)-mink(squeeze(FPT_mat(:,4,:)),2,2));
features_sync(:,2) = qrsd(:,2); % qrs duration
tdur = 2.*(maxk(squeeze(FPT_mat(:,12,:)),2,2)-mink(squeeze(FPT_mat(:,10,:)),2,2));
features_sync(:,3) = tdur(:,2); % t wave duration
pqi = 2.*(maxk(squeeze(FPT_mat(:,4,:)),2,2)-mink(squeeze(FPT_mat(:,1,:)),2,2));
features_sync(:,4) = pqi(:,2); % PQ interval
pri = 2.*(maxk(squeeze(FPT_mat(:,6,:)),2,2)-mink(squeeze(FPT_mat(:,1,:)),2,2));
features_sync(:,5) = pri(:,2); % PR interval
qti = 2.*(maxk(squeeze(FPT_mat(:,12,:)),2,2)-mink(squeeze(FPT_mat(:,4,:)),2,2));
features_sync(:,6) = qti(:,2); % QT interval
rr_multi = 2.*squeeze(diff(FPT_MultiChannel(:,6)));
rr_multi(end+1,1) = rr_multi(end,1); % RR interval
features_sync(:,7) = 1000.*(features_sync(:,6)./1000+0.154.*(1-0.001.*rr_multi)); % QTc (framingham corrections)
features_sync(:,8) = rr_multi;
end

