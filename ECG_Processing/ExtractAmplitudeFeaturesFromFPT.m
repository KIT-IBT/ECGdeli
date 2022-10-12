% -------------------------------------------------------
%
%    ExtractAmplitudeFeaturesFromFPT  - Calculate ECG amplitude features
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
% [feature_12leads] = ExtractAmplitudeFeaturesFromFPT(FPT, signal)
% Calculate amplitude features in ECG based on FPT 
%
%
% Inputs:
%       signal: ECG signal, NxL matrix (N: samples, L: nbr of leads)
%       FPT: Fiducial Point Table (FPT_Cell)
%
% Outputs:
%       feature_12leads: extracted amplitude features based on FPT
%       annotations in ECG signal (LxBx5 matrix, L: nbr of leads, B:
%       detected beats in FPT_Cell table, 5 amplitude features - P wave, Q
%       peak, R peak, S peak, T wave)
%
%
% Example Usage:
%       [Amplitude_feature_12leads] = ExtractAmplitudeFeaturesFromFPT(FPT_Cell, ecg_filtered_isoline);
%
% Revision history:
%
%

function [feature_12leads] = ExtractAmplitudeFeaturesFromFPT(FPT, signal)
%1st dimension: leads, 2nd dimension: beats, 3rd dimension: features (P
%wave, Q peak, R peak, S peak, T wave)

FPT_mat = cat(3,FPT{:});
%% extract features lead by lead
feature_12leads = zeros(12,size(FPT_mat,1), 5); % initialize feature vector
FPT_mat(FPT_mat>size(signal,1)) = size(signal,1); % correct FPT annotations exceeding the signal length
FPT_mat(FPT_mat<=0) = 1;

for leadID = 1:size(signal,2)
    feature_12leads(leadID,:,1) = signal(FPT_mat(:,2,leadID),leadID); % P wave amplitude
    feature_12leads(leadID,:,2) = signal(FPT_mat(:,5,leadID),leadID); % Q peak amplitude
    feature_12leads(leadID,:,3) = signal(FPT_mat(:,6,leadID),leadID); % R peak amplitude
    feature_12leads(leadID,:,4) = signal(FPT_mat(:,7,leadID),leadID); % S peak amplitude
    feature_12leads(leadID,:,5) = signal(FPT_mat(:,11,leadID),leadID); % T wave amplitude
end

end

