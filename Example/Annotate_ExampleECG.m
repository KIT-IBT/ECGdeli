%
%    Annotate_ExampleECG.m  - This script can be used as an example file
%    for annotating a 12 lead ECG
%
%    Ver. 1.0.0
%
%    Created:           Claudia Nagel (30.06.2020)
%    Last modified:     Claudia Nagel (12.10.2022)
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
% 
% This script can be used as an example file for testing and visualizing the functionalities of ECG deli.
% If filters an examplary 12 lead ECG (basline removal, frequency filtering
% and isoline correction). Subsequently, the annotation process is
% performed and the resulting wave boundaries for the P and T wave and the
% QRS complex are visulized with markers on the signal. 
%
% Required Input: Example ECG (arbitrary number of leads, use e.g. the one provided in the Example folder: 'Example/s0273lre_signal.mat')   
%
%
% The example signal is taken from PTB Diagnostic ECG Database [1], available on physionet [2].
% [1] Bousseljot R, Kreiseler D, Schnabel, A. Nutzung der EKG-Signaldatenbank CARDIODAT der PTB über das Internet. Biomedizinische Technik, Band 40, Ergänzungsband 1 (1995) S 317 
% [2] Goldberger A, Amaral L, Glass L, Hausdorff J, Ivanov PC, Mark R, Mietus JE, Moody GB, Peng CK, Stanley HE. PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.
%
%
% Revision history:
%
%
%%
clear 

addpath(genpath('..')); % add parent directory

% load 12 lead ECG signal

load('Example/s0273lre_signal.mat')

ecg = signal;
Fs = 1000; % specify sample rate

% visualize raw input signal
figure; 
ax1 = subplot(1,4,1);
plot(ecg(1:2000,1), 'LineWidth', 2); 
title('Unfiltered ECG Signal Lead I');
xlabel('time in ms'); ylabel('voltage in mV');
hold all; 


%% filtering of the 12-lead ECG

% Remove baseline wander
% usage: [filtered_signal,baseline]=ECG_Baseline_Removal(signal,samplerate,window_length,overlap)
[ecg_filtered_baseline,~] = ECG_Baseline_Removal(ecg,Fs,1,0.5);

% visualizing waveforms
ax2 = subplot(1,4,2);
plot(ecg_filtered_baseline(1:2000,1), 'LineWidth', 2); 
title('Baseline Removal');
xlabel('time in ms'); ylabel('voltage in mV');
linkaxes([ax1, ax2],'xy')

% filter noise frequencies
% frequencies are already optimized for ECG signals (literature values):
% Lowpass: 120 Hz, Highpass: 0.3 Hz, Bandstop (49-51 Hz)
[ecg_filtered_frq] = ECG_High_Low_Filter(ecg,Fs,1,40);
ecg_filtered_frq=Notch_Filter(ecg_filtered_frq,Fs,50,1);

% visualizing waveforms
ax3 = subplot(1,4,3);
plot(ecg_filtered_frq(1:2000,1), 'LineWidth', 2); 
title('Highpass, Lowpass, Bandstop');
xlabel('time in ms'); ylabel('voltage in mV');
linkaxes([ax1, ax2, ax3],'xy')

% isoline correction
% usage: [filteredsignal,offset,frequency_matrix,bins_matrix]=Isoline_Correction(signal,varargin)

[ecg_filtered_isoline,offset,~,~]=Isoline_Correction(ecg_filtered_frq);

ax4 = subplot(1,4,4);
plot(ecg_filtered_isoline(1:2000,1), 'LineWidth', 2); 
title('Isoline Correction');
xlabel('time in ms'); ylabel('voltage in mV');
linkaxes([ax1, ax2, ax3, ax4],'xy')


%% Feature calculation
% produce FPT Table
% usage: [FPT_MultiChannel,FPT_Cell]=Process_ECG_Multi(signal,samplerate,varargin)
[FPT_MultiChannel,FPT_Cell]=Annotate_ECG_Multi(ecg_filtered_isoline,Fs);

% extract FPTs for Channel 1 (Lead I):
FPT_LeadI = FPT_Cell{3,1};

Pwave_samples = reshape(FPT_LeadI(:,1:3), [1,size(FPT_LeadI(:,1:3),1)*size(FPT_LeadI(:,1:3),2)]);
QRS_samples = reshape([FPT_LeadI(:,4),FPT_LeadI(:,6), FPT_LeadI(:,8)] , [1,size(FPT_LeadI(:,1:3),1)*size(FPT_LeadI(:,1:3),2)]);
Twave_samples = reshape(FPT_LeadI(:,10:12), [1,size(FPT_LeadI(:,10:12),1)*size(FPT_LeadI(:,10:12),2)]);

% visualize fiducial points
figure; 
plot(ecg_filtered_isoline(:,3));
hold on; 
scatter(Pwave_samples, ecg_filtered_isoline(Pwave_samples,3), 'g', 'filled');
scatter(QRS_samples, ecg_filtered_isoline(QRS_samples,3), 'r', 'filled');
scatter(Twave_samples, ecg_filtered_isoline(Twave_samples,3), 'b', 'filled');
title('Filtered ECG');
xlabel('samples'); ylabel('voltage');
legend({'ECG signal', 'P wave', 'QRS complex', 'T wave'});

%% Calculate timing and amplitude features from FPT table
[Amplitude_feature_12leads] = ExtractAmplitudeFeaturesFromFPT(FPT_Cell, ecg_filtered_isoline);
[Timing_feature_12leads, Timing_feature_sync] = ExtractIntervalFeaturesFromFPT(FPT_Cell, FPT_MultiChannel);

% visualize R amplitude for lead II and V1
figure; hold all; 
nbrBeats = size(Amplitude_feature_12leads, 2);
scatter(1:1:nbrBeats, Amplitude_feature_12leads(2,:,3), 'r+'); % lead II (entry 2), all beats, R peak amplitude (entry 3)
scatter(1:1:nbrBeats, Amplitude_feature_12leads(7,:,3), 'ro'); % lead V1 (entry 7), all beats, R peak amplitude (entry 3)
xlabel('Beat ID')
ylabel('R peak amplitude in mV') 
legend({'R amplitude in II', 'R amplitude in V1'})

%% Determine P-wave morphology
[PMorph, Peaks] = Get_P_Morphology(ecg, Fs, FPT_Cell);

