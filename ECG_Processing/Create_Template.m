% -------------------------------------------------------
%
%    Create_Template  - Create a template of a whole beat or just the QRS complex
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
% [Template,posRpeak,ampRpeak,Template_matrix]=Create_Template(signal,samplerate,FPT,varargin)
% Create a template of a whole beat or just the QRS complex
%
% Function used to create a template of a whole beat or just the QRS complex. No
% previous classification of the beats is needed. varargin can be used to set which one should be created. Use 'ECG' for the
% whole beat or 'QRS' for the complex. In case varargin is empty a template
% for the whole ecg is generated.
%
% Inputs:
%       signal: A good looking ecg signal
%       samplerate: Sample Frequency of the ecg signal in Hz
%       FPT: Fidutial Point Table of the ecg signal
%       varargin: type of template to be created ('ECG' or 'QRS')
%
% Outputs:
%       Template: Template
%       posRpeak: Position of R peak
%       ampRpeak: Amplitude of R peak
%       Template_matrix: Matrix of all beats used to create the template
%
%
% Example Usage:
%       [Template,posRpeak,ampRpeak,Template_matrix]=Create_Template(signal,250,FPT,'QRS')
%
% Revision history:
%
%

function [Template,posRpeak,ampRpeak,Template_matrix]=Create_Template(signal,samplerate,FPT,varargin)

%Check if the given signal is a vector
if ~isvector(signal)
    error('The input ECG signal must be a vector.');
end

%Check if QRS complex or full ECG is required
if isempty(varargin)
    varargin='ECG';
else
    if length(varargin)>1
        L1=varargin{2}(1,1);
        L2=varargin{2}(1,2);
        flag_length=1;
    else
        L1=round(min(1/3*median(diff(FPT(:,6))),500/1000*samplerate));
        L2=round(min(0.6*median(diff(FPT(:,6))),1000/1000*samplerate));
        flag_length=0;
    end
end


if strcmp(varargin{1},'ECG') || strcmp(varargin{1},'ecg') || strcmp(varargin{1},'QRS') || strcmp(varargin{1},'qrs')
    case_var=0;
    if strcmp(varargin{1},'QRS') || strcmp(varargin{1},'qrs')
        case_var=1;
    end
else
    error('Template type not recognized');
end


if size(signal,2)>1
    if size(signal,1)==1
        signal=signal';
    elseif size(signal,1)<size(signal,2)
        signal=signal(1,:)';
    end
elseif size(signal,2)>size(signal,1)
    signal=signal(:,1);
end

%The size of the template is based on the size of the median RR interval
signal=Isoline_Correction(signal);

L=L1+L2+1;

%The beats used to create the template are selected first. Beats having
%strong deviations in their rhythmical properties are not considered. The
%deviations are measured using the poincare plot.
RR=diff(FPT(:,6));
X=[RR(1:end-1),RR(2:end)];
index=1:1:size(X,1);
for i=1:2
    SCORE=(X-[mean(X(index,1))*ones(size(X,1),1),mean(X(index,2))*ones(size(X,1),1)])*1/sqrt(2)*[1,-1;1,1];
    D1=abs(SCORE(:,1));
    D2=abs(SCORE(:,2));
    Thl1=2.5*std(D1);
    Thl2=0.7*std(D2);
    index=D1<Thl1 & D2<Thl2;
    if all(~index)
        break
    end
end
Preselected_Beats_Index=find(index)+1;

%In case too little beats are selected to build a template a new
%approach is taken
if length(Preselected_Beats_Index)<1/3*size(FPT,1)
    if all(~index)
        index=(SCORE(:,1)>=-Thl1 & SCORE(:,2)<=Thl2);
    else
        index=(SCORE(:,1)>=min(SCORE(index,1)) & SCORE(:,2)<=0);
    end
    Preselected_Beats_Index=find(index)+1;
    if ~flag_length
    L1=round(min(1/3*median(diff(FPT(:,6))),500/1000*samplerate));
    L2=round(min(2/3*median(X(index,2)),1000/1000*samplerate));
    L=L1+L2+1;
    end
    
end

%Check if Preselected_Beats_index is empty. In that case a random signal is
%used as Template
if isempty(Preselected_Beats_Index)
    %L1=round(300/1000*samplerate);
    %L2=round(700/1000*samplerate);
    Template=randn(L1+L2+1,1);
    posRpeak=L1+1;
    ampRpeak=1;
    disp('Warning: Not Enough QRS complexes to build a Template');
    return
end


%If the beats at the beginning and ending borders of ECG signal have
%smaller sizes than the template, they are removed
if FPT(Preselected_Beats_Index(1),6)-L1<1
    ind=find(FPT(Preselected_Beats_Index,6)-L1>=1,1,'first');
    Preselected_Beats_Index=Preselected_Beats_Index(ind:end);
end
if FPT(Preselected_Beats_Index(end),6)+L2>length(signal)
    ind=find(FPT(Preselected_Beats_Index,6)+L2<=length(signal),1,'last');
    Preselected_Beats_Index=Preselected_Beats_Index(1:ind);
end

%Beats having strong deviations in their morphology (maximal and minimal
%amplitudes) are not considered for now
Matrix_ecg=zeros(L,length(Preselected_Beats_Index));
MP=zeros(length(Preselected_Beats_Index),2);
for i=1:length(Preselected_Beats_Index)
    Matrix_ecg(:,i)=signal(FPT(Preselected_Beats_Index(i),6)-L1:FPT(Preselected_Beats_Index(i),6)+L2,1);
    MP(i,:)=[max(Matrix_ecg(:,i)),min(Matrix_ecg(:,i))];
end
Th11=quantile(MP(:,1),0.25);
Th12=quantile(MP(:,1),0.75);
Th21=quantile(MP(:,2),0.25);
Th22=quantile(MP(:,2),0.75);
Selected_Beats_Index=MP(:,1)>=Th11 & MP(:,1)<=Th12 & MP(:,2)>=Th21 & MP(:,2)<=Th22;
if ~sum(Selected_Beats_Index)
    Selected_Beats_Index=true(size(Preselected_Beats_Index));
end
Template=mean(Matrix_ecg(:,Selected_Beats_Index),2);

%Deviation in morphology is measured using the correlation coefficient.
%Beats having a maximum cross covariance lower than 0.75 are not considered.
%The procedure is repeated to ensure robustness
nit=2;
for j=1:nit
    T2=zeros(size(Matrix_ecg,2),1);
    for i=1:length(T2)
        T2(i)=l_operator(Matrix_ecg(:,i),Template);
    end
    Th=quantile(T2,0.15);
    
    if any(T2>=Th)
        Template=mean(Matrix_ecg(:,T2>=Th),2);
    else
        Template=mean(Matrix_ecg(:,Preselected_Beats_Index),2);
    end
    
    if j==nit
        Th=max(quantile(T2,0.2),0.75);
        if any(T2>=Th)
            Template=mean(Matrix_ecg(:,T2>=Th),2);
        end
    end
end

%Other important properties of the Template
Template_matrix=Matrix_ecg(:,T2>=Th);

if ~isempty(Template_matrix)
    %Offset removal
    [~,offset]=Isoline_Correction(Template_matrix(:));
    Template=Template-offset;
    Template_matrix=Template_matrix-offset;
    
    %Position of R peak
    ampRpeak=Template(L1+1);
    posRpeak=L1+1;
    
    %The QRS complex is segmented in case the user wnats only the QRS complex
    if case_var
        n1=round(samplerate*0.08);
        n2=round(samplerate*0.11);
        if L1-n1<1
            n1=L1;
        end
        if n2>L2
            n2=L2-1;
        end
        Template=Template(L1+1-n1:L1+1+n2);
        Template=Template-mean(Template);
        Template_matrix=Template_matrix(L1+1-n1:L1+1+n2,:);
        Template_matrix=Template_matrix-mean(Template);
        ampRpeak=Template(n1+1);
        posRpeak=n1+1;
    end
    
    if ampRpeak~=0
        %Template=Template/ampRpeak; %leave out the normalization
    else
        [~,posRpeak]=max(abs(Template(:,1)));
        ampRpeak=Template(posRpeak);
        %Template=Template/ampRpeak; %leave out the normalization
    end
else
    % If beat selection fails, a template is generated by averaging
    disp('Average template is used...')
    Template=mean(Matrix_ecg(:,Selected_Beats_Index),2);
    %Position of R peak
    ampRpeak=Template(L1+1);
    posRpeak=L1+1;
end