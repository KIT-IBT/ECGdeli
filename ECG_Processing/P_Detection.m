% -------------------------------------------------------
%
%    P_Detection  - Function to find the P wave and its boundaries in a ECG signal.
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
% FPT=P_Detection(signal,samplerate,FPT,varargin)
% Function to find the P wave and its boundaries in a ECG signal.
%
% Function to find the P wave and its boundaries in an ECG signal.
% Gives back a matrix FTP with the boundaries of the P waves and its peaks.
% The first column includes the position (sample values) of Pon, the second
% those for Ppeak and the third for Poff. 
%
% Inputs:
%       signal: ecg signal
%       samplerate: sampling rate (in Hz) of the signal
%       FPT: existing FPT table with detected QRS complexes. These shall be found in columns 5 (Q), 6 (R), 7 (S).
%       varargin: option to take the maximum of the ECG signal
%           as P wave peak (0) or just the position of the maximum of the wavelet
%           transform (1). By default, the maximum of the wavelet transform is taken.
%
% Outputs:
%       FPT: updated FPT table with P wave entries
%
%
% Example Usage:
%       FPT=P_Detection( signal,250,FPT,0 )
%
% Revision history:
%
function FPT=P_Detection( signal,samplerate,FPT,varargin )

disp('Detecting P Waves...');
% deciding which maximum is taken as P wave peak
if nargin==3
    take_wt_max=1;
elseif nargin==4
    take_wt_max=varargin{1};
else
    disp('Too many input arguments.')
    return;
end

% filtering the signal with a Gaussian filter
%filteredsignal=svmfftfilter_variable(signal,samplerate,0.1,15,'G');
% define the high-pass and low-pass frequencies

highpass_frequency=1;
lowpass_frequency=15;

% convert to double (is needed for the filtering operation)
if ~isa(signal,'double')
    filteredsignal=double(signal);
    flagsingle=1;
else
    filteredsignal=signal;
    flagsingle=0;
end

% Gaussian high-pass
sigmaf=highpass_frequency;
sigma=samplerate/(2*pi*sigmaf);
length_Gauss=2*round(4*sigma)+1;
h=-fspecial('gaussian',[length_Gauss,1],sigma);
h((length_Gauss+1)/2)=1+h((length_Gauss+1)/2);
for i=1:size(filteredsignal,2)
    filteredsignal(:,i)=conv(filteredsignal(:,i),h,'same');
end

% Gaussian low-pass
sigmaf=lowpass_frequency;
sigma=samplerate/(2*pi*sigmaf);
length_Gauss=round(8*sigma);
h=fspecial('gaussian',[length_Gauss,1],sigma);
for i=1:size(filteredsignal,2)
    filteredsignal(:,i)=conv(filteredsignal(:,i),h,'same');
end

% go back to single
if flagsingle==1
    filteredsignal=single(filteredsignal);
end

%% replacing QRST with sigmoid
[replaced_signal, ant] = Remove_QRST(filteredsignal,samplerate,FPT);

%% preparations for wavelet transformation

% frequency of P wave is about 7 Hz
x = floor(log2(samplerate/2/7));
if x<1
    disp('Very low sampling frequency');
    x=1;
end

% signal has to be extended for a discrete wavelet transformation
% desired length of the signal (power of 2)
if log2(length(replaced_signal))==nextpow2(length(replaced_signal))
    l=2^(nextpow2(length(replaced_signal))+1);
else
    l=2^nextpow2(length(replaced_signal));
end

% what is added on the right
l1=floor((l-length(replaced_signal))/2);

% what is added on the left
l2=l-length(signal)-l1;

% add values
ecg_w = wextend(1,'sp0',replaced_signal,l1,'r');
ecg_w = wextend(1,'sp0',ecg_w,l2,'l');

% definition of the quadratic spline wavelet
g=[0, -2, 2, 0]; %Highpass
h=[1/8, 3/8, 3/8, 1/8]; %Lowpass

%% Discrete wavelet transformation
% using the defined wavelet and the detail coefficients for level x
[~,swd]=swt(ecg_w,x+1,h,g);
%[~,swd]=swt(ecg_w,x+1,'db1');
Dx=swd(x,:)';
Dx=Dx((l2+1:end-l1));

%% Discrete wavelet transformation of the mirrored signal
% using the defined wavelet and the detail coefficients for level x
t=length(ecg_w):-1:1;
[~,swd2]=swt(ecg_w(t),x+1,h,g);
%[~,swd2]=swt(ecg_w(t),x+1,'db1');
Dx2=swd2(x,:)';
Dx2=Dx2(t);
Dx2=Dx2((l2+1:end-l1));

%% Prepare detection
% sum the two wavelet transformations
sum_sig=Dx+Dx2;

% derivative of some signals are needed later
d_sum_sig=diff(sum_sig);
d_filteredsignal=diff(filteredsignal);

% compute the RR intervals
rpeak=FPT(:,6);
r_before=circshift(rpeak,1);
RRint=rpeak-r_before;
RRint=RRint(2:end);

% set on and off points for search intervals
on=rpeak(1:end-1)+round(RRint*ant)-round(35*10^(-3)*samplerate);
off=FPT(:,5)-round(35*10^(-3)*samplerate);
off=off(2:end);
if off(1,1)<=0
    off(1,1)=1;
end

% catch case if on is bigger than off (in case of corrupt FPT)
x=find((on>off)==1);
for i=1:1:length(x)
    while on(x(i))>=off(x(i))
        off(x(i))=off(x(i))+round(10*10^(-3)*samplerate);
    end
end

% use a cell array for the detection; later the characteristic points of a
% P wave will be converted to a matrix.
p=cell(length(off),3);

%% search the peaks of the P wave for each RR interval
for i=1:1:length(off)-1
    % derivative of the signal in a search interval
    d_sum_sig_i=d_sum_sig(rpeak(i):rpeak(i+1));
    t=1:length(d_sum_sig_i)-1;
    % find max
    pos_wt_max=(find((d_sum_sig_i(t)>=0 & d_sum_sig_i(t+1)<0)));
    % find min
    neg_wt_max=(find((d_sum_sig_i(t)<0 & d_sum_sig_i(t+1)>=0)));
    all_wt_max=[pos_wt_max; neg_wt_max];
    % set indentifier to -1 (=min)
    all_wt_max(:,3)=-1; %Min
    % set indentifier to 1 (=max)
    all_wt_max((1:length(pos_wt_max)),3)=1; %Max
    % save results in all_wt_max
    if ~isempty(all_wt_max)
        all_wt_max(:,1)=all_wt_max(:,1)+rpeak(i);
        all_wt_max(:,2)=sum_sig(all_wt_max(:,1));
        all_wt_max(:,4)=on(i);
        all_wt_max(:,5)=off(i);
        all_wt_max(:,6)=abs(all_wt_max(:,2));
        % Sort values
        all_wt_max=sortrows(all_wt_max,[-6 -1]); %descending
    else % case not found anything: identifier=0
        all_wt_max=zeros(1,6);
        all_wt_max(1,4)=on(i);
        all_wt_max(1,5)=off(i);
    end
    % save in cell array
    p{i,2}=all_wt_max;
    p{i,3}=0;
    
    % search the P wave by looking at the peaks of the WT for a P wave
    if ~isempty(all_wt_max) && all_wt_max(1,1)~=0
        % actually search at the highest absolute value of the wavelet
        % transform
        p{i,1}=search_p(replaced_signal, all_wt_max(1,:), samplerate, take_wt_max);
        n=1;
        [j, ~]=size(all_wt_max);
        % if there is a possible position for a p wave continue with a check of
        % the found wave
        if p{i,1}(1,1)~=0
            % look if found p wave has a change in sign in first derivative
            if and(~(d_filteredsignal(p{i,1}(1,1)-1)<0 && d_filteredsignal(p{i,1}(1,1)+1)>=0), ~(d_filteredsignal(p{i,1}(1,1)-1)>0 && d_filteredsignal(p{i,1}(1,1)+1)<=0)) || p{i,1}(1,1)<=on(i)
                % while we have not found a correct maximum or minimum continue
                % search
                while xor((~(d_filteredsignal(p{i,1}(1,1)-1)<0 && d_filteredsignal(p{i,1}(1,1)+1)>=0)), ~(d_filteredsignal(p{i,1}(1,1)-1)>0 && d_filteredsignal(p{i,1}(1,1)+1)<=0)) && p{i,1}(1,1)>=on(i)
                    n=n+1;
                    % break condition to avoid infinity loop
                    if j < n
                        p{i,1}=zeros(1,5);
                        p{i,1}(1,1)=round((off(i)-on(i))/2)+on(i);
                        p{i,1}(1,4)=all_wt_max(1,4);
                        p{i,1}(1,5)=all_wt_max(1,5);
                        %p{i,1}(1,6)=0.5;
                        break;
                    end
                    % redo the search for a P wave
                    p{i,1}=search_p(replaced_signal, all_wt_max(n,:), samplerate, take_wt_max);
                end
            end
        end
        
        %if there is no maximum in the wavelet transform, there is no P wave
        %present. Set half of the search interval as P wave.
    else
        p{i,1}=zeros(1,5);
        p{i,1}(1,1)=round((off(i)-on(i))/2)+on(i);
        p{i,1}(1,4)=all_wt_max(1,4);
        p{i,1}(1,5)=all_wt_max(1,5);
        %p{i,1}(1,6)=0.5;
    end
end

% convert the found peaks to matrix
p=cell2mat(p(:,1));

% create a matrix for pon, ppeak and poff
FPTp=zeros(length(p(:,1)),3);
FPTp(:,2)=p(:,1);


%% Starting delination
% create template for every 50 beats; use the absolute value of the wavelet
% transform
numberofwaves=50;
numberit=ceil(length(p)/numberofwaves);
% define the length of the template on the left and on the right of the
% peak
segmint=round(0.1*samplerate/2*3);

for a=1:1:numberit
    %% Templates for delination
    % build a template for positive or negative waves
    local_FPT=FPT(2:end-1,:);
    if a==numberit
        p_pos=p((numberit-1)*numberofwaves+1:end,:);
        local_FPT=local_FPT((numberit-1)*numberofwaves+1:end,:);
    else
        p_pos=p((a-1)*numberofwaves+1:a*numberofwaves,:);
        local_FPT=local_FPT((a-1)*numberofwaves+1:a*numberofwaves,:);
    end
    
    [sp,~]=size(p_pos);
    if sp>1
        segments_pos=zeros(sp,segmint*2+1); %create matrix for all segments
        
        % start segmentation
        for i=1:1:sp
            if p_pos(i,1)-segmint<=0
                d=abs(p_pos(i,1)-segmint);
                segments_pos(i,(1:d+1))=sum_sig(d+1);
                segments_pos(i,(d+2:end))=(sum_sig(1:p_pos(i,1)+segmint));
            elseif p(i,1)+segmint>length(sum_sig)
                d=p(i,1)+segmint-length(sum_sig);
                segments_pos(i,:)=(sum_sig(p_pos(i,1)-segmint:end));
                segments_pos(i,(d+1:end))=sum_sig(d);
            else
                segments_pos(i,:)=(sum_sig(p_pos(i,1)-segmint:p_pos(i,1)+segmint));
            end
        end
        
        % create the template and discard outliers with the PCA
        [~,score] = pca(abs(segments_pos),'NumComponents',1);
        tsqreduced = mahal(score,score);
        acc=tsqreduced<=(mean(tsqreduced)+std(tsqreduced));
        segments_fil=abs(segments_pos(acc,:));
        meansegments_pos=mean(segments_fil);
        
        template_distribution=cumsum(abs(meansegments_pos))/sum(abs(meansegments_pos));
        pon_template=find(template_distribution>=0.05, 1, 'first');
        poff_template=find(template_distribution>=0.85, 1, 'first');
        chkvar=3;
        
        % catch possible errors: pon after ppeak and poff before ppeak; estimate the width with 50ms
        if pon_template>=segmint-1
            pon_template=segmint-round(50*10^(-3)*samplerate);
        end
        if poff_template<=segmint+1
            poff_template=segmint+round(50*10^(-3)*samplerate);
        end
        
        
        % compute parameters
        area_pwavetemp_plus=cumsum(abs(meansegments_pos(1,segmint+2:poff_template)));
        area_plus=cumsum(abs(meansegments_pos(1,segmint+2:end)));
        area_pwavetemp_minus=cumsum(abs(meansegments_pos(1,pon_template:segmint-1)));
        area_minus=cumsum(abs(meansegments_pos(1,1:segmint-1)));
        thresh_minus=area_pwavetemp_minus(end)/area_minus(end);
        thresh_plus=area_pwavetemp_plus(end)/area_plus(end);
        
        % delineate each wave
        for i=1:1:length(segments_pos(:,1))
            % calculate the needed areas
            area_wave_minus=cumsum(abs(segments_pos(i,1:segmint)));
            area_wave_plus=cumsum(abs(segments_pos(i,(segmint+2):end)));
            
            % find extrema of wavelet transform
            d_template_minus=diff(abs(segments_pos(i,1:segmint)));
            d_template_plus=diff(abs(segments_pos(i,(segmint+2):end)));
            
            %% start with pon (it is the 'minus area'): search extrema
            % gain the estimated position by the area
            area_pon=find(area_wave_minus>=area_wave_minus(end)-(thresh_minus*area_wave_minus(end)),1,'first');
            t=1:length(d_template_minus)-1;
            inflpoint_minus=find((d_template_minus(t)>=0 & d_template_minus(t+1)<0))';
            
            % correct pon to the next extremum, if there is one
            if ~isempty(inflpoint_minus)
                %take the nearest inflection point as pon
                areainfl_pon=dsearchn(inflpoint_minus,area_pon);
                % save the position
                p_pos(i,6)=p_pos(i,1)-segmint+inflpoint_minus(areainfl_pon);
            else
                % otherwise take the estimation based on the area for pon
                p_pos(i,6)=p_pos(i,1)-segmint+area_pon;
            end
            
            %% continue with poff (it is the 'plus area'): search extrema
            % gain the estimated position by the area
            area_poff=find(area_wave_plus>=thresh_plus*area_wave_plus(end),1,'first');
            t=1:length(d_template_plus)-1;
            inflpoint_plus=find((d_template_plus(t)>=0 & d_template_plus(t+1)<0))';
            
            % correct pon to the next estremum, if there is one
            if ~isempty(inflpoint_plus)
                % take the nearest inflection point as poff
                areainfl_poff=dsearchn(inflpoint_plus,area_poff);
                % save the position
                p_pos(i,7)=p_pos(i,1)+inflpoint_plus(areainfl_poff);
            else
                % otherwisse take the estimation based on the area for poff
                p_pos(i,7)=p_pos(i,1)+area_poff;
            end
            
            % catch some possible errors
            % beginning of the signal
            if p_pos(i,6)<1
                p_pos(i,6)=1;
            end
            % end of signal reached
            if p_pos(i,7)>length(signal)
                p_pos(i,6)=length(signal);
            end
            % possible error: poff after Q
            if p_pos(i,7)>local_FPT(i,5)
                p_pos(i,7)=local_FPT(i,5)-round(30*10^(-3)*samplerate);
            end
        end
        % clear the not needed variables to avoid errors
        clear coeff score latent tsquare segments_fil
        % fill the FPT table
        for i=1:1:length(p_pos(:,1))
            indpos=find(FPTp(:,2)==p_pos(i,1));
            FPTp(indpos,1)=p_pos(i,6);
            FPTp(indpos,3)=p_pos(i,7);
            FPTp(indpos,4)=chkvar;
        end
    end
end

%% now set points before the first and after the last detected R point
FPT(2:end-1,1:3)=FPTp(:,1:3);
medpr=round(median(FPT(2:end-1,7)-FPT(2:end-1,2)));
medpron=round(median(FPT(2:end-1,7)-FPT(2:end-1,1)));
medproff=round(median(FPT(2:end-1,7)-FPT(2:end-1,3)));
if FPT(1,7)-medpron < 0
    FPT(1,1)=1;
    FPT(1,2)=2;
    FPT(1,3)=3;
else
    FPT(1,1)=FPT(1,7)-medpron;
    FPT(1,2)=FPT(1,7)-medpr;
    FPT(1,3)=FPT(1,7)-medproff;
end

FPT(end,1)=FPT(end,7)-medpron;
FPT(end,2)=FPT(end,7)-medpr;
FPT(end,3)=FPT(end,7)-medproff;

disp('Done');
end

function [ p ] = search_p( filteredsignal, wt_max, samplerate, take_wt_max )
%% Function to find peaks of P wave
% search in a small search interval around the maximum of the wavelet
% transform (wt_max) for a maximum or minimum value in the signal
% (filteredsignal) with a defined samplerate.
p_width=round(0.1*samplerate/4);
p=zeros(1,5);
p(1,3)=wt_max(1,3);
p(1,4)=wt_max(1,4);
p(1,5)=wt_max(1,5);

% take the maximum of the wavelet transform directly
if take_wt_max==1
    p(1,1) = wt_max(1);
    p(1,2) = filteredsignal(p(1,1));
else
    %maximum
    if wt_max(1)~=0 && sign(wt_max(1,3))==1
        %segmentation max - gives back relative position
        if wt_max(1)-(p_width)<0
            [p(1,2),p(1,1)] = max(filteredsignal(1:wt_max(1)+(p_width),1));
            p(1,1) = wt_max(1)-abs(wt_max(1)-(p_width))+p(1,1)-1;
        elseif wt_max(1)+(p_width)>length(filteredsignal)
            [p(1,2),p(1,1)] = max(filteredsignal(wt_max(1)-(p_width):end,1));
            p(1,1) = wt_max(1)-(p_width)+p(1,1)-1;
        else
            [p(1,2),p(1,1)] = max(filteredsignal(wt_max(1)-(p_width):wt_max(1)+(p_width),1));
            p(1,1) = wt_max(1)-(p_width)+p(1,1)-1;
        end
        %absolute position in interval
        p(1,2) = filteredsignal(p(1,1));
        
        %minimum
    elseif wt_max(1)~=0 && sign(wt_max(1,3))==-1
        %segmentation min - gives back relative position
        if wt_max(1)-(p_width)<0
            [p(1,2),p(1,1)] = min(filteredsignal(1:wt_max(1)+(p_width),1));
            p(1,1) = wt_max(1)-abs(wt_max(1)-(p_width))+p(1,1)-1;
        elseif wt_max(1)+(p_width)>length(filteredsignal)
            [p(1,2),p(1,1)] = min(filteredsignal(wt_max(1)-(p_width):end,1));
            p(1,1) = wt_max(1)-(p_width)+p(1,1)-1;
        else
            [p(1,2),p(1,1)] = min(filteredsignal(wt_max(1)-(p_width):wt_max(1)+(p_width),1));
            p(1,1) = wt_max(1)-(p_width)+p(1,1)-1;
        end
        %set value
        p(1,2) = filteredsignal(p(1,1));
        
    elseif wt_max(1)==0
        p(1,1)=round((p(1,5)-p(1,4))/2)+p(1,4);
    end
end
end

