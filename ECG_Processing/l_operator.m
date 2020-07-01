% -------------------------------------------------------
%
%    lop  - calculate the l-operator
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
% lop=l_operator(x,y,varargin)
% calculate the l-operator
%
% Function used to calculate the l-operator between the
% one dimensional signals x and y. The two signals must have the same length.	
%
% Inputs:
%       x: The first signal
%       y: The second signal
%       varargin: Optional variable that can be set to center the two signals
%           being compared. Use 'non_centered' (default) to calcualte l_operator 
%           without centering. Used 'centered' to subtract the mean value of each 
%           signalbefore calcualting the l-operatior. 
%
% Outputs:
%       lop: The outcome for the l-operator based similarity measure
%
%
% Example Usage:
%       lop=l_operator(sig1, sig2)
%
% Revision history:
%
%

function lop=l_operator(x,y,varargin)

%Check for optional input
if isempty(varargin)
    centering=false; %Default: no modifications carried out on the signals
else
    if strcmp(varargin{1},'non_centered');
        centering=false; %Default: no modifications carried out on the signals
    elseif strcmp(varargin{1},'centered');
        centering=true; %The mean value of signal is subtracted from each signal
    else
        error('Wrong optional input. The optional input to use HRT rules must be either ''centered'' or ''non_centered''')
    end
end

%check length of signals
if length(x) ~=length(y)
    error('Signals must have the same length')  
end

%check signals to be standing vectors
if isvector(x)
    if size(x,1)<size(x,2)
        x=x';
    end
else
        error('Signal x must be a standing vector')
end
if isvector(y)
    if size(y,1)<size(y,2)
        y=y';
    end
else
        error('Signal y must be a standing vector')
end

%Remove mean value of each signal
if centering
    x=x-mean(x);
    y=y-mean(y);
end

%l-operator
if all(x==0) && all(y==0)
    lop=1;
else
    lop=(2*x'*y)/(x'*x+y'*y);
end
