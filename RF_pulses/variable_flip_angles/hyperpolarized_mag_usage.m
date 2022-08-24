function [mxy, mz] = hyperpolarized_mag_usage(flip, E1);
% [mxy, mz] = hyperpolarized_mag_usage(flip [, E1])
%
% Determines the transverse and longitudinal magnetization after each RF
% pulse for hyperpolarized magnetization.
%
% INPUTS:
%   flip - flip angles, in radians
%   (optional) E1 = exp(-TR/T1) accounts for T1 decay between
%      acquisitions.  Default is 1.
%
% OUTPUTS:
%   mxy,mz - expected  magnetization with given E1
%
% (c) 2008-2013 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson

if nargin < 2
    E1 = 1;
end

N = length(flip);

mxy = sin(flip).* cumprod(cos([0, flip(1:end-1)])) .* E1.^[0:N-1];
mz = cumprod(cos(flip)) .* E1.^[0:N-1] ;