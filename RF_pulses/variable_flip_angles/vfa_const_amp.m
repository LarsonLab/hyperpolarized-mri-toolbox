function [flip, mxy, mz] = vfa_const_amp(N, flip_end, E1);
% [flip, mxy, mz] = vfa_const_amp(N, flip_end [, E1])
%
% Calculates series of flip angles for a constant amplitude signal in 
% hyperpolarized MR, assuming negligible thermal magnetization.
%
% INPUTS:
%   N - number of pulses
%   flip_end - final pulse flip angle, in radians
%   (optional) E1 = exp(-TR/T1) accounts for T1 decay between
%      acquisitions.  Default is 1.
%
% OUTPUTS:
%   flip - flip angles
%   mxy,mz - expected  magnetization with given E1
%
% (c) 2008-2013 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson


if nargin < 3
    E1 = 1;
end


flip(N) = flip_end;
for n = N-1:-1:1
    flip(n) = atan(E1*sin(flip(n+1)));
end

[mxy, mz] = hyperpolarized_mag_usage(flip, E1);
