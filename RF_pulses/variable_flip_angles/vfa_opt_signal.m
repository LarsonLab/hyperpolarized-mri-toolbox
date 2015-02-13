function [flip, mxy, mz] = vfa_opt_signal(N, E1);
% [flip, mxy, mz] = vfa_opt_signal(N,E1)
%
% Calculates series of flip angles for a maximum summed signal across
% all RF pulses.  Derived in Nagashima, K. Optimum pulse flip angles for multi-scan acquisition of hyperpolarized nmr and mri. J Magn Reson 190, 2 (February 2008), 183–188.
%
% INPUTS:
%   N - number of pulses
%   E1 = exp(-TR/T1) accounts for T1 decay between acquisitions.
%
% OUTPUTS:
%   flip - flip angles
%   mxy,mz - expected  magnetization with given E1
%
% (c) 2008-2013 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson


if nargin < 2 || E1 == 1
    flip = vfa_const_amp(N, pi/2);
    return
end

flip = acos( sqrt((E1^2-E1.^(2.*(N-[1:N]+1))) ./ (1-E1.^(2*(N-[1:N]+1)))) );

[mxy, mz] = hyperpolarized_mag_usage(flip, E1);
