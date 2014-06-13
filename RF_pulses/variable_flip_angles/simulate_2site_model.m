function [Mxy, Mz] = simulate_2site_model(Tin, R1, k, flips, TR)
% [Mxy, Mz] = simulate_2site_model(Tin, R1, k, flips, TR)
%    OR 
% [Mxy, Mz] = simulate_2site_model(Mz0, R1, k, flips, TR)
%
% Simulates the magnetization evolution in a two-site exchange model with
% varying flip angles.  Applied disretely as flips(:,1), TR, flips(:,2), TR, etc...
%
% INPUTS:
%   Tin or Mz0 - if Tin is specified, then model will start from all the
%       magnetization in one compound, then evolve for Tin prior to applying
%       flip angles
%       - if Mz0 (1x2) is specified, then model starts from this
%       magnetization distribution
%   R1 - relaxation times for 2 sites (1/s)
%   k - conversion rates between sites (1/s)
%   flips - flip angles (radians)
%   TR - repetition time (s)
%
% OUTPUTS:
%   Mxy,Mz - expected magnetization evolution immediately after flips
%
% (c) 2013-2014 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson

if isempty(R1)
  R1 = [0 0];  % infinite relaxation time
end

if length(R1) == 1
  R1 = [R1 R1];
end

N = size(flips,2);
A = [-R1(1)-k(1) +k(2)
    +k(1) -R1(2)-k(2)];

if length(Tin) == 1
    Mz0 = expm(A*Tin) * [1;0];
else
    Mz0 = Tin;
end

Mxy(1:2,1) = Mz0 .* sin(flips(:,1));
Mz(1:2,1) = Mz0 .* cos(flips(:,1));
for n = 2:N
    Mz_m = expm(A*TR) * Mz(:,n-1);
    Mxy(:,n) =  Mz_m .* sin(flips(:,n));
    Mz(:,n) = Mz_m .* cos(flips(:,n));
end