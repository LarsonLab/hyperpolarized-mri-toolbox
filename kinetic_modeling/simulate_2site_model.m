function [Mxy, Mz] = simulate_2site_model(Tin, R1, k, flips, TR)
% simulate_2site_model - Simulates data using a 2-site kinetic model with no
% input function.
% Starts from M0 = [1 0], then evolves for Tin, [then flip, then TR]-repeat
% Mxy and Mz are computed right after flips
%
% [Mxy, Mz] = simulate_2site_model(Tin, R1, k, flips, TR)
%
% INPUTS
%	Tin - delay until acquisition starts (s)
%   R1 - relaxation rates (1/s)
%   k - forward and reverse (optional) conversion rates (1/s)
%   flips - all flip angles
%   TR - time between flip angles applied (s)
% OUTPUTS
%   Mxy, Mz - resulting magnetization after each flip
%
% Author: Peder E. Z. Larson
%
% (c)2013-2016 The Regents of the University of California.
% All Rights Reserved.

if length(k) == 1
    k = [k,0];
end

N = size(flips,2);
A = [-R1(1)-k(1) +k(2)
    +k(1) -R1(2)-k(2)];

Mz0 = expm(A*Tin) * [1;0];

Mxy(1:2,1) = Mz0 .* sin(flips(:,1));
Mz(1:2,1) = Mz0 .* cos(flips(:,1));
for n = 2:N
    Mz_m = expm(A*TR) * Mz(:,n-1);
    Mxy(:,n) =  Mz_m .* sin(flips(:,n));
    Mz(:,n) = Mz_m .* cos(flips(:,n));
end