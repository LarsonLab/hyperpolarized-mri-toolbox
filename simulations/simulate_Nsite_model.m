function [Mxy, Mz] = simulate_Nsite_model(Mz0, R1, k, flips, TR, input_function)
% [Mxy, Mz] = simulate_Nsite_model(Mz0, R1, k, flips, TR, [input_function])
%
% Simulates the magnetization evolution in a N-site exchange model with
% single substrate and multiple products for varying flip angles.
% Allows for up to N=4
% Applied disretely as flips(:,1), TR, flips(:,2), TR, etc...
%
% INPUTS:
%   Mz0 [1xN] -  model starts from this
%       magnetization distribution
%   R1 - relaxation times for N sites (1/s)
%   k [Nx2] - conversion rates between sites (1/s), first column is forward
%   rate, second column is reverse rate
%   flips - flip angles (radians)
%   TR - repetition time (s)
%   input_function (optional) - additional input function to add to
%   pyruvate for all time points.
%
% OUTPUTS:
%   Mxy,Mz - expected magnetization evolution immediately after flips
%
% (c) 2013-2014 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson

Nmets = size(flips,1);
N = size(flips,2);

if isempty(R1)
    R1 = zeros(1,Nmets);  % infinite relaxation time
end

if length(R1) == 1
    R1 = R1*ones(1,Nmets);
end

if nargin < 6 || isempty(input_function) % || all(input_function == 0)
    % optionally could extend to input for all metabolites
    input_function = zeros(1,N);
end

switch Nmets
    case 2
        A = [-R1(1)-k(1,1) +k(1,2)
            +k(1,1) -R1(2)-k(1,2)];
    case 3
        A = [-R1(1)-k(1,1)-k(2,1) +k(1,2) +k(2,2)
            +k(1,1) -R1(2)-k(1,2) 0
            +k(2,1) 0 -R1(3)-k(2,2)];
    case 4
        A = [-R1(1)-k(1,1)-k(2,1) +k(1,2) +k(2,2) +k(3,2)
            +k(1,1) -R1(2)-k(1,2) 0 0
            +k(2,1) 0 -R1(3)-k(2,2) 0
            +k(3,1) 0 0 -R1(4)-k(3,2)];
end

Mxy(1:Nmets,1) = Mz0 .* sin(flips(:,1));
Mz(1:Nmets,1) = Mz0 .* cos(flips(:,1));

for n = 2:N
    xstar = - inv(A)*[input_function(n-1),0,0,0].';
    
    % solve next time point under assumption of constant input during TR
    Mz_m = xstar + expm(A*TR) * (Mz(:,n-1) - xstar);

    Mxy(:,n) =  Mz_m .* sin(flips(:,n));
    Mz(:,n) = Mz_m .* cos(flips(:,n));
end