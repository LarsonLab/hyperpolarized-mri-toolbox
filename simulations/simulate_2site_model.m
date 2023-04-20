function [Mxy, Mz] = simulate_2site_model(Tin, R1, k, flips, TR, input_function)
% [Mxy, Mz] = simulate_2site_model(Tin, R1, k, flips, TR, [input_function])
%    OR 
% [Mxy, Mz] = simulate_2site_model(Mz0, R1, k, flips, TR, [input_function])
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

warning('simulate_2site_model() is now combined into simulate_Nsite_model() and maybe removed in a future toolbox release');

[Mxy, Mz] = simulate_Nsite_model(Tin, R1, k, flips, TR, input_function);

return

if isempty(R1)
  R1 = [0 0];  % infinite relaxation time
end

if length(R1) == 1
  R1 = [R1 R1];
end

N = size(flips,2);
A = [-R1(1)-k(1) +k(2)
    +k(1) -R1(2)-k(2)];
Ad = expm(A*TR);

if nargin < 6 || isempty(input_function) || all(input_function == 0)
    use_input_function = 0;
else
    use_input_function = 1;
    Bd = A\(Ad-eye(2));
end


if length(Tin) == 1  % remove Tin feature?
    Mz0 = expm(A*Tin) * [1;0];
else
    Mz0 = Tin(:);
end


Mxy(1:2,1) = Mz0 .* sin(flips(:,1));
Mz(1:2,1) = Mz0 .* cos(flips(:,1));
for n = 2:N
    if use_input_function
        Mz_m = Ad*(Mz(:,n-1)) + Bd(:,1).*input_function(n-1);
    else
        Mz_m = Ad*(Mz(:,n-1));
   
    end
    Mxy(:,n) =  Mz_m .* sin(flips(:,n));
    Mz(:,n) = Mz_m .* cos(flips(:,n));     
end