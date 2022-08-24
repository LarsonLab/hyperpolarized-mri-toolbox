function flips = mvfa_const_amp(N, TR, k12, R1, Minit);
% flips = mvfa_const_amp(N, TR, k12, R1, Minit)
%
% Calculates series of flip angles for substrate and product that
% are exchanging aiming for a constant amplitude signal in 
% hyperpolarized MR, assuming negligible thermal magnetization and back-reaction.
%
% INPUTS:
%   N - number of pulses
%   TR - repetition time, sec
%   k12 - conversion rate from substrate to product, 1/sec
%   R1 - relaxation rate of substrate and product (specify both)
%   Minit - relative initial magnetization levels
%
% OUTPUTS:
%   flips - flip angles
% 
% Example:  TR of 300 ms, for 30 s, with conversion rate of 0.05
%      1/s and T1s of 30 and 35 s, and initial magnetization level
%      of product at 10% of substrate
%   flips = mvfa_const_amp(100, 0.3, 0.05, [1/30 1/35], [1 0.1]);
%
% (c) 2012-2013 The Regents of the University of California
% All Rights Reserved.
%
% Authors: Peder E. Z. Larson, Yan Ann Xing 2012

if isempty(R1)
  R1 = [0 0];  % infinite relaxation time
end

if length(R1) == 1
  R1 = [R1 R1];
end


flips(1,:) = vfa_const_amp(N, pi/2, exp(-TR * (R1(1) + k12))); % substrate flip angle

Sbounds = [0 sum(Minit)];  % bound product signal

mag(:,1) = Minit;
A = [-R1(1)-k12 0;
    k12 -R1(2)];

while diff(Sbounds) > Minit(2)/1e12
    S_test = sum(Sbounds)/2;

    flips(2,1) = asin(S_test/mag(2,1));
    for k = 2:N
        mag(:,k)= expm(A*TR)*mag(:,k-1);
        flips(2,k)= asin(S_test/mag(2,k));
        mag(:,k) = mag(:,k) .* cos(flips(:,k)); % loss due to RF
    end
    
    if isreal(flips(2,N)) && flips(2,N) < pi/2
        Sbounds(1) = S_test;
    else
        Sbounds(2) = S_test;
    end
end

flips(2,N) = pi/2;
flips = real(flips);