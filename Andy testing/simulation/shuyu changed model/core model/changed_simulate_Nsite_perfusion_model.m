function [Mxy, Mxy_iv, Mxy_ev, Mz, Mz_iv, Mz_ev] = changed_simulate_Nsite_perfusion_model(Tin, R1, k, flips, TR, VIF, VIFScale, kve)
% [Mxy, Mz] = simulate_Nsite_model(Mz0, R1, k, flips, TR, [input_function])
%
% Simulates the magnetization evolution in a N-site exchange model with
% single substrate and multiple products for varying flip angles.
% Allows for up to N=4
% Applied disretely as flips(:,1), TR, flips(:,2), TR, etc...
%
% INPUTS:
%   Mz0 [1xN] -  model starts from this  magnetization distribution
%       (If a single value is specified for Mz0, then this value is a time (s) for which the model will evolve prior to playing flips.  Magnetization is assumed to start all in substrate (e.g. pyruvate)
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


% nargin < 6
% isempty(VIF)
% all(VIF == 0)
% nargin < 6 || isempty(VIF)
% isempty(VIF) || all(VIF == 0)

if nargin < 6 || isempty(VIF) || all(VIF == 0)
    use_input_function = 0;
else
    use_input_function = 1;
    Nsim = 100;
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

Ad_TR = expm(A*TR);
if use_input_function
    Ad_Nsim = expm(A*TR/Nsim);
    Bd = A\(Ad_TR-eye(size(A,1)));
end

if length(Tin) == 1  % remove Tin feature?
    Mz0 = expm(A*Tin) * [1;0];
else
    Mz0 = Tin(:);
end


Mxy(1:Nmets,1) = Mz0 .* sin(flips(:,1));
Mz_ev(1:Nmets,1) = Mz0 .* cos(flips(:,1));%newVar: Mz_ev [nMet,nT]

M_iv = VIF.*VIFScale; %newVar: Mz_iv[1,nT], newInput: VIF[1,nT], VIFScale
Mz_iv = M_iv.*cos(flips(1,:));
Mxy_iv = M_iv.*sin(flips(1,:));

for n = 2:N
    if use_input_function
        % Mz_m = Mz(:,n-1);
        % % more accurate to spread out input over a number of samples to
        % % avoid unrealistically large signal jumps
        % for ni = 1:Nsim
        %     Mz_m =  Ad_Nsim * (Mz_m + [input_function(n-1)/Nsim;zeros(Nmets-1,1)]);
        % end
        Mz_m = Ad_TR*(Mz_ev(:,n-1)) + kve*Bd(:,1).*Mz_iv(n-1);%newInput:kve

    else
        Mz_m = Ad_TR * (Mz_ev(:,n-1));
    end
    Mxy(:,n) =  Mz_m .* sin(flips(:,n));
    Mz_ev(:,n) = Mz_m .* cos(flips(:,n));
end
Mz = Mz_ev;
Mz(1,:) = Mz(1,:) + Mz_iv;
Mxy_ev = Mxy;
Mxy(1,:) = Mxy(1,:) + Mxy_iv;