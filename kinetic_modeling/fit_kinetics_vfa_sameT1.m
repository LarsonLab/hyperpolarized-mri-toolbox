function [R1fit KPXfit  Sfit Minitfit] = fit_kinetics_vfa_sameT1(S, TR, flips, R1est, KPXest)
% fit_kinetics_vfa_sameT1 - Kinetic Model fitting routine that accounts for 
% arbitrary RF flip angles over time and between metabolites.
% This uses the following assumptions to provide more robust fit:
%   - No fitting of input function, fitting signal starts after input is
%   complete
%   - uni-directional conversion from substrate to 1 or 2 metabolic
%   products (i.e. pyruvate to lactate and pyruvate to alanine)
%   - Same T1 relaxation rate for all metabolites
%
% [R1fit KPXfit  Sfit Minitfit] = fit_kinetics_vfa_sameT1(S, TR, flips, R1est, KPXest)
%
% INPUTS
%	S - signal dynamics [# of metabolites, # of time points]
%   TR - repetition time per time point
%   flips - all flip angles [# of metabolites, # of time points x # of phase encodes]
%	R1est (optional) - relaxation rate initial guess (1/s)
%	kPXest (optional) - pyruvate to metabolites conversion rate initial guess (1/s)
% OUTPUTS
%   R1fit - fit relaxation rates (1/s)
%   KPXfit - fit conversion rates (1/s)
%   Sfit - fit signal dynamics
%   Minitfit - fit initial magnetizations
%
% Author: Peder E. Z. Larson
%
% (c)2013-2016 The Regents of the University of California.
% All Rights Reserved.

Nt = size(S,2);
Nmets = size(S,1);
if Nmets < 2 || Nmets > 3
    error('Only fitting of substrate to 1 or 2 metabolic products supported');
end

[Sscale Mzscale] = flips_scaling_factors(flips, Nt);

Mz_m0 = S(:,1) ./ Sscale(:,1);

if nargin < 5
    KPXest = 0.02*ones(1,Nmets-1);
end

if nargin < 4
    R1est = 1/25; % 1/s
end

% fit vector - initial guess
X0 = [R1est KPXest Mz_m0.']; 

    function Sest = model_exchange(x)
        
        if Nmets == 3
            A = [-x(1)-x(2)-x(3) 0 0;
                +x(2) -x(1) 0;
                +x(3) 0 -x(1)];
        else
            A = [-x(1)-x(2) 0;
                +x(2) -x(1)]; 
        end
        
        Mz_m(1:Nmets,1) = x(end-Nmets+1:end) ;
%        Mz_m(1:2,1) = Mz_m0;
        for It = 1:Nt*Nflips
             Mz_p(:,It) = Mz_m(:,It) .* cos(flips(:,It));
             Mxy(:,It) = Mz_m(:,It) .* sin(flips(:,It));
             Mz_m(:,It+1) = expm(A * TR/Nflips) * Mz_p(:,It);
        end
        Sest = squeeze(sum(reshape(Mxy,Nmets,Nflips,Nt),2));
%             Mz_p(:,It) = Mz_m(:,It) .* Mzscale(:,It);
%             Sest(:,It) = Mz_m(:,It) .* Sscale(:,It);
%             Mz_m(:,It+1) = expm(A * TR) * Mz_p(:,It);
%             
%         for It = 1:Nt
%         end
    end

    function res = g(x)
        res = model_exchange(x) - S;
        res = res(:);
    end

% some refinement possible here
opts = optimset('Display', 'off');%, 'FinDiffRelStep', [1 1 1 .1 .1]*sqrt(eps)); 
% FinDiffRelStep - use for normalization? or to preferentially fit certain parameters
% ie scale step size especially  for MzM0.  Or else scale data
% FinDiffType
%'MaxFunEvals', 1e8, 'TolFun', min(abs(S(:,end)))/1e14,'TolX', 1e-14);

%X = fminsearch(@g, X0, opts);
% set some bounds on the fitting
% 1/65 < R1 < 1/10, 0 < kPX < 1, 0 < Minit
lb = [1/65, zeros(1, Nmets-1), zeros(1,Nmets)];
ub = [1/10, ones(1, Nmets-1), Inf*ones(1,Nmets)];

Xfit = lsqnonlin(@g, X0, lb, ub, opts);
R1fit = Xfit(1); KPXfit = Xfit([2:Nmets]);
Minitfit = Xfit(end-Nmets+1:end);
Sfit = model_exchange(Xfit);

end