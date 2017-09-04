function [Sscale, Mzscale] = flips_scaling_factors(flips, Nt)
% flips_scaling_factors - compute magnetization scaling factors for each timepoint dataset from a
% series of RF pulses.  Useful for kinetic modeling, especially with complex flip angle strategies
% INPUTS
%   flips - [# of metabolites, # of excitations]
%       If excitation is different from Nt, then it is assumed that there are Nt/# of 
%	excitations that contribute to each resulting timepoint dataset
%   Nt - number of timepoints across which to calculate scaling
% OUTPUTS
%   Sscale - Signal scaling factors for correcting signals
%   Mzscale - Longitudinal magnetization scaling to account for in dynamics
%
% Author:  Peder E. Z. Larson
%
% (c)2015 The Regents of the University of California. All Rights
% Reserved.

% Compute the magnetization scaling 

% # of flip angles per timepoint (e.g. phase encodes)
Nflips = size(flips,2)/Nt;

Nmets = size(flips,1);

% magnetization and signal scaling factors for each timepoint based on flip
% angles
Mzscale = zeros(Nmets,Nt); Sscale = zeros(Nmets,Nt);
for t= 1:Nt
    Iflips = [1:Nflips] + (t-1)*Nflips;
    Mzscale(:,t) = prod(cos(flips(:,Iflips)),2);
    for n = 1:Nflips
        Sscale(:,t) = Sscale(:,t) + sin(flips(:,Iflips(n))) .* prod(cos(flips(:,Iflips(1:n-1))),2);
    end
end
Sscale = Sscale/Nflips;
