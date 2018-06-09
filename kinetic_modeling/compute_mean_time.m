function [ mean_time ] = compute_mean_time( S, dT )
%COMPUTE_MEAN_TIME Calculates the "mean time" for a signal
%   The mean time is defined as the center of mass of the time curve.
% INPUTS
%	S - signal dynamics, time must be in the last dimension [voxels and/or # of metabolites, # of time points]
%   dT - spacing of time points
% OUTPUTS
%   mean_time - mean time computed for all voxels and/or metabolites 

Ns = size(S);
tdim = length(Ns);
Nt = Ns(tdim);

tmat = repmat(reshape( (0:Nt-1)*dT, [ones(1,Ns-1),Nt]), size(S));
mean_time = sum( S .* tmat,tdim) ./ sum(pyr_time_input,tdim);

end

