function [ TTP ] = compute_TTP( S, dT )
%COMPUTE_TTP Calculates the time-to-peak (TTP) for a signal
% INPUTS
%	S - signal dynamics, time must be in the last dimension [voxels and/or # of metabolites, # of time points]
%   dT - spacing of time points
% OUTPUTS
%   TTP - time-to-peak computed for all voxels and/or metabolites 

Ns = size(S);
tdim = length(Ns);
Nt = Ns(tdim);

[~, Imax] = max(S, [], tdim);

TTP = (Imax-1) * dT;  % basic

% quadratic peak fit
if Imax < Nt && Imax > 1
    Ipeak = 0.5 * (S(Imax-1) - S(Imax+1) ) ./ (S(Imax-1) - 2*S(Imax) + S(Imax+1) );
else
    Ipeak = 0;
end

TTP = (Imax+Ipeak-1) * dT;

end

