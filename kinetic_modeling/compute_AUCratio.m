function [ AUCratio ] = compute_AUCratio( S )
%COMPUTE_AUCratio Calculates the Area-under-curve (AUC) ratio
%   A model-free approach for estimating metabolism.
%   Computes Area-under-curve for metabolites 
%   Hill et al. PLoS One, doi: 10.1371/journal.pone.0071996
% INPUTS
%	S - signal dynamics, must have substrate then product in 2nd to last dimension,
%   and time must be in the last dimension [metabolites, time points]
% OUTPUTS
%   AUCratio - Area-under-curve (AUC) ratio
% EXAMPLES
%   Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);
%   % Mxy(1,:) is pyruvate signal, Mxy(2,:) is lactate signal 
%   AUCratio = compute_AUCratio(Mxy);

AUCratio = sum(S(2,:))/sum(S(1,:));

end

