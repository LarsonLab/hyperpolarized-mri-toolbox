function [ AUCratio ] = compute_AUCratio( S )
%COMPUTE_AUCratio Calculates the Area-under-curve (AUC) ratio
%   A model-free approach for estimating metabolism.
%   Computes Area-under-curve for metabolites 
%   Hill et al. PLoS One, doi: 10.1371/journal.pone.0071996
% INPUTS
%	S - signal dynamics, must have substrate then product in 2nd to last dimension,
%   and time must be in the last dimension [voxels(optional), metabolites, time points]
% OUTPUTS
%   AUCratio - Area-under-curve (AUC) ratio
% EXAMPLES
%   Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);
%   % Mxy(1,:) is pyruvate signal, Mxy(2,:) is lactate signal 
%   AUCratio = compute_AUCratio(Mxy);

Ns = size(S);
if Ns > 2
    Stemp = reshape(S, [prod(Ns(1:end-2)), 2, Ns(end)]);
    AUCratio = sum(Stemp(:,2,:),3)/sum(Stemp(:,1,:),3);
    if Ns > 3
        AUCratio = reshape(AUCratio, Ns(1:end-2));
    end
else
    AUCratio = sum(S(2,:))/sum(S(1,:));
end


end

