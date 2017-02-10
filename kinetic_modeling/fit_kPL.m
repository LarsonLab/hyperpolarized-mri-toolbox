function [kPLfit, Sfit, objective_val] = fit_kPL(S, TR, flips, R1_fixed, kPL_est, noise_level, plot_flag)
% fit_kPL - Simple kinetic model fitting of conversion rate by fitting of
% product (e.g. lactate) signal at each time point.  Substrate (e.g.
% pyruvate) signal taken as is, and not fit to any function, eliminating
% need to fit input function. This uses the following assumptions:
%   - uni-directional conversion from substrate to metabolic products (i.e.
%   pyruvate to lactate)
%   - initial lactate magnetization is zero
%   - Fixed T1 relaxation rates for  metabolites
%
% [kPL_hat, objective_val] = fit_kPL(S, TR, flips, R1_fixed, kPL_est, noise_level, plot_flag)
%
% INPUTS
%	S - signal dynamics [voxels, # of metabolites, # of time points]
%   TR - repetition time per time point flips - all flip angles [# of
%   metabolites, # of time points x # of phase encodes]
%	R1_fixed (optional) - fixed relaxation rates (1/s) 
%   kPL_est (optional) - pyruvate to metabolites conversion rate initial guess (1/s)
%   noise_level (optional) - estimate standard deviation of noise at each time point to use
%   maximum likelihood fit of magnitude data (with Rician noise
%   distribution)
%   plot_flag (optional) - plot fits
% OUTPUTS
%   KPLfit - fit conversion rate (1/s)
%   Sfit - fit conversion rate (1/s)
%   objective_val - measure of fit error
%
% Authors: John Maidens,  Peder E. Z. Larson
%
% (c)2015-2016 The Regents of the University of California. All Rights
% Reserved.

if nargin < 4 || isempty(R1_fixed)
    R1_fixed = [1/25 1/25];
end

R1P_fixed = R1_fixed(1);
if length(R1_fixed) == 1
    R1L_fixed = R1_fixed(1);
else
    R1L_fixed = R1_fixed(2);
end


if nargin < 5 || isempty(kPL_est)
    kPL_est = 0.01;
end

if nargin < 6 || isempty(noise_level)
    % no noise level provided, so use least-squares fit (best for Gaussian
    % zero-mean noise)
    fit_method = 'ls';
else
    % otherwise use maximum likelihood (good for Rician noise from
    % magnitudes)
    fit_method = 'ml';
end

if nargin < 7
    plot_flag = 0;
end


disp('==== Computing parameter map ====')

size_S = size(S);  ndimsx = length(size_S)-2;
Nt = size_S(end); t = [0:Nt-1]*TR;
Nx = size_S(1:ndimsx);
if isempty(Nx)
    Nx = 1;
end
S = reshape(S, [prod(Nx), 2, Nt]);  % put all spatial locations in first dimension

[Sscale, Mzscale] = flips_scaling_factors(flips, Nt);

kPLfit = zeros([1,prod(Nx)]);  objective_val = kPLfit;
Sfit = zeros([prod(Nx),Nt]); ufit = zeros([prod(Nx),Nt]);

for i=1:size(S, 1)
if length(Nx) > 1
    disp([num2str( floor(100*(i-1)/size(S, 1)) ) '% complete'])
end
    % observed magnetization
    y1 = reshape(S(i, 1, :), [1, Nt]); % pyr
    y2 = reshape(S(i, 2, :), [1, Nt]); % lac
    if any(y1 ~= 0)
        % % plot of observed data for debugging
        % figure(1)
        % plot(t, y1, t, y2)
        % xlabel('time (s)')
        % ylabel('measured magnetization (au)')
        % legend('pyruvate', 'lactate')
        
        % estimate state magnetization (MZ) based on scaling from RF pulses
        x1 = y1./Sscale(1, :);
        x2 = y2./Sscale(2, :);
        
        % fit best kPL to data
        options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton');
        switch(fit_method)
            case 'ls'
                obj = @(kPL) norm(x2 - trajectories(kPL, x1, Mzscale, R1P_fixed, R1L_fixed, TR));
                
            case 'ml'
                obj = @(kPL) negative_log_likelihood_rician(kPL, x1, x2, Mzscale, R1P_fixed, R1L_fixed, TR, noise_level./Sscale(2,:));
                
        end
        [kPLfit(i), objective_val(i)] = fminunc(obj, kPL_est, options);
        [Sfit(i,:), ufit(i,:)] = trajectories(kPLfit(i), x1, Mzscale, R1P_fixed, R1L_fixed, TR);
        Sfit(i,:) = Sfit(i,:)  .* Sscale(2, :);
        
        if plot_flag
            % plot of fit for debugging
            figure(99)
            plot(t, x1/10, t, x2, t, Sfit(i,:)./ Sscale(2, :),'--', t, ufit(i,:)./ Sscale(1, :), 'k:')
            xlabel('time (s)')
            ylabel('estimated state magnetization (au)')
            title(num2str(kPLfit(i)))
            legend('pyruvate', 'lactate', 'lactate fit', 'input estimate')
            drawnow, pause(0.2)
        end
    end
end

if length(Nx) > 1
    kPLfit = reshape(kPLfit, Nx);
    Sfit = reshape(Sfit, [Nx, Nt]);
    objective_val = reshape(objective_val, Nx);
    disp('100 % complete')
end

end



