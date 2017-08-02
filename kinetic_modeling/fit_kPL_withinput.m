function [params_fit, x1fit, x2fit, objective_val] = fit_kPL_withinput(S, TR, flips, params_fixed, params_est, noise_level, plot_flag)
% fit_kPL_withinput - kinetic model fitting of pyruvate to lactate conversion 
% This model requires estimation or fixing of a bolus arrival time ('Tarrival') and
% bolus end time ('Tend').
% It is described in Zierhut et al. J Magn Reson 202 (2010) 85–92.  
% doi:10.1016/j.jmr.2009.10.003
%
% [params_fit, x1fit, x2fit, objective_val] = fit_kPL(S, TR, flips, params_fixed, params_est, noise_level, plot_flag)
%
% All params_* values are structures, with possible fields of 'kPL', 'R1L',
% 'R1P', 'Rinj' (injection rate) - units of 1/s - and 'Tarrival', 'Tend'
% with units of s.
%
% INPUTS
%	S - signal dynamics [voxels, # of metabolites, # of time points]
%   TR - repetition time per time point flips - all flip angles [# of
%   metabolites, # of time points x # of phase encodes]
%	params_fixed - structure of fixed parameters and values.  parameters not in
%       this structure will be fit
%   params_est (optional) - structure of estimated values for fit parameters pyruvate to metabolites conversion rate initial guess (1/s)
%       Also can include upper and lower bounds on parameters as *_lb and
%       *_ub (e.g. R1L_lb, R1L_ub)
%   noise_level (optional) - estimate standard deviation of noise in data
%       to use maximum likelihood fit of magnitude data (with Rician noise
%       distribution)
%   plot_flag (optional) - plot fits
% OUTPUTS
%   params_fit - structure of fit parameters
%   x1fit, x2fit - fit curve for pyruvate and lactate, respectively
%   objective_val - measure of fit error
%
% EXAMPLES - see test_fit_kPL_fcn.m
%
% Authors: John Maidens,  Peder E. Z. Larson
%
% (c)2015-2017 The Regents of the University of California. All Rights
% Reserved.

params_all = {'kPL', 'R1L', 'R1P', 'Rinj', 'Tarrival', 'Tend'};
params_default_est = [0.02, 1/25, 1/25, 0.1, 0, 12];
params_default_lb = [0, 1/50, 1/50, 0, 0, 0];
params_default_ub = [Inf, 1/10, 1/10 Inf 30 30];

if nargin < 4 || isempty(params_fixed)
    params_fixed = struct([]);
end

if nargin < 5 || isempty(params_est)
    params_est = struct([]);
end

I_params_est = [];
for n = 1:length(params_all)
    if ~isfield(params_fixed, params_all(n))
        I_params_est = [I_params_est, n];
    end
end
Nparams_to_fit = length(I_params_est);

for n = 1:Nparams_to_fit
    param_name = params_all{I_params_est(n)};
    if isfield(params_est, param_name)
        params_est_vec(n) = params_est.(param_name);
    else
        params_est_vec(n) = params_default_est(I_params_est(n));
    end
    if isfield(params_est, [param_name '_lb'])
        params_lb(n) = params_est.([param_name '_lb']);
    else
        params_lb(n) = params_default_lb(I_params_est(n));
    end
    if isfield(params_est, [param_name '_ub'])
        params_ub(n) = params_est.([param_name '_ub']);
    else
        params_ub(n) = params_default_ub(I_params_est(n));
    end
end


if nargin < 6 || isempty(noise_level)
    % no noise level provided, so use least-squares fit (best for Gaussian
    % zero-mean noise)
    fit_method = 'ls';
else
    % otherwise use maximum likelihood estimator for Rician noise from
    % magnitudes
    fit_method = 'ml';
end

if nargin < 7
    plot_flag = 0;
end

if plot_flag
    disp('==== Computing parameter map ====')
end

size_S = size(S);  ndimsx = length(size_S)-2;
Nt = size_S(end); t = [0:Nt-1]*TR;
Nx = size_S(1:ndimsx);
if isempty(Nx)
    Nx = 1;
end
S = reshape(S, [prod(Nx), 2, Nt]);  % put all spatial locations in first dimension

[Sscale, Mzscale] = flips_scaling_factors(flips, Nt);

params_fit_vec = zeros([prod(Nx),Nparams_to_fit]);  objective_val = zeros([1,prod(Nx)]);
x1fit = zeros([prod(Nx),Nt]); x2fit = zeros([prod(Nx),Nt]);

for i=1:size(S, 1)
    if length(Nx) > 1 && plot_flag
        disp([num2str( floor(100*(i-1)/size(S, 1)) ) '% complete'])
    end
    % observed magnetization (Mxy)
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
        
        % Fit Pyruvate first (Tarrive, Rinj, Tend...), then full system ??
        
        % fit to data
        options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton');
        lsq_opts = optimset('Display','none','MaxIter', 500, 'MaxFunEvals', 500);
        switch(fit_method)
            case 'ls'
                obj = @(var) trajectory_difference(var, x1, x2, params_fixed, TR, Mzscale);
                [params_fit_vec(i,:),objective_val(i)] = lsqnonlin(obj, params_est_vec, params_lb, params_ub, lsq_opts);
                
            case 'ml'
                obj = @(var) negative_log_likelihood_rician(var, x1, x2, Mzscale, params_fixed, TR, noise_level.*(Sscale(2,:).^2));
                [params_fit_vec(i,:), objective_val(i)] = fminunc(obj, params_est_vec, options);
                
        end
        [x1fit(i,:), x2fit(i,:)] = trajectories_withinput(params_fit_vec(i,:), params_fixed, TR, Nt, Mzscale);
        x1fit(i,:) = x1fit(i,:)  .* Sscale(1, :);
        x2fit(i,:) = x2fit(i,:)  .* Sscale(2, :);
        
        if plot_flag
            % plot of fit for debugging
            figure(99)
            plot(t, x1, t, x2, t, x1fit(i,:)./ Sscale(1, :),'--', t, x2fit(i,:)./ Sscale(2, :), 'k:')
            xlabel('time (s)')
            ylabel('estimated state magnetization (au)')
            title(num2str(params_fit_vec(i,:)))
            legend('pyruvate', 'lactate', 'pyruvate fit', 'lactate fit')
            drawnow, pause(0.5)
        end
    end
end


params_fit = struct([]);
nfit = 0;
for n = 1:length(params_all)
    if ~isfield(params_fixed, params_all(n))
        nfit = nfit+1;
        params_fit(1).(params_all{n})= params_fit_vec(:,nfit);
    end
end

if length(Nx) > 1
    for n = 1:Nparams_to_fit
        params_fit.(params_est_fields{n}) = reshape(params_fit.(params_est_fields{n}), Nx);
    end
    
    
    x1fit = reshape(x1fit, [Nx, Nt]);
    x2fit = reshape(x2fit, [Nx, Nt]);
    objective_val = reshape(objective_val, Nx);
    disp('100 % complete')
end

end

function diff_all = trajectory_difference(params_fit, x1, x2,  params_fixed, TR, Mzscale)
[x1fit, x2fit] = trajectories_withinput(params_fit, params_fixed, TR, length(x1), Mzscale) ;
diff_all = [ x1(:)-x1fit(:) ; x2(:)-x2fit(:)];
end

function [ l1 ] = negative_log_likelihood_rician(params_fit, x1, x2, Mzscale, params_fixed, TR, noise_level)
%FUNCTION NEGATIVE_LOG_LIKELIHOOD_RICIAN Computes log likelihood for
%    compartmental model with Rician noise
% noise level is scaled for state magnetization (Mz) domain

N = size(x1,2);

% compute trajectory of the model with parameter values
[x1fit, x2fit] = trajectories_withinput(params_fit, params_fixed, TR, length(x1), Mzscale);

x = [x1; x2];
xfit = [x1fit; x2fit];
%compute negative log likelihood
l1 = 0;
for t = 1:N
    for k = 1:2
        l1 = l1 - (...
            log(x(k, t)) - log(noise_level(t)) ...
            - (x(k, t)^2 + xfit(k, t)^2)/(2*noise_level(t)) ...
            + x(k, t)*xfit(k, t)/noise_level(t) ...
            + log(besseli(0, x(k, t)*xfit(k, t)/noise_level(t), 1))...
            );
    end
end
end

