function [params_fit, Sfit, ufit, objective_val] = fit_kPL(S, TR, flips, params_fixed, params_est, noise_level, plot_flag)
% fit_kPL - Simple kinetic model fitting of conversion rate by fitting of
% product (e.g. lactate) signal at each time point.  Substrate (e.g.
% pyruvate) signal taken as is, and not fit to any function, eliminating
% need to make any assumptions about the input function.
% This uses the following assumptions:
%   - uni-directional conversion from substrate to metabolic products (i.e.
%   pyruvate to lactate)
%   - initial lactate magnetization is zero (need to add)
% It also allows for fixing of parameters. Based on simulations, our
% current recommendation is to fix pyruvate T1, as it doesn't impact kPL substantially.
%
% [params_fit, Sfit, ufit, objective_val] = fit_kPL(S, TR, flips, params_fixed, params_est, noise_level, plot_flag)
%
% All params_* values are structures, with possible fields of 'kPL', 'R1L',
% and 'R1P', and units of 1/s.
% INPUTS
%	S - signal dynamics [voxels, # of metabolites, # of time points]
%   TR - repetition time per time point flips - all flip angles [# of
%   metabolites, # of time points x # of phase encodes]
%	params_fixed - structure of fixed parameters and values (1/s).  parameters not in
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
%   Sfit - fit curve for lactate
%   ufit - derived input function (unitless)
%   objective_val - measure of fit error
%
% EXAMPLES - see test_fit_kPL_fcn.m
%
% Authors: John Maidens,  Peder E. Z. Larson
%
% (c)2015-2017 The Regents of the University of California. All Rights
% Reserved.

params_all = {'kPL', 'R1L', 'R1P'};
params_default_est = [0.02, 1/25, 1/25];
params_default_lb = [0, 1/60, 1/60];
params_default_ub = [Inf, 1/10, 1/10];

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
    % otherwise use maximum likelihood (good for Rician noise from
    % magnitudes)
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
Sfit = zeros([prod(Nx),Nt]); ufit = zeros([prod(Nx),Nt]);

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
        
        % fit to data
        options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton');
        lsq_opts = optimset('Display','none','MaxIter', 500, 'MaxFunEvals', 500);
        switch(fit_method)
            case 'ls'
                obj = @(var) (x2 - trajectories_frompyr(var, x1, Mzscale, params_fixed, TR)).*Sscale(2,:);  % perform least-squares in signal domain
                [params_fit_vec(i,:),objective_val(i)] = lsqnonlin(obj, params_est_vec, params_lb, params_ub, lsq_opts);
                
            case 'ml'
                obj = @(var) negative_log_likelihood_rician_frompyr(var, x1, x2, Mzscale, params_fixed, TR, noise_level.*(Sscale(2,:).^2));
                [params_fit_vec(i,:), objective_val(i)] = fminunc(obj, params_est_vec, options);
                
        end
        [Sfit(i,:), ufit(i,:)] = trajectories_frompyr(params_fit_vec(i,:), x1, Mzscale, params_fixed, TR);
        Sfit(i,:) = Sfit(i,:)  .* Sscale(2, :);
        ufit(i,:) = ufit(i,:)  .* Sscale(1, :);
        
        if plot_flag
            % plot of fit for debugging
            figure(99)
            subplot(2,1,1)
            plot(t, x1, t, x2, t, Sfit(i,:)./ Sscale(2, :),'--', t, ufit(i,:)./ Sscale(1, :), 'k:')
            xlabel('time (s)')
            ylabel('state magnetization (au)')
            subplot(2,1,2)
            plot(t, y1, t, y2, t, Sfit(i,:),'--', t, ufit(i,:), 'k:')
            xlabel('time (s)')
            ylabel('signal (au)')
            title(num2str(params_fit_vec(i,:),4))
            legend('pyruvate', 'lactate', 'lactate fit', 'input estimate')
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
        param_name = params_all{I_params_est(n)};
        params_fit.(param_name) = reshape(params_fit.(param_name), Nx);
    end
    
    
    Sfit = reshape(Sfit, [Nx, Nt]);
    ufit = reshape(ufit, [Nx, Nt]);
    objective_val = reshape(objective_val, Nx);
    disp('100 % complete')
end

end

function [ l1 ] = negative_log_likelihood_rician_frompyr(params_fit, x1, x2, Mzscale, params_fixed, TR, noise_level)
%FUNCTION NEGATIVE_LOG_LIKELIHOOD_RICIAN Computes log likelihood for
%    compartmental model with Rician noise
% noise level is scaled for state magnetization (Mz) domain

N = size(x1,2);

% compute trajectory of the model with parameter values
x2fit = trajectories_frompyr(params_fit, x1, Mzscale, params_fixed, TR);

% compute negative log likelihood
l1 = 0;
for t = 1:N
    for k = 1
        l1 = l1 - (...
            log(x2(k, t)) - log(noise_level(t)) ...
            - (x2(k, t)^2 + x2fit(k, t)^2)/(2*noise_level(t)) ...
            + x2(k, t)*x2fit(k, t)/noise_level(t) ...
            + log(besseli(0, x2(k, t)*x2fit(k, t)/noise_level(t), 1))...
            );
    end
end
end

