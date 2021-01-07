function [params_fit, Sfit, ufit, err_metrics] = fit_kPL(S, TR, flips, params_fixed, params_est, noise_level, plot_flag)
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
% [params_fit, Sfit, ufit, err_metrics] = fit_kPL(S, TR, flips, params_fixed, params_est, noise_level, plot_flag)
%
% INPUTS
%	S - signal dynamics [voxels, # of metabolites, # of time points]
%   TR - repetition time per time point 
%   flips (radians) - all flip angles [# of metabolites, # of time points x # of phase encodes]
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
%   err_metrics - structure of fit error metrics (kPL error, R^2, Chi^2)
%
% EXAMPLES - see test_fit_kPL_fcn.m
%
% Authors: John Maidens, Peder E. Z. Larson, Daniele Mammoli, Natalie Korn
%
% (c)2015-2018 The Regents of the University of California. All Rights
% Reserved.

warning('fit_kPL() is now combined into fit_pyr_kinetics() and maybe removed in a future toolbox release');

if nargin < 4 || isempty(params_fixed)
    params_fixed = struct([]);
end

if nargin < 5 || isempty(params_est)
    params_est = struct([]);
end

if nargin < 6 
    noise_level = [];
end

if nargin < 7
    plot_flag = 0;
end



[params_fit, Sfit, ufit] = fit_pyr_kinetics(S, TR, flips, params_fixed, params_est, noise_level, plot_flag);


% err_metrics!! 
return


params_all = {'kPL', 'R1L', 'R1P', 'L0_start'};
err_all = {'ub', 'lb', 'err', 'Rsq', 'CHIsq'};
params_default_est = [0.01, 1/25, 1/25, 0];
params_default_lb = [-Inf, 1/60, 1/60, -Inf];
params_default_ub = [Inf, 1/10, 1/10, Inf];

if nargin < 4 || isempty(params_fixed)
    params_fixed = struct([]);
end

if nargin < 5 || isempty(params_est)
    params_est = struct([]);
end

if nargin < 6 
    noise_level = [];
end

if nargin < 7
    plot_flag = 0;
end

if isempty(noise_level)
    % no noise level provided, so use least-squares fit (best for Gaussian
    % zero-mean noise)
    fit_method = 'ls';
else
    % otherwise use maximum likelihood (good for Rician noise from
    % magnitudes)
    fit_method = 'ml';
    params_default_lb(4) = 0;  % set default lower bound for initial lactate to be non-negative
end

% Option to propogate fitting to model from different points in time
Istart = 1;


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
Ssize=size(S);
[Sscale, Mzscale] = flips_scaling_factors(flips, Nt);

params_fit_vec = zeros([prod(Nx),Nparams_to_fit]);  objective_val = zeros([1,prod(Nx)]);
err_vec=zeros([prod(Nx),5]);
Sfit = zeros([prod(Nx),Nt]); ufit = zeros([prod(Nx),Nt]);
fitsize=size(params_fit_vec);

for i=1:size(S, 1)
    if length(Nx) > 1 && plot_flag
        disp([num2str( floor(100*(i-1)/size(S, 1)) ) '% complete'])
    end
    % observed magnetization (Mxy)
    y1 = reshape(S(i, 1, :), [1, Nt]); % pyr
    y2 = reshape(S(i, 2, :), [1, Nt]); % lac
    
    if any(y1(:) ~= 0)
        % % plot of observed data for debugging
        % figure(1)
        % plot(t, y1, t, y2)
        % xlabel('time (s)')
        % ylabel('measured magnetization (au)')
        % legend('pyruvate', 'lactate')
        
        % estimate state magnetization (MZ) based on scaling from RF pulses
        x1 = y1./(Sscale(1, :)+eps);  % add eps to eliminate divide by zero errors
        x2 = y2./(Sscale(2, :)+eps);
        
        % normalize state magentization (MZ) so L0_start fitting parameter
        % are on a similar scale to kPL, R1P, etc otherwise fitting may
        % fail
        x_scale = max([x1(:);x2(:)]);
        x1_scaled = x1 / x_scale;
        x2_scaled = x2 / x_scale;
        
        % fit to data
        options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton');
        lsq_opts = optimset('Display','none','MaxIter', 500, 'MaxFunEvals', 500);
        switch(fit_method)
            case 'ls'
                obj = @(var) (x2_scaled - trajectories_frompyr(var, x1_scaled, Mzscale, params_fixed, TR, Istart)) .* Sscale(2,:);  % perform least-squares in signal domain
                [params_fit_vec(i,:),objective_val(i),resid,~,~,~,J] = lsqnonlin(obj, params_est_vec, params_lb, params_ub, lsq_opts);
                
                % extract 95% confidence interval on lactate timecourse fitting
                sigma = nlparci(params_fit_vec(i,1),resid,'jacobian',J(:,1));
                %sigma

            case 'ml'
                obj = @(var) negative_log_likelihood_rician_frompyr(var, x1_scaled, x2_scaled, Mzscale, params_fixed, TR, Istart, noise_level .* (Sscale(2,:).^2)/x_scale);
                [params_fit_vec(i,:), objective_val(i)] = fminunc(obj, params_est_vec, options);
                    
        end
        [Mzfit(i,:), ufit(i,:)] = trajectories_frompyr(params_fit_vec(i,:), x1_scaled, Mzscale, params_fixed, TR, Istart);
        Sfit(i,:) = Mzfit(i,:)*x_scale  .* Sscale(2, :);
        ufit(i,:) = ufit(i,:)*x_scale  .* Sscale(1, :);
        
        % export goodness of fit parameters (ub, lb, total error, R^2, chi^2)
        if exist('sigma', 'var')
            % export goodness of fit parameters (ub, lb, total error, R^2, chi^2)
            err_vec(i,1)=sigma(1,2); % lb
            err_vec(i,2)=sigma(1,1); % ub
            err_vec(i,3)=sigma(1,2)-sigma(1,1); % err
        end

        err_vec(i,4)=1-sum( (y2-Sfit(i,:)).^2 ) ./ sum( (y2 - mean(y2)).^2 );  % Rsquared
        err_vec(i,5)= sum( (y2-Sfit(i,:)).^2 ./ (y2 + Sfit(i,:)) ) /2; % CHI squared distance
        
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
for n = 1:length(params_all)-1  % don't output L0_start
    if ~isfield(params_fixed, params_all(n))
        nfit = nfit+1;
        params_fit(1).(params_all{n})= params_fit_vec(:,nfit);
    end
end

if length(Nx) > 1
    for n = 1:Nparams_to_fit-1 % don't output L0_start
        param_name = params_all{I_params_est(n)};
        params_fit.(param_name) = reshape(params_fit.(param_name), Nx);
    end
    
    Sfit = reshape(Sfit, [Nx, Nt]);
    ufit = reshape(ufit, [Nx, Nt]);
    objective_val = reshape(objective_val, Nx);
	if plot_flag
    	disp('100 % complete')
    end
end

err_metrics=struct([]);
nfit=0;
for n = 1:length(err_all)
    nfit=nfit+1;
    err_metrics(1).(err_all{n})= err_vec(:,nfit);
end

if length(Nx) > 1
    for n = 1:length(err_all)
        param_name = err_all{(n)};
        err_metrics.(param_name) = reshape(err_metrics.(param_name), Nx);
    end
end

end

function [ l1 ] = negative_log_likelihood_rician_frompyr(params_fit, x1, x2, Mzscale, params_fixed, TR, Istart, noise_level)
%FUNCTION NEGATIVE_LOG_LIKELIHOOD_RICIAN Computes log likelihood for
%    compartmental model with Rician noise
% noise level is scaled for state magnetization (Mz) domain

N = size(x1,2);

% compute trajectory of the model with parameter values
x2fit = trajectories_frompyr(params_fit, x1, Mzscale, params_fixed, TR, Istart);

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

