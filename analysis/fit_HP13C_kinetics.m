function [params_fit, Sfit, ufit, error_metrics] = fit_HP13C_kinetics(S, TR, flips, params_fixed, params_est, fitting_options)
% fit_HP13C_kinetics - Pharmacokinetic model fitting function for HP 13C MRI.
%
% Fits precursor-product kinetic model, assuming origination from a single
% substrate/precursor to multiple products
% It is nominally setup for modeling pyruvate to lactate, (optional) bicarbonate and alanine.
% An "input-less" method is used, eliminating need to make any assumptions about the input function.
%
% It also allows for fixing of parameters. Based on simulations, our
% current recommendation is to fix substrate relaxation rate, as it doesn't 
% impact kPX substantially.
%
% [params_fit, Sfit, ufit, error_metrics] = fit_HP13C_kinetics(S, TR, flips, params_fixed, params_est, fitting_options)
%
% All params_* values are structures, including possible fields of 'k', (1/s),
% 'R1P', 'R1L', 'R1B', 'R1A' (1/s).
% INPUTS
%   S - signal dynamics [voxels, # of metabolites, # of time points]
%		Substrate (e.g. Pyruvate) should be the first metabolite, followed by each product
%	TR (s) - [float] OR [array] - repetition time per time point
%	flips (radians) - all flip angles [# of metabolites, # of time points x # of phase encodes]
%	params_fixed (optional) - structure of fixed parameters and values (1/s).  parameters not in
%       this structure will be fit
%   params_est (optional) - structure of estimated values for fit parameters pyruvate to metabolites conversion rate initial guess (1/s)
%       Also can include upper and lower bounds on parameters as *_lb and
%       *_ub (e.g. R1L_lb, R1L_ub)
%   fitting_options (optional):
%       plot_flag (0/1) - plot fit results;
%       fit_method (string) - 'ls' for least-squares, 'ml' for maximum likelihood
%       noise_level (float) - noise level for maximum likelihood fitting
% OUTPUTS
%   params_fit - structure of fit parameters
%   Sfit - fit curves
%   ufit - derived input function (unitless)
%   error_metrics - measurements of fit error, including upper/lower
%       bounds, estimated kPL error, Rsquared, and Chisquared (untested)
%
% EXAMPLES - see test_fit_pyr_kinetics.m
%
% Authors: John Maidens,  Peder E. Z. Larson
%
% (c)2015-2025 The Regents of the University of California. All Rights
% Reserved.

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;


size_S = size(S);  ndimsx = length(size_S)-2;
Nt = size_S(end);
Nx = size_S(1:ndimsx);
Nmets = size_S(end-1);

if size(TR) == 1
    t = [0:Nt-1]*TR;
    TR = repmat(TR,[1,Nt]);
else
    t = cumsum(TR)-TR(1);
end

if isempty(Nx)
    Nx = 1;
end

% parse flip angles, assume phase_encodes=1
if length(flips) == 1 % only a single flip angle for all metabolites
    flips = repmat(flips, [Nmets Nt]);
elseif length(flips) == Nmets % only one flip angle for each metabolite, but same across time
    flips = repmat(flips(:), [1 Nt]);
end

[size_flips] = size(flips);
if all(size_flips == [Nt Nmets])
    flips = flips.';
elseif size_flips(1) ~= Nmets
    error('Flip angles should be of size [Nmets Nt*#phase_encodes], [Nmets 1], or [1]')
end

params_all = {'k',...
    'R1', ...
    'Mz0'};
params_size = [Nmets-1, Nmets, Nmets-1];
params_default_est = {ones(1,Nmets-1)*0.01, ...
    1/25*ones(1,Nmets), ...
    zeros(1,Nmets-1)};
params_default_lb = {ones(1,Nmets-1)*-Inf, ...
    1/50*ones(1,Nmets), ...
    -Inf*ones(1,Nmets-1)};
params_default_ub = {ones(1,Nmets-1)*Inf, ...
    1/10*ones(1,Nmets), ...
    Inf*ones(1,Nmets-1)};


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
params_est_vec = []; params_lb = []; params_ub = [];

for n = 1:Nparams_to_fit
    param_names{n} = params_all{I_params_est(n)};

    if isfield(params_est, param_names{n})
        params_est_vec = vertcat(params_est_vec, params_est.(param_names{n})(:));
    else
        params_est_vec = vertcat(params_est_vec, params_default_est{I_params_est(n)}(:));
    end
    if isfield(params_est, [param_names{n} '_lb'])
        params_lb = vertcat(params_lb, params_est.([param_names{n} '_lb'])(:));
    else
        params_lb = vertcat(params_lb, params_default_lb{I_params_est(n)}(:));
    end
    if isfield(params_est, [param_names{n} '_ub'])
        params_ub = vertcat(params_ub, params_est.([param_names{n} '_ub'])(:));
    else
        params_ub = vertcat(params_ub, params_default_ub{I_params_est(n)}(:));
    end
end

% detect whether user has provided fitting options, if not insert default values such as plot_flag
if nargin < 6 || isempty(fitting_options)
    fitting_options = struct([]);
end
if ~isfield(fitting_options, 'plot_flag')
    fitting_options.plot_flag = 0;
end
if ~isfield(fitting_options,'fit_method')
    fitting_options.fit_method = 'ls';
    noise_level = 0;
end

if fitting_options.plot_flag
    disp('==== Computing parameter map ====')
end

Sreshape = reshape(S, [prod(Nx), Nmets, Nt]);  % put all spatial locations in first dimension
% if Nmets < 4
%     Sreshape = cat(2, Sreshape, zeros([prod(Nx) 4-Nmets, Nt]));  % add zero data for unused metabolites
%     flips = cat(1, flips, ones([4-Nmets, size(flips,2)]));
% end

[Sscale, Mzscale] = flips_scaling_factors(flips, Nt);

params_fit_vec = zeros([prod(Nx),Nparams_to_fit]);  objective_val = zeros([1,prod(Nx)]);
lb = zeros([prod(Nx),Nparams_to_fit]); ub = zeros([prod(Nx),Nparams_to_fit]); err = zeros([prod(Nx),Nparams_to_fit]);
Sfit = zeros([prod(Nx),Nmets-1,Nt]); ufit = zeros([prod(Nx),Nt]);
Rsq = zeros([prod(Nx),Nmets-1]); CHIsq = zeros([prod(Nx),Nmets-1]);

for i=1:size(Sreshape, 1)
    % observed magnetization (Mxy)
    Mxy = reshape(Sreshape(i, :, :), [Nmets, Nt]);
    
    if any(Mxy(:) ~= 0)

        if prod(Nx) > 1
            disp([num2str( floor(100*(i-1)/size(Sreshape, 1)) ) '% complete'])
        end


        % estimate state magnetization (MZ) based on scaling from RF pulses
        Mz = Mxy./Sscale;
        
        % option to propogate inputless model from various points in time
        Istart = 1;
        
        % fit to data
        % options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton'); % not supported in Octave
        options = optimset('Display','none','Algorithm','quasi-newton');
        lsq_opts = optimset('Display','none','MaxIter', 500, 'MaxFunEvals', 500);
        
        switch(fitting_options.fit_method)
            case 'ls'
                obj = @(var) difference_inputless(var, params_fixed, TR, Mzscale, Sscale, Mz, Istart, Nmets) ;  % perform least-squares in signal domain
                [params_fit_vec(i,:),objective_val(i),resid,~,~,~,J] = lsqnonlin(obj, params_est_vec, params_lb, params_ub, lsq_opts);
    
                if ~isOctave
                    % extract 95% confidence interval on lactate timecourse fitting
                    CI = nlparci(params_fit_vec(i,:),resid,'jacobian',J);
                    lb(i,:) = CI(:,1);
                    ub(i,:) = CI(:,2);
                    err(i,:) = CI(:,2)-CI(:,1);
                end
            case 'ml'
                obj = @(var) negative_log_likelihood_rician_inputless(var, params_fixed, TR, Mzscale, Mz, noise_level.*(Sscale).^2, Istart, Nmets);
                [params_fit_vec(i,:), objective_val(i)] = fminunc(obj, params_est_vec, options);
        end
        
        [Mzfit, ufit(i,:)] = trajectories_inputless(params_fit_vec(i,:), params_fixed, TR,  Mzscale, Mz(1,:), Istart);
        
        Sfit(i,:,:) = Mzfit(2:Nmets,:)  .* Sscale(2:Nmets, :);
        ufit(i,:) = ufit(i,:)  .* Sscale(1, :) .* TR;
        
        I_flip = isfinite(Mz);
 
        % Note that Sfit only inlcudes products signals, not pyruvate,
        % hence the difference in indexing here and in plots below
        for I = 2:Nmets
            Rsq(i,I-1) = 1 - sum( (Sreshape(i,I,I_flip(I,:))-Sfit(i,I-1,I_flip(I,:))).^2 ,3 ) ...
                ./ sum( (Sreshape(i,I,I_flip(I,:)) - repmat(mean(Sreshape(i,I,I_flip(I,:)),3),[1 1 length(find(I_flip(I,:)))])).^2, 3);
            % chi squared distance - also pdist2 in Octave
            CHIsq(i,I-1) = sum( (Sreshape(i,I,I_flip(I,:))-Sfit(i,I-1,I_flip(I,:))).^2 ./ ...
                (Sreshape(i,I,I_flip(I,:))+Sfit(i,I-1,I_flip(I,:))) ,3) / 2;
        end
        
        if fitting_options.plot_flag
            % plot of fit for debugging
           
            figure(99)
            subplot(2,1,1)
            plot(t(I_flip(1,:)), Mz(1,I_flip(1,:)), '-x',t(I_flip(1,:)), ufit(i,I_flip(1,:))./ Sscale(1, I_flip(1,:)), 'k:')
            hold on
            for I = 2:Nmets
                plot(t(I_flip(I,:)), Mz(I,I_flip(I,:)), '-x', t, Mzfit(I,:), '--');
            end
            hold off
            xlabel('time (s)')
            ylabel('state magnetization (au)')

            subplot(2,1,2)
            plot(t(I_flip(1,:)), Mxy(1,I_flip(1,:)), '-x',t(I_flip(1,:)), ufit(i,I_flip(1,:)), 'k:')
            hold on
            for I = 2:Nmets
                plot(t(I_flip(I,:)), Mxy(I,I_flip(I,:)), '-x', t(I_flip(I,:)),  squeeze(Sfit(i,I-1,I_flip(I,:))), '--')
            end
            hold off
            xlabel('time (s)')
            ylabel('signal (au)')
            
            fit_results_string = [];
            for n = 1:Nparams_to_fit
                fit_results_string = [fit_results_string, param_names{n} ' = ' num2str(params_fit_vec(i,n),2) ' '];
            end            
            title(fit_results_string)
            disp(fit_results_string)
            
            products_legend{1} = 'substrate';
            products_legend{2} = 'input estimate';
            for n = 1:Nmets-1
                products_legend{2*n+1} = ['product ' num2str(n)];
                products_legend{2*n+2} = ['product ' num2str(n) ' fit'];
            end    
            legend( products_legend)
            drawnow, pause(0.5)
        end
    end
end

error_metrics=struct('objective_val', objective_val);
error_metrics.Rsq = Rsq;
error_metrics.CHIsq = CHIsq;

% below will have bugs
params_fit = struct([]);
for n = 1:Nparams_to_fit
    params_fit(1).(param_names{n})= params_fit_vec(:,n);
    if strcmp(fitting_options.fit_method, 'ls') && ~isOctave
        error_metrics.(param_names{n}).lb = lb(:,n);
        error_metrics.(param_names{n}).ub = ub(:,n);
        error_metrics.(param_names{n}).err = err(:,n);
    end
end

return

% fix below?

if length(Nx) > 1
    for n = 1:Nparams_to_fit
        params_fit.(param_names{n}) = reshape(params_fit.(param_names{n}), Nx);
        if strcmp(fitting_options.fit_method, 'ls') && ~isOctave
            error_metrics.(param_names{n}).lb = reshape(error_metrics.(param_names{n}).lb, Nx);
            error_metrics.(param_names{n}).ub = reshape(error_metrics.(param_names{n}).ub, Nx);
            error_metrics.(param_names{n}).err = reshape(error_metrics.(param_names{n}).err, Nx);
        end
    end
    
    Sfit = reshape(Sfit, [Nx, Nmets-1, Nt]);
    ufit = reshape(ufit, [Nx, Nt]);

    error_metrics.objective_val = reshape(error_metrics.objective_val, Nx);
    for n = 1:length(products_string)
        error_metrics.(products_string{n}).Rsq = reshape(error_metrics.(products_string{n}).Rsq, Nx);
        error_metrics.(products_string{n}).CHIsq =  reshape(error_metrics.(products_string{n}).CHIsq, Nx);
    end
    disp('100 % complete')
end

end

function convert_params_to_params_vec
end

function convert_params_vec_to_params
end

function diff_products = difference_inputless(params_fit, params_fixed, TR, Mzscale, Sscale, Mz, Istart, Nmets)
Mzfit = trajectories_inputless(params_fit, params_fixed, TR,  Mzscale, Mz(1,:), Istart) ;
temp_diff = (Mz - Mzfit) .* Sscale; % compute difference in the signal (Mxy) domain
diff_products = temp_diff(2:Nmets,:); % only differences are in the metabolic products
diff_products = diff_products(isfinite(diff_products));
end

function [ l1 ] = negative_log_likelihood_rician_inputless(params_fit, params_fixed, TR, Mzscale, Mz, noise_level, Istart, Nmets)
%FUNCTION NEGATIVE_LOG_LIKELIHOOD_RICIAN Computes log likelihood for
%    compartmental model with Rician noise
% noise level is scaled for state magnetization (Mz) domain

%https://doi.org/10.1109/42.712125

N = size(Mzscale,2);

% compute trajectory of the model with parameter values
Mzfit = trajectories_inputless(params_fit, params_fixed, TR,  Mzscale, Mz(1,:), Istart) ;

% compute negative log likelihood
l1 = 0;
for t = 1:N
    for k = 2:Nmets
        if isfinite(Mz(k,t))
            l1 = l1 - (...
                log(Mz(k, t)) - log(noise_level(k,t)) ...
                - (Mz(k, t)^2 + Mzfit(k, t)^2)/(2*noise_level(k,t)) ...
                + Mz(k, t)*Mzfit(k, t)/noise_level(k,t) ...
                + log(besseli(0, Mz(k, t)*Mzfit(k, t)/noise_level(k,t), 1))...
                );
        end
    end
end
end

function [Mz_all, u] = trajectories_inputless( params_fit, params_fixed, TR, Mzscale , Mz_substrate, Istart )
% Compute product magnetizations using a uni-directional two-site model
% Uses substrate magnetization measurements, estimated relaxation and
% conversion rates

if nargin < 6
    Istart = 1;
end

Nmets = size(Mzscale,1); N = size(Mzscale,2);
Mz_all = zeros(Nmets, N);
u = zeros(1,N);

params_all = {'k', 'R1', 'Mz0'};
params_size = [Nmets-1, Nmets, Nmets-1];

nfit = 0;
for n = 1:length(params_all)
    if isfield(params_fixed, params_all(n))
        eval([params_all{n} '= params_fixed.(params_all{n});']);
    else
        eval([params_all{n} '= params_fit(nfit + [1:params_size(n)]);']);
        nfit = nfit+params_size(n);
    end
end

switch Nmets
    case 2
        A = [-R1(1)-k(1) 0
            +k(1) -R1(2)];
    case 3
        A = [-R1(1)-k(1)-k(2) 0 0
            +k(1) -R1(2) 0
            +k(2) 0 -R1(3)];
    case 4
        A = [-R1(1)-k(1)-k(2)-k(3) 0 0 0
            +k(1) -R1(2) 0 0
            +k(2) 0 -R1(3) 0
            +k(3) 0 0 -R1(4)];
end


% account for timepoints where no pyruvate flip applied
I_substrate_flip = isfinite(Mz_substrate);
Mz_substrate = interp1(find(I_substrate_flip), Mz_substrate(I_substrate_flip), 1:length(Mz_substrate),'linear',0);

% inputless - force Mz substrate based on signal
Mz_all(1,:) = Mz_substrate;
Mz_all(2:end,Istart) = Mz0;

for It=Istart:N-1
    
    Mz_init = Mz_all(:,It) .* Mzscale(:, It);
    
    % estimate input, assuming this is constant during TR interval
    % This calculation could be improved for noise stability?
    u(It) = ( Mz_substrate(It+1) - Mz_init(1)*exp(A(1,1)*TR(It+1)) ) * -A(1,1) / (1 - exp(A(1,1))*TR(It+1));
    
    xstar = - inv(A)*[u(It),zeros(1,Nmets-1)].';
    
    % solve next time point under assumption of constant input during TR
    Mz_all(:,It+1) = xstar + expm(A*TR(It+1)) * (Mz_init - xstar);
    
    
end

% reverse in time
for It=Istart:-1:2
    
    Mz_init = Mz_all(:,It);
    
    % estimate input, assuming this is constant during TR interval
    % This calculation could be improved for noise stability?
    u(It-1) = ( Mz_substrate(It-1)*Mzscale(1,It-1) - Mz_init(1)*exp(A(1,1)*-TR(It-1)) ) * -A(1,1) / (1 - exp(A(1,1)*-TR(It-1)));
    
    xstar = - inv(A)*[u(It-1),zeros(1,Nmets-1)].';
    
    % solve previous time point under assumption of constant input during TR
    Mz_plus = xstar + expm(A*-TR(It-1)) * (Mz_init - xstar);
    
    
    Mz_all(:,It-1) = Mz_plus ./ Mzscale(:, It-1);
    
    
end


end
