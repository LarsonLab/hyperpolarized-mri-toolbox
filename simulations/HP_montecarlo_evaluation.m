function [ hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, exp )
% [ hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, exp );
%
% Evaluate hyperpolarized carbon-13 MRI experiment using Monte Carlo
% simulations.
% Evaluation is performed considering a given set of acquisition
% parameters and kinetic modeling method.  The Area-under-curve ratio
% (AUCratio) derived in Hill et al. PLoS One, doi:
% 10.1371/journal.pone.0071996 , is always computed as a reference.
%
% INPUTS:
%   acq - structure containing acquisition parameters, including
%       TR, flips, N (number of timepoints)
%   fitting - structure containing fitting parameters, including
%       fit_fcn, params_est, params_fixed
%       (for use with fit_kPL* functions)    
%   exp - structure containing experimental parameters (optional, not yet implemen)
%
% OUTPUTS:
%   hdata, hsim - handles to figures from function
%
% (c) 2017 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson

% simulation and plotting parameters
NMC = 250;  % 100 maybe ok
Nexp_values = 8;
ratio_limits = [.8 1.2];

% fitting
fit_fcn = fitting.fit_fcn;
params_fixed = fitting.params_fixed;
params_est = fitting.params_est;


% experimental parameters
if nargin < 3 || isempty(exp)
    exp = struct([]);
end
params_all = {'kPL', 'R1L', 'R1P', 'std_noise', 'Tbolus'};
params_default = [0.02, 1/25, 1/25, 0.01, 12];

for n = 1:length(params_all)
    param_name = params_all{n};
    if ~isfield(params_fixed, param_name)
       exp.(param_name) = params_default(n);
    end
end

% default experiment values
R1 = [exp.R1P, exp.R1L]; kPL = exp.kPL; std_noise = exp.std_noise;
Tbolus = exp.Tbolus;

% experiment simulation ranges
exp.kPL_min = 0; exp.kPL_max = 0.03;    % approx kpL max in human studies
exp.std_noise_min = 0; exp.std_noise_max = 0.02;
exp.Tarrive_min = -5; exp.Tarrive_max = 5;
exp.Tbolus_min = 8; exp.Tbolus_max = 18;
exp.R1L_min = 1/15; exp.R1L_max = 1/35;
exp.R1P_min = 1/15; exp.R1P_max = 1/40;

% parameters to plot/fit: kPL, noise, arrival, duration, T1L, T1P
Nplot1 = 3; Nplot2 = 2;

t = [0:acq.N-1]*acq.TR;
Mz0 = [0,0];
input_function = gampdf(t,Tbolus/2,1);
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1

Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);
AUC_predicted = sum(Mxy(2,:))/sum(Mxy(1,:));

%% sample data
for Iflips = 1:size(acq.flips,3)
    [Mxy Mz] = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);
    Sn(1:size(Mxy,1), 1:size(Mxy,2),  Iflips) = Mxy + randn(size(Mxy))*std_noise;
end

hdata = figure;
subplot(121) , plot(t, squeeze(Sn(1,:,:))), title('sample simulated pyruvate')
xlabel('time (s)'), ylabel('Signal')
subplot(122) , plot(t, squeeze(Sn(2,:,:))), title('sample simulated lactate')
xlabel('time (s)'), ylabel('Signal')

%% setup for plots
hsim = figure;
Iplot = 1;

%% KPL test

kPL_test = linspace(exp.kPL_min, exp.kPL_max, Nexp_values).';

kPL_fit = zeros(length(kPL_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(kPL_test)
    Mxy = simulate_2site_model(Mz0, R1, [kPL_test(Itest) 0], acq.flips, acq.TR, input_function);
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
    Mxy = simulate_2site_model(Mz0, R1, [kPL_test(Itest) 0], acq.flips, acq.TR, input_function);
    AUC_predicted_test(Itest) = sum(Mxy(2,:))/sum(Mxy(1,:));
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_with_mean_and_std(kPL_test, kPL_fit./repmat(kPL_test(:),[1,NMC]),AUC_fit./repmat(AUC_predicted_test(:), [1, NMC]));
xlabel('k_{PL}'),  ylim(ratio_limits)

%% SNR

std_noise_test = linspace(exp.std_noise_min, exp.std_noise_max, Nexp_values).';

kPL_fit = zeros(length(std_noise_test), NMC); AUC_fit = kPL_fit;
Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);

for Itest = 1:length(std_noise_test)
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise_test(Itest), params_fixed, params_est);
    
end


subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_with_mean_and_std(std_noise_test, kPL_fit./kPL, AUC_fit./AUC_predicted);
xlabel('\sigma^2'),  ylim(ratio_limits)


%% bolus tests: arrival time

Tarrive_test= linspace(exp.Tarrive_min, exp.Tarrive_max, Nexp_values);

kPL_fit = zeros(length(Tarrive_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(Tarrive_test)
    t_test = t + Tarrive_test(Itest);
    input_function_test = interp1(t, input_function, t_test, 'linear', 0);
    Mz0_test = [sum(input_function)-sum(input_function_test) 0];
    
    Mxy = simulate_2site_model(Mz0_test, R1, [kPL 0], acq.flips, acq.TR, input_function_test);
    %     figure(99),  plot(t, Mxy), pause
    %
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_with_mean_and_std(Tarrive_test, kPL_fit./kPL, AUC_fit./AUC_predicted);
xlabel('Tarrive'), ylim(ratio_limits)


%% bolus tests: duration

Tbolus_test = linspace(exp.Tbolus_min, exp.Tbolus_max, Nexp_values);


kPL_fit = zeros(length(Tbolus_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(Tbolus_test)
    t_test = t*Tbolus/Tbolus_test(Itest);  % stretch/squeeze to modulate bolus
    input_function_test = interp1(t, input_function, t_test, 'linear', 0);
    input_function_test = input_function_test/sum(input_function_test) * sum(input_function); % normalized
    
    Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function_test);
    %     figure(99),  plot(t, Mxy), pause
    %
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_with_mean_and_std(Tbolus_test, kPL_fit./kPL, AUC_fit./AUC_predicted);
ylim(ratio_limits), xlabel('Tbolus')



%% T1 tests - remaining flaw...

R1L_test = linspace(exp.R1L_min, exp.R1L_max, Nexp_values);

kPL_fit = zeros(length(R1L_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(R1L_test)
    
    Mxy = simulate_2site_model(Mz0, [R1(1), R1L_test(Itest)], [kPL 0], acq.flips, acq.TR, input_function);
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_with_mean_and_std(R1L_test, kPL_fit/kPL, AUC_fit/AUC_predicted);
ylim(ratio_limits), xlabel('R_{1L}')

%% T1 tests

R1P_test = linspace(exp.R1P_min, exp.R1P_max, Nexp_values);

kPL_fit = zeros(length(R1P_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(R1P_test)
    
    Mxy = simulate_2site_model(Mz0, [R1P_test(Itest), R1(2)], [kPL 0], acq.flips, acq.TR, input_function);
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_with_mean_and_std(R1P_test, kPL_fit/kPL, AUC_fit/AUC_predicted);
ylim(ratio_limits), xlabel('R_{1P}')

end

function h=plot_with_mean_and_std(x, y1, y2)

h = gca;

DELTA1 = std(y1, [], 2);
DELTA2 = std(y2, [], 2);
Y1 = mean(y1,2);
Y2 = mean(y2,2);
plot(x, Y1, 'b-', x, Y1+DELTA1, 'b--', x, Y1-DELTA1, 'b--', ...
    x, Y2, 'g-', x, Y2+DELTA2, 'g--', x, Y2-DELTA2, 'g--')
end


function [kPL_fit, AUC_fit] = fitting_simulation(fit_fcn, Mxy, TR, flips, NMC, std_noise, params_fixed, params_est);

kPL_fit = zeros(1,NMC); AUC_fit = zeros(1,NMC);
parfor n = 1:NMC
    Sn = Mxy + randn(size(Mxy))*std_noise;
    
    params_fit = fit_fcn(Sn, TR, flips, params_fixed, params_est, [], 0);

    kPL_fit(n) = params_fit.kPL;
    
    AUC_fit(n) = sum(Sn(2,:))/sum(Sn(1,:));
end

end
