function [results, hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, exp )
% [ results, hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, exp );
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
%   exp - structure containing experimental parameters (optional)
%
% OUTPUTS:
%   results - structure containing summary of results
%   hdata, hsim - handles to figures from function
%
% (c) 2017 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson

% simulation and plotting parameters
NMC = 250;  % 100 maybe ok
Nexp_values = 8;
ratio_limits = [-.2 .2];

% fitting
fit_fcn = fitting.fit_fcn;
params_fixed = fitting.params_fixed;
params_est = fitting.params_est;
if isfield(fitting, 'NMC')
    NMC = fitting.NMC;
end


% experimental parameters
if nargin < 3 || isempty(exp)
    exp = struct([]);
end
exp_params_all = {'kPL', 'R1L', 'R1P', 'std_noise', 'Tbolus', 'Tarrival'};
exp_params_default = [0.02, 1/25, 1/30, 0.005, 12, 4];

for n = 1:length(exp_params_all)
    param_name = exp_params_all{n};
    if ~isfield(exp, param_name)
       exp.(param_name) = exp_params_default(n);
    end
end

% default experiment values
R1 = [exp.R1P, exp.R1L]; kPL = exp.kPL; std_noise = exp.std_noise;
Tbolus = exp.Tbolus;

% experiment simulation ranges
exp.kPL_min = 0; exp.kPL_max = 0.04;    % approx kpL max in human studies
exp.std_noise_min = 0; exp.std_noise_max = 0.01;
exp.Tarrival_min = 0; exp.Tarrival_max = 8;
exp.Tbolus_min = 10; exp.Tbolus_max = 14;
exp.R1L_min = 1/35; exp.R1L_max = 1/15;
exp.R1P_min = 1/40; exp.R1P_max = 1/20;
exp.B1error_min = -.2; exp.B1error_max = .2;
exp.B1diff_min = -.2; exp.B1diff_max = .2;

% parameters to plot/fit: kPL, noise, arrival, duration, T1L, T1P, B1,
% vascular parameters?
Nplot1 = 4; Nplot2 = 2;

t = [0:acq.N-1]*acq.TR;
Mz0 = [0,0];
t_input = t+acq.TR-exp.Tarrival;
if isfield(exp, 'input_function')
    input_function = exp.input_function;
else
    input_function = gampdf(t_input,4,Tbolus/4);  % gives a full-width half-max of the bolus of ~ Tbolus sec
end
input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
results.input_function = input_function;

Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);
AUC_predicted = sum(Mxy(2,:))/sum(Mxy(1,:));

%% sample data

[Mxy Mz] = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);
results.sample_data = Mxy + randn(size(Mxy))*std_noise;
results.sample_data_time = t;

hdata = figure;
subplot(121) , plot(t, squeeze(results.sample_data(1,:))), title('sample simulated pyruvate')
xlabel('time (s)'), ylabel('Signal')
subplot(122) , plot(t, squeeze(results.sample_data(2,:))), title('sample simulated lactate')
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
%[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(kPL_test, kPL_fit./repmat(kPL_test(:),[1,NMC]),AUC_fit./repmat(AUC_predicted_test(:), [1, NMC]));
%ylim(ratio_limits)
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(kPL_test, (kPL_fit - repmat(kPL_test(:),[1,NMC]))/kPL,(AUC_fit - repmat(AUC_predicted_test(:), [1, NMC]))/AUC_predicted);
ylim(ratio_limits)
xlabel('k_{PL}'),  xlim([exp.kPL_min, exp.kPL_max])

results.kPL_test.kPL_avg_error = mean(kPL_std) ;  % precision measurement - normalized for comparison with other parameters
results.kPL_test.kPL_avg_bias = mean(abs(kPL_mean)) ;  % accuracy measurement
results.kPL_test.kPL_std_bias = std(kPL_mean) ;  % accuracy measurement

results.kPL_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.kPL_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.kPL_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement

%% SNR

std_noise_test = linspace(exp.std_noise_min, exp.std_noise_max, Nexp_values).';

kPL_fit = zeros(length(std_noise_test), NMC); AUC_fit = kPL_fit;
Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);

for Itest = 1:length(std_noise_test)
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise_test(Itest), params_fixed, params_est);
    
end


subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(std_noise_test, kPL_fit./kPL-1, AUC_fit./AUC_predicted-1);
xlabel('\sigma'),  xlim([exp.std_noise_min, exp.std_noise_max])
ylim(ratio_limits)


results.noise_test.kPL_avg_error = mean(kPL_std);  % precision measurement
results.noise_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
results.noise_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement

results.noise_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.noise_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.noise_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement


%% bolus tests: arrival time

Tarrival_test= linspace(exp.Tarrival_min, exp.Tarrival_max, Nexp_values);

kPL_fit = zeros(length(Tarrival_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(Tarrival_test)
    t_test = t_input - (Tarrival_test(Itest)-exp.Tarrival);
    input_function_test = interp1(t_input, input_function, t_test, 'linear', 0);
    Mz0_test = [sum(input_function)-sum(input_function_test) 0];
    
    Mxy = simulate_2site_model(Mz0_test, R1, [kPL 0], acq.flips, acq.TR, input_function_test);
    %     figure(99),  plot(t, Mxy), pause
    %
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(Tarrival_test, kPL_fit./kPL-1, AUC_fit./AUC_predicted-1);
xlabel('Tarrival'), xlim([exp.Tarrival_min, exp.Tarrival_max]), ylim(ratio_limits)

results.Tarrival_test.kPL_avg_error = mean(kPL_std);  % precision measurement
results.Tarrival_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
results.Tarrival_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement

results.Tarrival_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.Tarrival_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.Tarrival_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement



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
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(Tbolus_test, kPL_fit./kPL-1, AUC_fit./AUC_predicted-1);
ylim(ratio_limits), xlim([exp.Tbolus_min, exp.Tbolus_max]), xlabel('Tbolus')

results.Tbolus_test.kPL_avg_error = mean(kPL_std);  % precision measurement
results.Tbolus_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
results.Tbolus_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement

results.Tbolus_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.Tbolus_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.Tbolus_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement




%% T1 lactate tests - this is hard

R1L_test = linspace(exp.R1L_min, exp.R1L_max, Nexp_values);

kPL_fit = zeros(length(R1L_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(R1L_test)
    
    Mxy = simulate_2site_model(Mz0, [R1(1), R1L_test(Itest)], [kPL 0], acq.flips, acq.TR, input_function);
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(1./R1L_test, kPL_fit/kPL-1, AUC_fit/AUC_predicted-1);
ylim(ratio_limits), xlim(1./[exp.R1L_max, exp.R1L_min]), xlabel('T_{1L}')

results.R1L_test.kPL_avg_error = mean(kPL_std);  % precision measurement
results.R1L_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
results.R1L_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement

results.R1L_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.R1L_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.R1L_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement


%% T1 pyruvate tests

R1P_test = linspace(exp.R1P_min, exp.R1P_max, Nexp_values);

kPL_fit = zeros(length(R1P_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(R1P_test)
    
    Mxy = simulate_2site_model(Mz0, [R1P_test(Itest), R1(2)], [kPL 0], acq.flips, acq.TR, input_function);
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(1./R1P_test, kPL_fit/kPL-1, AUC_fit/AUC_predicted-1);
ylim(ratio_limits), xlim(1./[exp.R1P_max, exp.R1P_min]), xlabel('T_{1P}')

results.R1P_test.kPL_avg_error = mean(kPL_std);  % precision measurement
results.R1P_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
results.R1P_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement

results.R1P_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.R1P_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.R1P_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement


%% B1 error tests
    % simulate inaccurate B1, & unknown

B1error_test = linspace(exp.B1error_min, exp.B1error_max, Nexp_values);

kPL_fit = zeros(length(B1error_test), NMC); AUC_fit = kPL_fit;

for Itest = 1:length(B1error_test)
    Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips * (1+B1error_test(Itest)), acq.TR, input_function);
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips, NMC, std_noise, params_fixed, params_est);
    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(B1error_test, kPL_fit/kPL-1, AUC_fit/AUC_predicted-1);
ylim(ratio_limits), xlim([exp.B1error_min, exp.B1error_max]), xlabel('% B_{1} error')

results.B1error_test.kPL_avg_error = mean(kPL_std);  % precision measurement
results.B1error_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
results.B1error_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement

results.B1error_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.B1error_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.B1error_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement


%% B1 difference tests
    % simulate inaccurate B1, but known

B1diff_test = linspace(exp.B1diff_min, exp.B1diff_max, Nexp_values);

kPL_fit = zeros(length(B1diff_test), NMC); AUC_fit = kPL_fit;
clear AUC_predicted_test

for Itest = 1:length(B1diff_test)
    Mxy = simulate_2site_model(Mz0, R1, [kPL 0], acq.flips * (1+B1diff_test(Itest)), acq.TR, input_function);
    
    [kPL_fit(Itest,:), AUC_fit(Itest,:)] = fitting_simulation(fit_fcn,Mxy, acq.TR, acq.flips* (1+B1diff_test(Itest)), NMC, std_noise, params_fixed, params_est);
    AUC_predicted_test(Itest) = sum(Mxy(2,:))/sum(Mxy(1,:));
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
[~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(B1diff_test, kPL_fit/kPL-1, AUC_fit./repmat(AUC_predicted_test(:), [1, NMC]) -1);
ylim(ratio_limits), xlim([exp.B1diff_min, exp.B1diff_max]), xlabel('% B_{1} difference')

results.B1diff_test.kPL_avg_error = mean(kPL_std);  % precision measurement
results.B1diff_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
results.B1diff_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement

results.B1diff_test.AUC_avg_error = mean(AUC_std);  % precision measurement
results.B1diff_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
results.B1diff_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement


end

function [h,Y1,Y2,DELTA1,DELTA2]=plot_with_mean_and_std(x, y1, y2)

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
