function [results, hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, experiment )
% [ results, hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, experiment );
%
% Evaluate hyperpolarized carbon-13 MRI experiment using Monte Carlo
% simulations.
% Evaluation is performed considering a given set of acquisition
% parameters and kinetic modeling method.  
% A calibrated Area-under-curve ratio
% (cAUCratio) derived in Hill et al. PLoS One, doi:
% 10.1371/journal.pone.0071996 , is always computed as a reference.
%
% INPUTS:
%   acq - structure containing acquisition parameters, must include
%       TR, flips, N (number of timepoints)
%   fitting - structure containing fitting/modeling function and relevant
%       parameters, including
%       fit_fcn, params_est, params_fixed
%       (for use with fit_pyr_kinetics* functions)    
%       To compare multiple
%   experiment - structure containing experimental parameters and ranges for simulations (optional)
%       USE THIS TO SPECIFY WHICH SIMULATIONS TO RUN?
%
% OUTPUTS:
%   results - structure containing summary of results
%   hdata, hsim - handles to figures from function
%
% (c) 2017 The Regents of the University of California
% All Rights Reserved.
%
% Author: Peder E. Z. Larson

global isOctave;
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;

% simulation and plotting parameters
Nexp_values = 8;
ratio_limits = [-.2 .2];

% fitting
N_fitting_methods = length(fitting);

% experimental parameters
if nargin < 3 || isempty(experiment)
    experiment = struct([]);
end

% default experiment values - chosen based on UCSF Human Prostate Cancer
% studies
exp_params_all = {'kPL', 'R1L', 'R1P', 'std_noise', 'Tbolus', 'Tarrival', ... % nominal values
    'kPL_min', 'kPL_max', 'R1L_min', 'R1L_max', 'R1P_min', 'R1P_max', ...
    'std_noise_min', 'std_noise_max', 'Tbolus_min', 'Tbolus_max', 'Tarrival_min', 'Tarrival_max', ...
    'B1error_min', 'B1error_max', 'B1diff_min', 'B1diff_max', 'NMC'}; % experiment simulation ranges
exp_params_default = [0.02, 1/25, 1/30, 0.005, 8, 4 ... % nominal values
    0.000, 0.04, 1/35, 1/15, 1/40, 1/20, ...
    0, 0.01, 6, 10, 0, 8, ...
    -.2, .2, -.2, .2, 250]; % experiment simulation ranges

for n = 1:length(exp_params_all)
    param_name = exp_params_all{n};
    if ~isfield(experiment, param_name)
       experiment.(param_name) = exp_params_default(n);
    end
end

R1 = [experiment.R1P, experiment.R1L]; kPL = experiment.kPL;
Tbolus = experiment.Tbolus;

Nplot1 = 4; Nplot2 = 2;

% default input function and sample data
t = [0:acq.N-1]*acq.TR;
Mz0 = [0,0];
if isfield(experiment, 'input_function')
    input_function = experiment.input_function;
    input_function = input_function/sum(input_function); % normalize so total input magnetization = 1
    t_input = t+acq.TR-experiment.Tarrival;
else
    [input_function, t_input] = realistic_input_function(acq.N, acq.TR, experiment.Tarrival, experiment.Tbolus);  % gives a full-width half-max of the bolus of ~ Tbolus sec
end
results.input_function = input_function;

%% sample data at nominal parameter values

[Mxy Mz] = simulate_Nsite_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);
results.sample_data = Mxy + randn(size(Mxy))*experiment.std_noise;
results.sample_data_time = t;

hdata = figure;
subplot(121) , plot(t, squeeze(results.sample_data(1,:))/experiment.std_noise), title('sample simulated pyruvate')
xlabel('time (s)'), ylabel('Signal (SNR)')
subplot(122) , plot(t, squeeze(results.sample_data(2,:))/experiment.std_noise), title('sample simulated lactate')
xlabel('time (s)'), ylabel('Signal (SNR)')

%% setup for plots
for f = 1:N_fitting_methods
    legend_description{f} = fitting(f).fit_description;
    fitting(f).nominal_metric_value = compute_nominal_metric_value(fitting(f).fit_fcn, Mxy, kPL);
    nominal_metric(f) = compute_nominal_metric_value(fitting(f).fit_fcn, Mxy, kPL);
    
    switch func2str(fitting(f).fit_fcn)
        case 'compute_AUCratio'
            fitting(f).input_arguments = {};
        case {'compute_mean_time', 'compute_TTP'} % will use lactate mean time or TTP
            fitting(f).input_arguments = {acq.TR};
        otherwise
            fitting(f).input_arguments = {acq.TR, acq.flips, fitting(f).params_fixed, fitting(f).params_est, [], 0};
    end
end

nominal_metric_matrix=  repmat( nominal_metric(:), [1, experiment.NMC Nexp_values]);

warnStruct = warning('off','all');

%% KPL test
clear metric_mean metric_std metric_fits

kPL_test = linspace(experiment.kPL_min, experiment.kPL_max, Nexp_values).';

for Itest = 1:length(kPL_test)

    Mxy = simulate_Nsite_model(Mz0, R1, [kPL_test(Itest) 0], acq.flips, acq.TR, input_function);
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips, experiment.NMC, experiment.std_noise);

    for f = 1:N_fitting_methods
        nominal_metric_kPL_test(f,Itest) = compute_nominal_metric_value(fitting(f).fit_fcn, Mxy, kPL_test(Itest));
    end
    
end

hsim{1} = figure;
plot_fit_results(fitting, kPL_test, metric_fits, 'k_{PL} (1/s)');

hsim{2} = figure;
Iplot = 1;

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, kPL_test(2:end), ...
    metric_fits(:,:,2:end)./repmat( reshape(nominal_metric_kPL_test(:,2:end), [N_fitting_methods 1 Nexp_values-1]), [1, experiment.NMC 1])-1, 'k_{PL} (1/s)')
xlim([kPL_test(2), experiment.kPL_max])

% add legend

% In Octave I wasn't able to handle the legend with so many subplots 
if ~isOctave
    legh = legend;
    legh.Position = [0.35 0.01 0.3 0.1];
end

% 
results.kPL_test.kPL_test = kPL_test;
results.kPL_test.metric_mean = metric_mean;
results.kPL_test.metric_std = metric_std;

%% SNR
clear metric_mean metric_std metric_fits

std_noise_test = linspace(experiment.std_noise_min, experiment.std_noise_max, Nexp_values).';

Mxy = simulate_Nsite_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function);

for Itest = 1:length(std_noise_test)
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips, experiment.NMC, std_noise_test(Itest));    
end

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, std_noise_test, metric_fits./nominal_metric_matrix-1, '\sigma');
xlim([experiment.std_noise_min, experiment.std_noise_max])

results.noise_test.std_noise_test = std_noise_test;
results.noise_test.metric_mean = metric_mean;
results.noise_test.metric_std = metric_std;


%% bolus tests: arrival time
clear metric_mean metric_std metric_fits

Tarrival_test= linspace(experiment.Tarrival_min, experiment.Tarrival_max, Nexp_values);

for Itest = 1:length(Tarrival_test)
    t_test = t_input - (Tarrival_test(Itest)-experiment.Tarrival);
    input_function_test = interp1(t_input, input_function, t_test, 'linear', 0);
    Mz0_test = [sum(input_function)-sum(input_function_test) 0];
    
    Mxy = simulate_Nsite_model(Mz0_test, R1, [kPL 0], acq.flips, acq.TR, input_function_test);
    %     figure(99),  plot(t, Mxy), pause
    %
    
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips, experiment.NMC, experiment.std_noise);    
end

% figure
% plot_fit_results(fitting, Tarrival_test, metric_fits, 'Tarrival (s)');

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, Tarrival_test, metric_fits./nominal_metric_matrix-1, 'Tarrival (s)');
xlim([experiment.Tarrival_min, experiment.Tarrival_max])

results.Tarrival_test.Tarrival_test = Tarrival_test;
results.Tarrival_test.metric_mean = metric_mean;
results.Tarrival_test.metric_std = metric_std;

%% bolus tests: duration
clear metric_mean metric_std metric_fits

Tbolus_test = linspace(experiment.Tbolus_min, experiment.Tbolus_max, Nexp_values);


for Itest = 1:length(Tbolus_test)
    t_test = t*Tbolus/Tbolus_test(Itest);  % stretch/squeeze to modulate bolus
    input_function_test = interp1(t, input_function, t_test, 'linear', 0);
    input_function_test = input_function_test/sum(input_function_test) * sum(input_function); % normalized
    
    Mxy = simulate_Nsite_model(Mz0, R1, [kPL 0], acq.flips, acq.TR, input_function_test);
    %     figure(99),  plot(t, Mxy), pause
    %
    
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips, experiment.NMC, experiment.std_noise);    
    
end

% figure
% plot_fit_results(fitting, Tbolus_test, metric_fits, 'Tbolus (s)');
% 
subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, Tbolus_test, metric_fits./nominal_metric_matrix-1, 'Tbolus (s)');
xlim([experiment.Tbolus_min, experiment.Tbolus_max])

results.Tbolus_test.Tbolus_test = Tbolus_test;
results.Tbolus_test.metric_mean = metric_mean;
results.Tbolus_test.metric_std = metric_std;

% 
% subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
% [~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(Tbolus_test, kPL_fit./kPL-1, AUC_fit./AUC_predicted-1);
% ylim(ratio_limits), xlim([experiment.Tbolus_min, experiment.Tbolus_max]), xlabel('Tbolus (s)')
% 
% results.Tbolus_test.kPL_avg_error = mean(kPL_std);  % precision measurement
% results.Tbolus_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
% results.Tbolus_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement
% 
% results.Tbolus_test.AUC_avg_error = mean(AUC_std);  % precision measurement
% results.Tbolus_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
% results.Tbolus_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement
% 
% 
% 

%% T1 lactate tests - this is hard
clear metric_mean metric_std metric_fits

R1L_test = linspace(experiment.R1L_min, experiment.R1L_max, Nexp_values);

for Itest = 1:length(R1L_test)
    
    Mxy = simulate_Nsite_model(Mz0, [R1(1), R1L_test(Itest)], [kPL 0], acq.flips, acq.TR, input_function);
    
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips, experiment.NMC, experiment.std_noise);    
    
end

% figure
% plot_fit_results(fitting, 1./R1L_test, metric_fits, 'T_{1L} (s)');
subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, 1./R1L_test, metric_fits./nominal_metric_matrix-1, 'T_{1L} (s)');
xlim([1./experiment.R1L_max, 1./experiment.R1L_min])

results.R1L_test.R1L_test = R1L_test;
results.R1L_test.metric_mean = metric_mean;
results.R1L_test.metric_std = metric_std;


% subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
% [~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(1./R1L_test, kPL_fit/kPL-1, AUC_fit/AUC_predicted-1);
% ylim(ratio_limits), xlim(1./[experiment.R1L_max, experiment.R1L_min]), xlabel('T_{1L} (s)')
% 
% results.R1L_test.kPL_avg_error = mean(kPL_std);  % precision measurement
% results.R1L_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
% results.R1L_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement
% 
% results.R1L_test.AUC_avg_error = mean(AUC_std);  % precision measurement
% results.R1L_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
% results.R1L_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement


%% T1 pyruvate tests
clear metric_mean metric_std metric_fits

R1P_test = linspace(experiment.R1P_min, experiment.R1P_max, Nexp_values);

for Itest = 1:length(R1P_test)
    
    Mxy = simulate_Nsite_model(Mz0, [R1P_test(Itest), R1(2)], [kPL 0], acq.flips, acq.TR, input_function);
    
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips, experiment.NMC, experiment.std_noise);    
    
end

% figure
% plot_fit_results(fitting, 1./R1P_test, metric_fits, 'T_{1P} (s)');

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, 1./R1P_test, metric_fits./nominal_metric_matrix-1, 'T_{1P} (s)');
xlim([1./experiment.R1P_max, 1./experiment.R1P_min])

results.R1P_test.R1P_test = R1P_test;
results.R1P_test.metric_mean = metric_mean;
results.R1P_test.metric_std = metric_std;

% 
% subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
% [~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(1./R1P_test, kPL_fit/kPL-1, AUC_fit/AUC_predicted-1);
% ylim(ratio_limits), xlim(1./[experiment.R1P_max, experiment.R1P_min]), xlabel('T_{1P} (s)')
% 
% results.R1P_test.kPL_avg_error = mean(kPL_std);  % precision measurement
% results.R1P_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
% results.R1P_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement
% 
% results.R1P_test.AUC_avg_error = mean(AUC_std);  % precision measurement
% results.R1P_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
% results.R1P_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement
% 

%% B1 error tests
    % simulate inaccurate B1, & unknown

    clear metric_mean metric_std metric_fits

B1error_test = linspace(experiment.B1error_min, experiment.B1error_max, Nexp_values);

for Itest = 1:length(B1error_test)
    Mxy = simulate_Nsite_model(Mz0, R1, [kPL 0], acq.flips * (1+B1error_test(Itest)), acq.TR, input_function);
    
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips, experiment.NMC, experiment.std_noise);    
    
end

% figure
% plot_fit_results(fitting, B1error_test, metric_fits, '% B_1 error');

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, B1error_test, metric_fits./nominal_metric_matrix-1, '% B_1^+ error');
xlim([experiment.B1error_min, experiment.B1error_max])

results.B1error_test.B1error_test = B1error_test;
results.B1error_test.metric_mean = metric_mean;
results.B1error_test.metric_std = metric_std;


% 
% subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
% [~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(B1error_test, kPL_fit/kPL-1, AUC_fit/AUC_predicted-1);
% ylim(ratio_limits), xlim([experiment.B1error_min, experiment.B1error_max]), xlabel('% B_{1} error')
% 
% results.B1error_test.kPL_avg_error = mean(kPL_std);  % precision measurement
% results.B1error_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
% results.B1error_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement
%
% results.B1error_test.AUC_avg_error = mean(AUC_std);  % precision measurement
% results.B1error_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
% results.B1error_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement
%

%% B1 difference tests
%     % simulate inaccurate B1, but known
% clear metric_mean metric_std metric_fits
%
B1diff_test = linspace(experiment.B1diff_min, experiment.B1diff_max, Nexp_values);

for Itest = 1:length(B1diff_test)
    Mxy = simulate_Nsite_model(Mz0, R1, [kPL 0], acq.flips * (1+B1diff_test(Itest)), acq.TR, input_function);
    
    for f = 1:N_fitting_methods
        switch func2str(fitting(f).fit_fcn)
            case {'compute_AUCratio', 'compute_mean_time', 'compute_TTP'} % will use lactate mean time or TTP
                
            otherwise
                fitting(f).input_arguments = {acq.TR, acq.flips* (1+B1diff_test(Itest)), fitting(f).params_fixed, fitting(f).params_est, [], 0};
        end
        nominal_metric_B1diff_test(f,Itest) = compute_nominal_metric_value(fitting(f).fit_fcn, Mxy, kPL);
    end
    
    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, acq.TR, acq.flips* (1+B1diff_test(Itest)), experiment.NMC, experiment.std_noise);
    
end

% reset fitting.input_arguments
for f = 1:N_fitting_methods
    switch func2str(fitting(f).fit_fcn)
        case {'compute_AUCratio', 'compute_mean_time', 'compute_TTP'} % will use lactate mean time or TTP
            
        otherwise
            fitting(f).input_arguments = {acq.TR, acq.flips, fitting(f).params_fixed, fitting(f).params_est, [], 0};
    end
end
    

% figure
% plot_fit_results(fitting, B1diff_test, metric_fits, '% B_1 change');

subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
plot_fit_results_normalized(fitting, B1diff_test, ...
    metric_fits./repmat( reshape(nominal_metric_B1diff_test, [N_fitting_methods 1 Nexp_values]), [1, experiment.NMC 1])-1, '% B_1^+ change')
xlim([experiment.B1diff_min, experiment.B1diff_max])

results.B1diff_test.B1diff_test = B1diff_test;
results.B1diff_test.metric_mean = metric_mean;
results.B1diff_test.metric_std = metric_std;

% 
% subplot(Nplot1, Nplot2, Iplot); Iplot = Iplot+1;
% [~,kPL_mean,AUC_mean,kPL_std,AUC_std]=plot_with_mean_and_std(B1diff_test, kPL_fit/kPL-1, AUC_fit./repmat(AUC_predicted_test(:), [1, NMC]) -1);
% ylim(ratio_limits), xlim([experiment.B1diff_min, experiment.B1diff_max]), xlabel('% B_{1} difference')
% 
% 
% 
% results.B1diff_test.kPL_avg_error = mean(kPL_std);  % precision measurement
% results.B1diff_test.kPL_avg_bias = mean(abs(kPL_mean));  % accuracy measurement
% results.B1diff_test.kPL_std_bias = std(kPL_mean);  % accuracy measurement
% 
% results.B1diff_test.AUC_avg_error = mean(AUC_std);  % precision measurement
% results.B1diff_test.AUC_avg_bias = mean(abs(AUC_mean));  % accuracy measurement
% results.B1diff_test.AUC_std_bias = std(AUC_mean);  % accuracy measurement
% 
%

warning(warnStruct);

%% HELPER FUNCTIONS

    function [h,Y1,Y2,DELTA1,DELTA2]=plot_with_mean_and_std(x, y1, y2)
        
        h = gca;
        
        DELTA1 = std(y1, [], 2);
        DELTA2 = std(y2, [], 2);
        Y1 = mean(y1,2);
        Y2 = mean(y2,2);
        plot(x, Y1, 'b-', x, Y2, 'g-', ...  % means
            x, Y1+DELTA1, 'b--', x, Y1-DELTA1, 'b--', ... % stds
            x, Y2+DELTA2, 'g--', x, Y2-DELTA2, 'g--')
    end

    function plot_fit_results(fitting, xvalues, metric_fits, xlabel_string)
        for f = 1:N_fitting_methods
            subplot(1, N_fitting_methods, f)
            shadedErrorBar(xvalues, squeeze(metric_fits(f,:,:)), {@mean, @std})
            xlabel(xlabel_string), ylabel(fitting(f).metric), title(fitting(f).fit_description)
        end
        
    end

    function plot_fit_results_normalized(fitting, xvalues, metric_fits, xlabel_string)
        
        lineprops = {'b', 'g', 'r', 'c', 'm', 'y'};
        for f = 1:N_fitting_methods
            %            shadedErrorBar(xvalues, squeeze(metric_fits(f,:,:)), {@mean, @std}, 'lineprops', lineprops{f})
            ymean = mean(squeeze(metric_fits(f,:,:)));
            ystd = std(squeeze(metric_fits(f,:,:)));
            plot(xvalues, ymean, lineprops{f}, 'DisplayName', fitting(f).fit_description)
            if f == 1
                ylim(ratio_limits)
                hold on
            end
                plot(xvalues, ymean+ystd, [lineprops{f} '--'], ...
                xvalues, ymean-ystd, [lineprops{f} '--'],'HandleVisibility','off')

        end
        
        xlabel(xlabel_string)
        ylabel('relative error')
        %legend?
        hold off
    end

    function nominal_metric_value = compute_nominal_metric_value(fit_fcn, Mxy, kPL)
        switch func2str(fit_fcn)
            case 'compute_AUCratio'
                nominal_metric_value = compute_AUCratio(Mxy);
            case 'compute_mean_time'
                nominal_metric_value = compute_mean_time(Mxy(2,:), acq.TR);
            case 'compute_TTP'
                nominal_metric_value = compute_TTP(Mxy(2,:), acq.TR);
            otherwise
                nominal_metric_value = kPL;
        end
        
    end

    function [metric_mean, metric_std, metric_fits] = fitting_simulation(fitting, Mxy, TR, flips, NMC, std_noise);
        
        metric_fits = zeros(N_fitting_methods,NMC);
        
        for n = 1:NMC
            Sn = Mxy + randn(size(Mxy))*std_noise;
            
            for f = 1:N_fitting_methods
                params_fit = fitting(f).fit_fcn(Sn, fitting(f).input_arguments{:});
                switch func2str(fitting(f).fit_fcn)
                    case 'compute_AUCratio'
                        metric_fits(f,n) = params_fit;
                    case {'compute_mean_time', 'compute_TTP'} % will use lactate mean time or TTP
                        metric_fits(f,n) = params_fit(2);
                    otherwise  % for fit_* functions
                        metric_fits(f,n) = params_fit.(fitting(f).metric);
                end
            end
        end
        
        metric_mean = mean(metric_fits, 2);
        metric_std = std(metric_fits, [], 2);
        
    end

end