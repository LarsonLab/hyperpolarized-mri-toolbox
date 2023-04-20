clear all

experiment.NMC = 500;  % less for quicker testing

flip_scheme = 4;  % see below

% default experiment values

experiment.R1P = 1/25;  experiment.R1L =1/25;  experiment.kPL = 0.02; experiment.std_noise = 0.015;
experiment.Tarrival = 4;  % negative means bolus arrival before acquisition starts
experiment.Tbolus = 10;

        disp('Running Monte Carlo Simulation')
        fit_description = [];
        
        clear params_est params_fixed acq fitting
        
        % default fitting parameters
        R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;
        params_fixed.R1P = R1P_est;
        params_est.kPL = kPL_est;
            fit_description = [fit_description, '  Fixed T1L'];
            params_fixed.R1L = R1L_est;

 %       params_fixed.L0_start = 0;  % ok to be free parameter?

        % allowing for fit of R1L increases variability substantially
        % R1P minimal change
        % constraining R1L to narrow range maybe reasonable compromise
        
           fit_description = [fit_description, '  Inputless fitting'];
           fitting(1).fit_fcn = @fit_pyr_kinetics;
        fitting(1).metric = 'kPL';
        
                Tin = 0; Tacq = 60; acq.TR = 3; acq.N = Tacq/acq.TR;
                acq.flips = repmat([10*pi/180; 40*pi/180], [1 acq.N]);
        
        
        fitting(1).params_est = params_est; fitting(1).params_fixed = params_fixed;
        fitting(1).fit_description = fit_description;
        

%%
% simulation and plotting parameters
Nexp_values = 8;
ratio_limits = [-.2 .2];

% fitting
N_fitting_methods = length(fitting);

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

% hdata = figure;
% subplot(121) , plot(t, squeeze(results.sample_data(1,:))/experiment.std_noise), title('sample simulated pyruvate')
% xlabel('time (s)'), ylabel('Signal (SNR)')
% subplot(122) , plot(t, squeeze(results.sample_data(2,:))/experiment.std_noise), title('sample simulated lactate')
% xlabel('time (s)'), ylabel('Signal (SNR)')
% 
%% setup for plots
for f = 1:N_fitting_methods
    legend_description{f} = fitting(f).fit_description;
    fitting(f).nominal_metric_value = kPL; %compute_nominal_metric_value(fitting(f).fit_fcn, Mxy, kPL);
    nominal_metric(f) = kPL; %compute_nominal_metric_value(fitting(f).fit_fcn, Mxy, kPL);
    
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

%% TR test
clear metric_mean metric_std metric_fits

TR_test = linspace(1, 10, Nexp_values).';
        metric_fits = zeros(N_fitting_methods,experiment.NMC, Nexp_values);

for Itest = 1:length(TR_test)

                    Tin = 0; Tacq = 60; acq.TR = TR_test(Itest); acq.N = ceil(Tacq/acq.TR);
                    switch flip_scheme
                        case 1
                acq.flips = [ 10*pi/180*ones(1,acq.N); 40*pi/180*ones(1,acq.N)];
                flip_description = 'pyr=10, lac=40';
                        case 2
                acq.flips = [ 10*pi/180*ones(1,acq.N); vfa_opt_signal(acq.N, exp(-acq.TR/25))];
                flip_description = 'pyr=10, lac=optimal SNR';
                        case 3
                acq.flips = [ acos(cos(10*pi/180)^(20/acq.N))*ones(1,acq.N); vfa_opt_signal(acq.N, exp(-acq.TR/25))];
                flip_description = 'pyr=equal to TR=3, 10, lac=optimal SNR';
                        case 4
                acq.flips = [ acos(cos(10*pi/180)^(20/acq.N))*ones(1,acq.N); acos(cos(40*pi/180)^(20/acq.N))*ones(1,acq.N)];
                flip_description = 'pyr=equal to 10, lac=equal to 40 for TR=3';
                    end
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

    Mxy = simulate_Nsite_model(Mz0, R1, [kPL 0], acq.flips, TR_test(Itest), input_function);
%    [metric_mean(:,Itest), metric_std(:, Itest) metric_fits(:,:,Itest)] = fitting_simulation(fitting,Mxy, TR_test(Itest), acq.flips, experiment.NMC, experiment.std_noise);

fitting(1).input_arguments = {TR_test(Itest), acq.flips, fitting(f).params_fixed, fitting(f).params_est, [], 0};
        
        for n = 1:experiment.NMC
            Sn = Mxy + randn(size(Mxy))*experiment.std_noise;
            
            for f = 1:N_fitting_methods
                params_fit = fitting(f).fit_fcn(Sn, fitting(f).input_arguments{:});
                switch func2str(fitting(f).fit_fcn)
                    case 'compute_AUCratio'
                        metric_fits(f,n,Itest) = params_fit;
                    case {'compute_mean_time', 'compute_TTP'} % will use lactate mean time or TTP
                        metric_fits(f,n,Itest) = params_fit(2);
                    otherwise  % for fit_* functions
                        metric_fits(f,n,Itest) = params_fit.(fitting(f).metric);
                end
            end
        end
        
        metric_mean = mean(metric_fits, 2);
        metric_std = std(metric_fits, [], 2);
 
    
end

hsim{1} = figure;
%plot_fit_results(fitting, TR_test, metric_fits, 'TR (s)');
        for f = 1:N_fitting_methods
            subplot(1, N_fitting_methods, f)
            shadedErrorBar(TR_test, squeeze(metric_fits(f,:,:)), {@mean, @std})
            xlabel('TR (s)'), ylabel(fitting(f).metric), title([ fitting(f).fit_description ' ' flip_description])
        end

        return
        
hsim{2} = figure;
%plot_fit_results_normalized(fitting, TR_test, metric_fits./nominal_metric_matrix-1, 'TR (s)');
metric_fits_plot = metric_fits./nominal_metric_matrix-1;
        lineprops = {'b', 'g', 'r', 'c', 'm', 'y'};
        for f = 1:N_fitting_methods
            %            shadedErrorBar(xvalues, squeeze(metric_fits(f,:,:)), {@mean, @std}, 'lineprops', lineprops{f})
            ymean = mean(squeeze(metric_fits_plot(f,:,:)));
            ystd = std(squeeze(metric_fits_plot(f,:,:)));
            plot(TR_test, ymean, lineprops{f}, 'DisplayName', fitting(f).fit_description)
            if f == 1
                ylim(ratio_limits)
                hold on
            end
                plot(TR_test, ymean+ystd, [lineprops{f} '--'], ...
                TR_test, ymean-ystd, [lineprops{f} '--'],'HandleVisibility','off')

        end
        
        xlabel('TR (s)')
        ylabel('relative error')
        %legend?
        hold off


    legh = legend;
    legh.Position = [0.35 0.01 0.3 0.1];
