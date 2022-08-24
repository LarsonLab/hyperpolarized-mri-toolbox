% The goal of this notebook is to demonstrate simulation tools for 
% comparing metrics of hyperpolarized MR metabolic data.  
% It is based on simulations originally presented in the following citation(s)

%IT WORKS 
%%  Sample Data for model

acq_sample.TR = 3;
acq_sample.N = 20;
acq_sample.flips = repmat([25*pi/180; 25*pi/180], [1 acq_sample.N]);
experiment_sample.R1P = 1/25;
experiment_sample.R1L = 1/25;
experiment_sample.KPL = 0.02;
experiment_sample.Tarrival = 0; 
experiment_sample.Tbolus = 8;

[input_function, t_input] = realistic_input_function(acq_sample.N, acq_sample.TR,experiment_sample.Tarrival,experiment_sample.Tbolus);

[Mxy, Mz] = simulate_Nsite_model([0,0], [experiment_sample.R1P experiment_sample.R1L],[experiment_sample.KPL 0], acq_sample.flips,acq_sample.TR,input_function);

plot(t_input,Mxy)
xlabel('time (s)'),ylabel('Signal')
legend('pyruvate','lactate')
title('Simulated Data(without noise)')

%% Monte Carlo Sim
clear experiment acq

%default experiment values
experiment.R1P = 1/25;
experiment.R1L = 1/25;
experiment.KPL = 0.02;
experiment.std_noise = 0.005;
experiment.Tarrival = 0; 
experiment.Tboulus = 8;

%Acquisition variables 
tacq = 48;
acq.TR = 3;
acq.N = tacq/acq.TR;

clear params_est params_fixed 

%default fitting parameters variables
R1P_est = 1/25; 
R1L_est = 1/25;
KPL_est = .02;
Tarrival_est = experiment.Tarrival;
Tbolus_est = experiment.Tboulus;
Rinj_est = 0.1;

params_fixed.R1P = R1P_est;
params_fixed.R1L = R1L_est;
params_est.KPL = KPL_est;

%THese strucutures choose the fitting methods to test
clear fitting

fitting(1).fit_fcn = @fit_pyr_kinetics;
fitting(1).params_fixed = params_fixed;
fitting(1).params_est = params_est;
fitting(1).fit_description = 'Inputless fitting';
fitting(1).metric = 'kPL';

fitting(2).fit_fcn = @fit_pyr_kinetics;

clear params_est params_fixed
params_fixed.R1P = R1P_est;
params_est.R1L = R1L_est;
params_est.KPl = KPL_est;

fitting(2).params_fixed = params_fixed;
fitting(2).params_est = params_est;
fitting(2).fit_description = 'Inputless fitting with T1 fit';
fitting(2).metric = 'kPL';

fitting(3).fit_fcn = @compute_AUCratio;
fitting(3).metric = 'AUCratio';
fitting(3).fit_description = 'AUC Ratio';

%Run the sim
experiment.NMC = 25;

% Options to choose from three different flip angle schemes
flip_scheme = 2;
switch flip_scheme
    case 1
        acq.flips = repmat([25*pi/180;25*pi/180], [1 acq.N]);
        flip_description = 'constant 25-degrees';
    case 2
        acq.flips = repmat([20*pi/180; 30*pi/180], [1 acq.N]);
        flip_description = 'metabolite-specific 20(pyruvate)/30(lactate)-degrees';
    case 3
        k12 = 0.05; % for variable flip angle designs
        acq.flips = [vfa_const_amp(acq.N, pi/2, exp(-acq.TR * ( k12))); ... % T1-effective pyruvate variable flip angles
            vfa_opt_signal(acq.N, exp(-acq.TR * ( R1L_est)))]; % max lactate SNR variable flip angle
        flip_description = 'metabolite-specific variable flip angle';
end

disp(['flip angle scheme: ' flip_description])

[results, hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, experiment );

