% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 36; TR = 2; N = Tacq/TR;
R1_aKG_C1 = 1/20; R1_aKG_C5 = 1/20; R1_2HG = 1/25;
k_aKG_2HG = 0.01 ;

ratio_aKG_C1toC5 = 10;
k_aKG_C1toC5 = 10;  % relatively fast exchange, this is a guess
k_aKG_C5toC1 = k_aKG_C1toC5 *ratio_aKG_C1toC5;   %

k_all = [k_aKG_C1toC5 k_aKG_2HG  % aKG C1 rates to C5 and 2HG
    k_aKG_C5toC1  k_aKG_2HG  % aKG C5 rates to C1 and 2HG
    0 0];  % 2HG rates to aKG

std_noise = 0.002;
Nmets = 3;  % aKG_C1, aKG_C5, 2HG

input_condition = 1; % choose from various simulated starting conditions
switch input_condition
    case 1
        % gamma variate input function - most realistic
        Tarrival = 0;  Tbolus = 8;
        input_function = realistic_input_function(N, TR, Tarrival, Tbolus);
        Mz0 = zeros(1,Nmets);
    case 2
        % boxcar input function
        Tbolus = 8;  Tarrival = 0;
        Ibolus = [1:round(Tbolus/TR)] + round(Tarrival/TR);
        Rinj = 1/Tbolus; % normalize for a total magnetization input = 1
        input_function = zeros(1,N);
        Mz0 = zeros(1,Nmets);  input_function(Ibolus) =  Rinj*TR;
end

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = ones(2,N)*30*pi/180;  flip_descripton{1} = 'constant, single-band';
flips(1:2,1:N,2) = repmat([10;35]*pi/180,[1 N]);  flip_descripton{2} = 'constant, multi-band';

k12 = 0.05; % for variable flip angle designs
flips(1:2,1:N,3) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective pyruvate variable flip angles
    vfa_opt_signal(N, exp(-TR * ( R1_2HG)))]; % max lactate SNR variable flip angle
flip_descripton{3} = 'max product SNR variable flip, multi-band';
% flips(1:2,1:N,4) = repmat(vfa_const_amp(N, pi/2), [2,1]);  % RF compensated variable flip angle
% flips(1:2,1:N,5) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective variable flip angle
%     vfa_const_amp(N, pi/2, exp(-TR * ( - k12)))];
% flips(1:2,1:N,6) = [vfa_const_amp(N, pi/2, exp(-TR * (k12))); ... % saturation recovery
%     ones(1,N)*pi/2];

% dimension 1: aKG_C1, aKG_C5, 2HG
% put same flips for C5 and 2HG
flips = cat(1, flips, flips(2,:,:));  % duplicate

% input includes both C1 and C5 label
input_function = [input_function; input_function/ratio_aKG_C1toC5 ];

N_flip_schemes = size(flips,3);

t = [0:N-1]*TR + Tin;

figure
subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('Product (Lactate, Alanine, Bicarbonate) flips')
legend(flip_descripton)

flip_description_array = [repmat('    ',N_flip_schemes,1),  char(flip_descripton(:))];


% generate simulated data
plot_simulated_data = 1;
noise_S = randn([Nmets-1 N])*std_noise;  % same noise for all flip schedules
for Iflips = 1:N_flip_schemes
    [Mxy(1:Nmets, 1:N, Iflips), Mz] = simulate_aKG_model(Mz0, [R1_aKG_C1 R1_aKG_C5, R1_2HG], k_all , flips(:,:,Iflips), TR, input_function);

    % data has overlap of aKG C5 and 2HG resonances
    Mxy_overlapped(1:Nmets-1, 1:N, Iflips) = [Mxy(1, 1:N, Iflips); ...
        Mxy(2, 1:N, Iflips) + Mxy(3, 1:N, Iflips)];
    % add noise
    Sn(1:Nmets-1, 1:size(Mxy,2),  Iflips) = Mxy_overlapped(:,:,Iflips) + noise_S;

    if plot_simulated_data
        figure(99)
        subplot(311)
        plot(t, Mxy(:,:,Iflips))
        subplot(312)
        plot(t, Mxy_overlapped(:,:,Iflips))
        subplot(313)
        plot(t, Sn(:,:,Iflips))
        pause
    end
end




%%
% initial parameter guesses
R1A_est = 1/25; R1H_est = 1/25;
kA2_est = .02;
ratio_aKG_C1toC5_est = 10;
plot_fits = 1;
fit_function = @fit_aKG_kinetics;

%% Test fitting - 2-site
disp('2-site model: pyruvate -> lactate')
disp('Fitting kPL with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1A_est; params_fixed.R1L = R1H_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1,1:size(Mxy,2),  Iflips)] = fit_function(Mxy(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1,1:size(Mxy,2),  Iflips)] = fit_function(Sn(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1,1:size(Mxy,2),  Iflips)] = fit_function(abs(Sn(1:2,:,Iflips)), TR, flips(1:2,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL  '])
disp(num2str(k_aKG_2HG,2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([1;1 + Nparams_fit;1 + 2*Nparams_fit]); 
disp(['KPL   '])
disp([num2str(k_fit.',2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([1;1 + Nparams_fit;1 + 2*Nparams_fit]); 
disp(['KPL   '])
disp([num2str(k_fit.',2), flip_description_array])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([1;1 + Nparams_fit;1 + 2*Nparams_fit]); 
disp(['KPL   '])
disp([num2str(k_fit.',2), flip_description_array])

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)

disp('Press any key to continue')
disp(' ')

pause

%% Test fitting - 2-site plus T1 lactate fitting
disp('2-site model: pyruvate -> lactate')
disp('Fitting kPL  and T1 of lactate')
disp('Fitting both parameters leads to increases in variability, which can be')
disp('alleviated to some extent by constraints on the parameter values')
disp('(Pyruvate T1 values have very small effect on inputless fitting approach)')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1A_est;
params_est.kPL = kPL_est;  params_est.R1L = R1H_est;
% set constraints on lactate T1:
params_est.R1L_lb = 1/40;
params_est.R1L_ub = 1/15;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1,1:size(Mxy,2),  Iflips)] = fit_function(Mxy(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1,1:size(Mxy,2),  Iflips)] = fit_function(Sn(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1,1:size(Mxy,2),  Iflips)] = fit_function(abs(Sn(1:2,:,Iflips)), TR, flips(1:2,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL      R1L '])
disp(num2str([k_aKG_2HG R1_2HG],2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      R1L '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      R1L '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      R1L '])
disp([num2str(k_fit,2), flip_description_array])

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)

disp('Press any key to continue')
disp(' ')

pause

