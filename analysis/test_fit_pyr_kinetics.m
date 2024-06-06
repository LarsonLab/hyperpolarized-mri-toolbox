% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; R1B = 1/15; R1A = 1/25;
kPL = 0.05; kPB = 0.03; kPA = 0.02;
std_noise = 0.005;

input_condition = 1; % choose from various simulated starting conditions
switch input_condition
    case 1
        % gamma variate input function - most realistic
        Tarrival = 0;  Tbolus = 12;
        input_function = realistic_input_function(N, TR, Tarrival, Tbolus);
        Mz0 = [0,0,0,0]; 
    case 2
        % boxcar input function
        Tbolus = 12;  Tarrival = 0;
        Ibolus = [1:round(Tbolus/TR)] + round(Tarrival/TR);
        Rinj = 1/Tbolus; % normalize for a total magnetization input = 1
        Mz0 = [0,0,0,0]; 
        input_function = zeros(1,N);
        input_function(Ibolus) =  Rinj*TR;
    case 3
        Mz0 = [1.5,0,0,0]; % no input function
        input_function = [];
    case 4
        Tin = 6; Mz0 = Tin; % no input function, delayed start
        input_function = [];
end

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = ones(2,N)*30*pi/180;  flip_descripton{1} = 'constant, single-band';
flips(1:2,1:N,2) = repmat([20;35]*pi/180,[1 N]);  flip_descripton{2} = 'constant, multi-band';

k12 = 0.05; % for variable flip angle designs
flips(1:2,1:N,3) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective pyruvate variable flip angles
    vfa_opt_signal(N, exp(-TR * ( R1L)))]; % max lactate SNR variable flip angle
flip_descripton{3} = 'max product SNR variable flip, multi-band';
% flips(1:2,1:N,4) = repmat(vfa_const_amp(N, pi/2), [2,1]);  % RF compensated variable flip angle
% flips(1:2,1:N,5) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective variable flip angle
%     vfa_const_amp(N, pi/2, exp(-TR * ( - k12)))];
% flips(1:2,1:N,6) = [vfa_const_amp(N, pi/2, exp(-TR * (k12))); ... % saturation recovery
%     ones(1,N)*pi/2];

flips = cat(1, flips, repmat( flips(2,:,:), [2 1 1]));

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
noise_S = randn([4 N])*std_noise;  % same noise for all flip schedules
for Iflips = 1:N_flip_schemes
    [Mxy(1:4, 1:N, Iflips), Mz] = simulate_Nsite_model(Mz0, [R1P R1L R1B R1A], [kPL 0; kPB 0; kPA 0], flips(:,:,Iflips), TR, input_function);
    % add noise
    Sn(1:size(Mxy,1), 1:size(Mxy,2),  Iflips) = Mxy(:,:,Iflips) + noise_S;
end

% initial parameter guesses
R1P_est = 1/25; R1L_est = 1/25; R1B_est = 1/15; R1A_est = 1/25;
kPL_est = .02; kPB_est = .02; kPA_est = .02;
plot_fits = 1;
fit_function = @fit_pyr_kinetics;

%% Test fitting - 2-site
disp('2-site model: pyruvate -> lactate')
disp('Fitting kPL with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
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
disp(num2str(kPL,2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/N_flip_schemes;
params_index = [0:N_flip_schemes-1]*Nparams_fit + 1;
k_fit = k_fit(params_index); 
disp(['KPL   '])
disp([num2str(k_fit.',2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit(params_index); 
disp(['KPL   '])
disp([num2str(k_fit.',2), flip_description_array])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit(params_index); 
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
params_fixed.R1P = R1P_est;
params_est.kPL = kPL_est;  params_est.R1L = R1L_est;
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
disp(num2str([kPL R1L],2))

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

%% Test fitting - fit kPX only
disp('4-site model: pyruvate -> lactate, bicarbonate, alanine')
disp('Fitting kPL, kPB, and kPA, with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:3,1:size(Mxy,2),  Iflips)] = fit_function(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:3,1:size(Mxy,2),  Iflips)] = fit_function(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
 %   [params_fitn_mag(:,Iflips) Snfit_mag(1:3,1:size(Mxy,2),  Iflips)] = fit_function(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL      KPB      KPA    '])
disp(num2str([kPL, kPB, kPA],2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:3];[1:3] + Nparams_fit;[1:3] + 2*Nparams_fit]); 
disp(['KPL      KPB      KPA    '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:3];[1:3] + Nparams_fit;[1:3] + 2*Nparams_fit]); 
disp(['KPL      KPB      KPA    '])
disp([num2str(k_fit,2), flip_description_array])
% disp('Noisy magnitude fit results:')
% k_fit = struct2array(params_fitn_mag);
% k_fit = k_fit([[1:3];[1:3] + Nparams_fit;[1:3] + 2*Nparams_fit]); 
% disp(['KPL      KPB      KPA    '])
% disp([num2str(k_fit,2), flip_description_array])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
%plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits')% (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)
subplot(223) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
%plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')
subplot(224) , plot(t, squeeze(Sn(4,:,:)))
hold on, plot(t, squeeze(Snfit_complex(3,:,:)),':')
%plot(t, squeeze(Snfit_mag(3,:,:)),'--'), hold off
title('Alanine signals and fits')

disp('Press any key to continue')
disp(' ')

pause
%% Test fitting - 3-site
disp('3-site model: pyruvate -> lactate, bicarbonate')
disp('Fitting kPL, and kPB, with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:2,1:size(Mxy,2),  Iflips)] = fit_function(Mxy(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:2,1:size(Mxy,2),  Iflips)] = fit_function(Sn(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:2,1:size(Mxy,2),  Iflips)] = fit_function(abs(Sn(1:3,:,Iflips)), TR, flips(1:3,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL      KPB   '])
disp(num2str([kPL, kPB],2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      KPB   '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      KPB   '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['KPL      KPB   '])
disp([num2str(k_fit,2), flip_description_array])

figure
subplot(131) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(132) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)
subplot(133) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')

return
disp('Press any key to continue')
disp(' ')

pause

