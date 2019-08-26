% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; R1B = 1/15; R1A = 1/25;
kPL = 0.05; kPB = 0.03; kPA = 0.02;
std_noise = 0.005;

% gamma variate input function - most realistic
t = [1:N]*TR;
Tarrival = 0;
A = 4; B = 1*TR;
input_function = gampdf(t-Tarrival,A,B);  % gamma distribution -continued input
input_function = input_function/sum(input_function);% normalize for a total magnetization input = 1
Mz0 = zeros(4,1);

% Test over multiple combinations of flip angle schemes
k12 = 0.03; % for variable flip angle designs
flips(1:2,1:N,1) = ones(2,N)*30*pi/180;  % constant, single-band
flips(1:2,1:N,2) = repmat([20;35]*pi/180,[1 N]);  % constant, multi-band
flips(1:2,1:N,3) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % max lactate SNR variable flip angle
    vfa_opt_signal(N, exp(-TR * ( R1L)))];

flips = cat(1, flips, repmat( flips(2,:,:), [2 1 1]));

N_flip_schemes = size(flips,3);

t = [0:N-1]*TR + Tin;

figure
subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('Product (Lactate, Alanine, Bicarbonate) flips')
legend('constant','multiband',  'max lactate SNR vfa'); %, 'vfa', 'T1-effective vfa', 'Saturation Recovery')

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

%% Test fitting - fit kPX only
disp('Fitting kPL, kPB, and kPA, with fixed relaxation rates:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
 %   [params_fitn_mag(:,Iflips) Snfit_mag(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL  = ' num2str(kPL)])
disp(['KPB  = ' num2str(kPB)])
disp(['KPA  = ' num2str(kPA)])

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:3],[1:3] + Nparams_fit,[1:3] + 2*Nparams_fit]); 
disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:3],[1:3] + Nparams_fit,[1:3] + 2*Nparams_fit]); 
disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])
% disp('Noisy magnitude fit results:')
% k_fit = struct2array(params_fitn_mag);
% k_fit = k_fit([[1:3],[9:11],[17:19]]); 
% disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
%plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits')% (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')
subplot(223) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
%plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')
subplot(224) , plot(t, squeeze(Sn(4,:,:)))
hold on, plot(t, squeeze(Snfit_complex(3,:,:)),':')
%plot(t, squeeze(Snfit_mag(3,:,:)),'--'), hold off
title('Alanine signals and fits')

return
%% Test fitting - 3-site

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:2,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Mxy(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:2,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Sn(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:2,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(abs(Sn(1:3,:,Iflips)), TR, flips(1:3,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL  = ' num2str(kPL)])
disp(['KPB  = ' num2str(kPB)])

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:2],[1:2] + Nparams_fit,[1:2] + 2*Nparams_fit]); 
disp([['KPL  = ';'KPB  = '],num2str(reshape(k_fit,[2 3]),2)])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:2],[1:2] + Nparams_fit,[1:2] + 2*Nparams_fit]); 
disp([['KPL  = ';'KPB  = '],num2str(reshape(k_fit,[2 3]),2)])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1:2],[1:2] + Nparams_fit,[1:2] + 2*Nparams_fit]); 
disp([['KPL  = ';'KPB  = '],num2str(reshape(k_fit,[2 3]),2)])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')
subplot(223) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')

%% Test fitting - 2-site

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Mxy(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Sn(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(abs(Sn(1:2,:,Iflips)), TR, flips(1:2,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL  = ' num2str(kPL)])

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([1,1 + Nparams_fit,1 + 2*Nparams_fit]); 
disp([['KPL  = '],num2str(reshape(k_fit,[1 3]),2)])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([1,1 + Nparams_fit,1 + 2*Nparams_fit]); 
disp([['KPL  = '],num2str(reshape(k_fit,[1 3]),2)])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([1,1 + Nparams_fit,1 + 2*Nparams_fit]); 
disp([['KPL  = '],num2str(reshape(k_fit,[1 3]),2)])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')
