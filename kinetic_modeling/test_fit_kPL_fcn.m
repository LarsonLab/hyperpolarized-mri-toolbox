% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; KPL = 0.05; std_noise = 0.01;
k12 = 0.05; % for variable flip angle designs
input_function = zeros(1,N);
Mz0 = [0,0];  input_function(1:6) = gampdf([1:6],4,1)*3;  % gamma variate input function
%Mz0 = [0,0]; input_function(1:4) =  1; % boxcar input function
%Mz0 = [1,0]; % no input function


% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = ones(2,N)*30*pi/180;  % constant, single-band
flips(1:2,1:N,2) = repmat([20;35]*pi/180,[1 N]);  % constant, multi-band
flips(1:2,1:N,3) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % max lactate SNR variable flip angle
    vfa_opt_signal(N, exp(-TR * ( R1L)))];
% flips(1:2,1:N,4) = repmat(vfa_const_amp(N, pi/2), [2,1]);  % RF compensated variable flip angle
% flips(1:2,1:N,5) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective variable flip angle
%     vfa_const_amp(N, pi/2, exp(-TR * ( - k12)))];
% flips(1:2,1:N,6) = [vfa_const_amp(N, pi/2, exp(-TR * (k12))); ... % saturation recovery
%     ones(1,N)*pi/2];

t = [0:N-1]*TR + Tin;

figure
subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('Lactate flips')
legend('constant','multiband',  'max lactate SNR vfa'); %, 'vfa', 'T1-effective vfa', 'Saturation Recovery')

% generate simulated data
for Iflips = 1:size(flips,3)
    [Mxy(1:2, 1:N, Iflips), Mz] = simulate_2site_model(Mz0, [R1P R1L], [KPL 0], flips(:,:,Iflips), TR, input_function);
    % add noise
    Sn(1:size(Mxy,1), 1:size(Mxy,2),  Iflips) = Mxy(:,:,Iflips) + randn([size(Mxy,1), size(Mxy,2)])*std_noise;
end

%% Test fitting - fit kPL and T1 of lactate
R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;
params_fixed.R1P = R1P_est;
params_est.kPL = kPL_est; params_est.R1L = R1L_est;

plot_fits = 0;

for Iflips = 1:size(flips,3)
    % no noise
    [params_fit(:,Iflips) Sfit(1:size(Mxy,2),  Iflips)] = fit_kPL(Mxy(:, :, Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);

    % noise, complex-valued fitting
    [params_fitn_complex(:,Iflips) Snfit_complex(1:size(Mxy,2),  Iflips)] = fit_kPL(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);

    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:size(Mxy,2),  Iflips)] = fit_kPL(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P, R1L, KPL))
disp('Noiseless fit results:')
disp(['KPL         R1L   = ']); disp(num2str(reshape(struct2array(params_fit), 2, size(flips,3)).'))
disp('Noisy complex fit results:')
disp(['KPL         R1L   = ']); disp(num2str(reshape(struct2array(params_fitn_complex), 2, size(flips,3)).'))
disp('Noisy magnitude fit results:')
disp(['KPL         R1L   = ']); disp(num2str(reshape(struct2array(params_fitn_mag), 2, size(flips,3)).'))

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(:,:)),':')
plot(t, squeeze(Snfit_mag(:,:)),'--')
title('Lactate signals and fits (dashed)')
legend('constant','multiband', 'vfa', 'T1-effective vfa', 'max lactate SNR vfa', 'Saturation Recovery')

pause(2)
%% Test fitting - fit kPL only
clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est; 

for Iflips = 1:size(flips,3)
    % no noise
    [params_fit(:,Iflips) Sfit(1:size(Mxy,2),  Iflips)] = fit_kPL(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);

    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:size(Mxy,2),  Iflips)] = fit_kPL(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);

    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:size(Mxy,2),  Iflips)] = fit_kPL(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P, R1L, KPL))
disp('Noiseless fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fit).'))
disp('Noisy complex fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fitn_complex).'))
disp('Noisy magnitude fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fitn_mag).'))

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(:,:)),':')
plot(t, squeeze(Snfit_mag(:,:)),'--')
title('Lactate signals and fits (dashed)')
legend('constant','multiband', 'vfa', 'T1-effective vfa', 'max lactate SNR vfa', 'Saturation Recovery')
