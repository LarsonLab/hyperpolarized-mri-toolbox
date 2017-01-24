% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1 = [1/25 1/25]; KPL = 0.05; std_noise = 0.01;
k12 = 0.05; % for variable flip angle designs
input_function = zeros(1,N); 
Mz0 = [0,0];  input_function(1:6) = gampdf([1:6],4,1)*3;  % gamma variate input function
%Mz0 = [0,0]; input_function(1:4) =  1; % boxcar input function
%Mz0 = [1,0]; % no input function


% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = ones(2,N)*30*pi/180;  % constant, single-band
flips(1:2,1:N,2) = repmat([20;35]*pi/180,[1 N]);  % constant, multi-band
flips(1:2,1:N,3) = repmat(vfa_const_amp(N, pi/2), [2,1]);  % RF compensated variable flip angle
flips(1:2,1:N,4) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective variable flip angle
		    vfa_const_amp(N, pi/2, exp(-TR * ( - k12)))];
flips(1:2,1:N,5) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % max lactate SNR variable flip angle
		    vfa_opt_signal(N, exp(-TR * ( R1(2))))];
flips(1:2,1:N,6) = [vfa_const_amp(N, pi/2, exp(-TR * (k12))); ... % saturation recovery
		    ones(1,N)*pi/2];

         t = [0:N-1]*TR + Tin;

        figure
 subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
 subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('Lactate flips')
legend('constant','multiband', 'vfa', 'T1-effective vfa', 'max lactate SNR vfa', 'Saturation Recovery')
        
%% Test fitting

for Iflips = 1:size(flips,3)
     [Mxy Mz] = simulate_2site_model(Mz0, R1, [KPL 0], flips(:,:,Iflips), TR, input_function);
     % no noise
    [KPLfit(Iflips) Sfit(1:size(Mxy,2),  Iflips)] = fit_kPL(Mxy, TR, flips(:,:,Iflips), R1);
    % add noise
    Sn(1:size(Mxy,1), 1:size(Mxy,2),  Iflips) = Mxy + randn(size(Mxy))*std_noise;
    [KPLfitn_complex(Iflips) Snfit_complex(1:size(Mxy,2),  Iflips)] = fit_kPL(Sn(:,:,Iflips), TR, flips(:,:,Iflips),R1);
    % magnitude fitting with noise
    [KPLfitn_mag(Iflips) Snfit_mag(1:size(Mxy,2),  Iflips)] = fit_kPL(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),R1,[],std_noise/10);
 end

disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1(1), R1(2), KPL))
disp('Noiseless fit results:')
disp(['KPLfit = ' num2str(KPLfit)])
disp('Noisy complex fit results:')
disp(['KPLfitn_complex = ' num2str(KPLfitn_complex)])
disp('Noisy magnitude fit results:')
disp(['KPLfitn_mag = ' num2str(KPLfitn_mag)])
 
  figure
 subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
 subplot(122) , plot(t, squeeze(Sn(2,:,:)))
    hold on, plot(t, squeeze(Snfit_complex(:,:)),':')
    plot(t, squeeze(Snfit_mag(:,:)),'--')
title('Lactate signals and fits (dashed)')
legend('constant','multiband', 'vfa', 'T1-effective vfa', 'max lactate SNR vfa', 'Saturation Recovery')
