% Script for testing fit_kPL kinetic model fitting functions

clear all

% choose fitting function to test
fit_function = @fit_kPL_perfused_voxel;
plot_fits = 1;

% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; kPL = 0.2; std_noise = 0.001;
kve = 0.05; vb = 0.1;

input_function = zeros(1,N);

input_condition = 1; % choose from various simulated starting conditions
switch input_condition
    case 1
        % gamma variate input function - most realistic
        t = [1:N]*TR;
        Tarrival = 0;  Tbolus = 12;
        input_function = realistic_input_function(N, TR, Tarrival, Tbolus);
        Mz0 = [0,0]; 
    case 2
        % boxcar input function
        Tbolus = 12;  Tarrival = 0;
        Ibolus = [1:round(Tbolus/TR)] + round(Tarrival/TR);
        Rinj = 1/Tbolus; % normalize for a total magnetization input = 1
        Mz0 = [0,0]; input_function(Ibolus) =  Rinj*TR;

end
VIF = [input_function; zeros(1,N)];
VIFscale = 1;

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = ones(2,N)*30*pi/180;  flip_descripton{1} = 'constant, single-band';
flips(1:2,1:N,2) = repmat([20;35]*pi/180,[1 N]);  flip_descripton{2} = 'constant, multi-band';
k12 = 0.05; % for variable flip angle designs
flips(1:2,1:N,3) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective pyruvate variable flip angles
    vfa_opt_signal(N, exp(-TR * ( R1L)))]; % max lactate SNR variable flip angle
flip_descripton{3} = 'max lactate SNR variable flip, multi-band';
% flips(1:2,1:N,4) = repmat(vfa_const_amp(N, pi/2), [2,1]);  % RF compensated variable flip angle
% flips(1:2,1:N,5) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective variable flip angle
%     vfa_const_amp(N, pi/2, exp(-TR * ( - k12)))];
% flips(1:2,1:N,6) = [vfa_const_amp(N, pi/2, exp(-TR * (k12))); ... % saturation recovery
%     ones(1,N)*pi/2];
N_flip_schemes = size(flips,3);

t = [0:N-1]*TR + Tin;

figure
subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('Lactate flips')
legend(flip_descripton)

flip_description_array = [repmat('    ',N_flip_schemes,1),  char(flip_descripton(:))];

% generate simulated data
noise_S = randn([2 N])*std_noise;  % same noise for all flip schedules
for Iflips = 1:N_flip_schemes
    [Mxy_ev Mz_ev Mxy_iv Mz_iv] = simulate_Nsite_perfused_voxel_model(Mz0, [R1P R1L], [kPL 0], kve, vb, flips(:,:,Iflips), TR, VIF*VIFscale);
    S(:,:,Iflips) = Mxy_ev + Mxy_iv;
    % add noise
    Sn(:, :,  Iflips) = S(:,:,Iflips) + noise_S;

    figure(Iflips)
    subplot(221)
    plot(t, Mz_iv), title('M_Z intravascular')
    subplot(222)
    plot(t, Mz_ev), title('M_Z extravascular')
    subplot(223)
    plot(t, Mxy_iv), title('M_{XY} intravascular')
    subplot(224)
    plot(t, Mxy_ev), title('M_{XY} extravascular')
    
end

% initial parameter guesses
R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;  kve_est = 0.02; vb_est = 0.1;



%% Test fitting - fit kPL only
disp('Fitting kPL, with fixed relaxation rates:')
disp('Fixing relaxation rates improves the precision of fitting, but potential')
disp('for bias in fits when incorrect relaxation rate is used')
disp(' ')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
params_fixed.kve = kve_est; params_fixed.vb = vb_est;
params_est.VIFscale = 1;
params_est.kPL = kPL_est;
% add IV function

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:2,1:size(S,2),  Iflips)] = fit_function(S(:,:,Iflips), VIF, TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:2,1:size(S,2),  Iflips)] = fit_function(Sn(:,:,Iflips), VIF, TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:2,1:size(S,2),  Iflips)] = fit_function(abs(Sn(:,:,Iflips)), VIF, TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P, R1L, kPL))
disp('Noiseless fit results:')
disp(['kPL  = ']); disp([num2str(struct2array(params_fit).'), flip_description_array])
disp('Noisy complex fit results:')
disp(['kPL  = ']); disp([num2str(struct2array(params_fitn_complex).'), flip_description_array])
disp('Noisy magnitude fit results:')
disp(['kPL  = ']); disp([num2str(struct2array(params_fitn_mag).'), flip_description_array])

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(:,:)),':')
plot(t, squeeze(Snfit_mag(:,:)),'--')
title('kPL fit: Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')

