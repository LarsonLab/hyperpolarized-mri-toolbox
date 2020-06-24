% Script for testing fit_k_aKG_2HG kinetic model fitting functions

clear all

% Test values
Tin = 0; Tacq = 36; TR = 2; N = Tacq/TR;
T1_factor_invivo = .6;  % estimated ratio between in vitro and in vivo T1 
R1_aKG_C1 = 1/(52*T1_factor_invivo); R1_aKG_C5 = 1/(41*T1_factor_invivo); R1_2HG = 1/(26*T1_factor_invivo);  %T1 values from Chaumeil et al
k_aKG_2HG = 0.002;

% Should add Glutamate, update model (should it really be exchange between C1/C5???  or just ensure they start with fixed ratio and have fixed ratio in input?)

ratio_aKG_C1toC5 = 99./1.109;  % 13C enriched to 99% vs natural abundance
k_aKG_C1toC5 = 0;  %  this is a complete guess... do we expect much exchange between C1/C5?
k_aKG_C5toC1 = k_aKG_C1toC5 *ratio_aKG_C1toC5;   %

k_all = [k_aKG_C1toC5 k_aKG_C5toC1 % aKG C1 to C5 rates
    k_aKG_2HG 0]; % aKG rates to 2HG, 2HG rates to aKG

std_noise = 0.001;
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
        Tbolus = 6;  Tarrival = 0;
        Ibolus = [1:round(Tbolus/TR)] + round(Tarrival/TR);
        Rinj = 1/Tbolus; % normalize for a total magnetization input = 1
        input_function = zeros(1,N);
        Mz0 = zeros(1,Nmets);  input_function(Ibolus) =  Rinj*TR;
end

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = ones(2,N)*20*pi/180;  flip_descripton{1} = 'constant, single-band';
flips(1:2,1:N,2) = repmat([3;35]*pi/180,[1 N]);  flip_descripton{2} = 'constant, multi-band';

k12 = 0.05; % for variable flip angle designs
flips(1:2,1:N,3) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective pyruvate variable flip angles
    vfa_opt_signal(N, exp(-TR * ( R1_2HG)))]; % max lactate SNR variable flip angle
flip_descripton{3} = 'max product SNR variable flip, multi-band';

delay_start = 0;
if delay_start
    Tdelay = 10;
    Ndelay = round(Tdelay/TR);
    flips(:,1:Ndelay,:) = eps;
end

% dimension 1: aKG_C1, aKG_C5, 2HG
% put same flips for C5 and 2HG
flips = cat(1, flips, flips(2,:,:));  % duplicate

% input includes both C1 and C5 label - now in simulate function
%input_function = [input_function; input_function/ratio_aKG_C1toC5 ];

N_flip_schemes = size(flips,3);

t = [0:N-1]*TR + Tin;

figure
subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('aKG C1 flips')
subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('aKG C5 and 2HG flips')
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
        ylabel('M_{XY}'), legend('aKG C1', 'aKG C5', '2-HG')
        title(flip_descripton{Iflips})
        subplot(312)
        plot(t, Mxy_overlapped(:,:,Iflips))
        ylabel('Signal'), legend('aKG C1', 'aKG C5 + 2-HG')
        subplot(313)
        plot(t, Sn(:,:,Iflips))
        ylabel('Signal with noise'), legend('aKG C1', 'aKG C5 + 2-HG')
        xlabel('time (s)')
        pause
    end
end


if delay_start
    Mxy_overlapped = Mxy_overlapped(:,Ndelay+1:end,:);
    Sn = Sn(:,Ndelay+1:end,:);
    flips = flips(:,Ndelay+1:end,:);
    t= t(Ndelay+1:end);
    N = N-Ndelay;
end



%%
% initial parameter guesses
R1_aKG_C1_est = R1_aKG_C1; R1_aKG_C5_est = R1_aKG_C5; R1_2HG_est = R1_2HG;  % same used in simulation
k_aKG_2HG_est = 0.0 ;

k_aKG_C1toC5_est = k_aKG_C1toC5; % same used in simulation
k_aKG_C5toC1_est = k_aKG_C5toC1;

plot_fits = 1;
fit_function = @fit_aKG_kinetics;

%% Test fitting
disp('2-site model: aKG -> 2HG')
disp('Fitting k_aKG_2HG with fixed relaxation rates and fixed C1/C5 ratio:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1_aKG_C1 = R1_aKG_C1_est;
params_fixed.R1_aKG_C5 = R1_aKG_C5_est;
params_fixed.R1_2HG = R1_2HG_est;
% known ratio/rates for C1/C5 exchange
params_fixed.k_aKG_C1toC5 = k_aKG_C1toC5_est;
params_fixed.k_aKG_C5toC1 = k_aKG_C5toC1_est;

% just fit aKG 2HG rate
params_est.k_aKG_2HG = k_aKG_2HG_est;

for Iflips = 1:N_flip_schemes
    % no noise
    if delay_start
        params_est.S0_aKG_C1 = Mxy_overlapped(1,1,Iflips)./sin(flips(1,1,Iflips));
        params_est.S0_2HG =  Mxy_overlapped(2,1,Iflips)./sin(flips(2,1,Iflips)) - params_est.S0_aKG_C1/ratio_aKG_C1toC5;
    end
    [params_fit(:,Iflips) Sfit(1,1:N,  Iflips)] = fit_function(Mxy_overlapped(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    if delay_start
        params_est.S0_aKG_C1 = Sn(1,1,Iflips)./sin(flips(1,1,Iflips));
        params_est.S0_2HG =  Sn(2,1,Iflips)./sin(flips(2,1,Iflips)) - params_est.S0_aKG_C1/ratio_aKG_C1toC5;
    end
    [params_fitn_complex(:,Iflips) Snfit_complex(1,1:N,  Iflips)] = fit_function(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
end

disp('Input:')
disp(['k_aKG_2HG  '])
disp(num2str(k_aKG_2HG,2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([1;1 + Nparams_fit;1 + 2*Nparams_fit]); 
disp(['k_aKG_2HG   '])
disp([num2str(k_fit.',2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([1;1 + Nparams_fit;1 + 2*Nparams_fit]); 
disp(['k_aKG_2HG   '])
disp([num2str(k_fit.',2), flip_description_array])

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('aKG C1 signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':'), hold off
title('2HG+aKG C5 signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)

if delay_start
    % not compatible with glutamate simulation below
    return
end

disp('Press any key to continue')
disp(' ')

pause

%% add in glutamate

R1_Glu = 1/20;
k_aKG_Glu = 0.001;
k_all = [k_aKG_C1toC5 k_aKG_C5toC1 % aKG C1 to C5 rates
    k_aKG_2HG 0 % aKG rates to 2HG, 2HG rates to aKG
    k_aKG_Glu 0]; % aKG rates to Glu, Glu rates to aKG


% dimension 1: aKG_C1, aKG_C5, 2HG, Glu
% put same flips for 2HG and Glu (could change this
flips = cat(1, flips, flips(2,:,:));  % duplicate
Nmets = 4;

% generate simulated data
plot_simulated_data = 1;
noise_S = randn([Nmets-1 N])*std_noise;  % same noise for all flip schedules
Mz0 = zeros(1,Nmets);

for Iflips = 1:N_flip_schemes
    [Mxy(1:Nmets, 1:N, Iflips), Mz] = simulate_aKG_model(Mz0, [R1_aKG_C1 R1_aKG_C5, R1_2HG, R1_Glu], k_all , flips(:,:,Iflips), TR, input_function);

    % data has overlap of aKG C5 and 2HG resonances
    Mxy_overlapped(1:Nmets-1, 1:N, Iflips) = [Mxy(1, 1:N, Iflips); ...
        Mxy(2, 1:N, Iflips) + Mxy(3, 1:N, Iflips); Mxy(4, 1:N, Iflips)];
    % add noise
    Sn(1:Nmets-1, 1:size(Mxy,2),  Iflips) = Mxy_overlapped(:,:,Iflips) + noise_S;

    if plot_simulated_data
        figure(99)
        subplot(311)
        plot(t, Mxy(:,:,Iflips))
        ylabel('M_{XY}'), legend('aKG C1', 'aKG C5', '2-HG', 'Glu')
        title(flip_descripton{Iflips})
        subplot(312)
        plot(t, Mxy_overlapped(:,:,Iflips))
        ylabel('Signal'), legend('aKG C1', 'aKG C5 + 2-HG', 'Glu')
        subplot(313)
        plot(t, Sn(:,:,Iflips))
        ylabel('Signal with noise'), legend('aKG C1', 'aKG C5 + 2-HG', 'Glu')
        xlabel('time (s)')
        pause
    end
end

% initial parameter guesses
R1_aKG_C1_est = R1_aKG_C1; R1_aKG_C5_est = R1_aKG_C5; R1_2HG_est = R1_2HG;  R1_Glu_est = R1_Glu; % same used in simulation
k_aKG_2HG_est = 0.0 ;
k_aKG_Glu_est = 0.0 ;

k_aKG_C1toC5_est = k_aKG_C1toC5; % same used in simulation
k_aKG_C5toC1_est = k_aKG_C5toC1;

plot_fits = 1;
fit_function = @fit_aKG_kinetics;

% Test fitting
disp('3-site model: aKG -> 2HG, aKG -> Glutamate')
disp('Fitting k_aKG_2HG and k_aKG_Glu with fixed relaxation rates and fixed C1/C5 ratio:')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1_aKG_C1 = R1_aKG_C1_est;
params_fixed.R1_aKG_C5 = R1_aKG_C5_est;
params_fixed.R1_2HG = R1_2HG_est;
params_fixed.R1_Glu = R1_Glu_est;
% known ratio/rates for C1/C5 exchange
params_fixed.k_aKG_C1toC5 = k_aKG_C1toC5_est;
params_fixed.k_aKG_C5toC1 = k_aKG_C5toC1_est;

% just fit aKG 2HG rate
params_est.k_aKG_2HG = k_aKG_2HG_est;
params_est.k_aKG_Glu = k_aKG_Glu_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:2,1:N,  Iflips)] = fit_function(Mxy_overlapped(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:2,1:N,  Iflips)] = fit_function(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
end

disp('Input:')
disp(['k_aKG_2HG   k_aKG_Glu '])
disp(num2str([k_aKG_2HG k_aKG_Glu],2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['k_aKG_2HG   k_aKG_Glu '])
disp([num2str(k_fit,2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:2];[1:2] + Nparams_fit;[1:2] + 2*Nparams_fit]); 
disp(['k_aKG_2HG   k_aKG_Glu '])
disp([num2str(k_fit,2), flip_description_array])

figure
subplot(131) , plot(t, squeeze(Sn(1,:,:)))
title('aKG C1 signals')
subplot(132) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':'), hold off
title('2HG+aKG C5, ')
subplot(133) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':'), hold off
title('Glu signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)


return 
%% Test fitting
disp('2-site model: aKG -> 2HG')
disp('Fitting k_aKG_2HG with fixed relaxation rates and but fitting C1/C5 exchange rate:')
disp('C1/C5 ratio should be fixed (natural abundance C13)')
disp('')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1_aKG_C1 = R1_aKG_C1_est;
params_fixed.R1_aKG_C5 = R1_aKG_C5_est;
params_fixed.R1_2HG = R1_2HG_est;

% estimate ratio/rates for C1/C5 exchange

% perturb slightly from actual values to look for convergence
k_aKG_C1toC5_est = .5;
k_aKG_C5toC1_est = k_aKG_C1toC5_est *ratio_aKG_C1toC5;  
% VERY NICE seems to get this correctly estimated without much issue!

params_est.k_aKG_C1toC5 = k_aKG_C1toC5_est;
params_est.k_aKG_C5toC1 = k_aKG_C5toC1_est;

% just fit aKG 2HG rate
params_est.k_aKG_2HG = k_aKG_2HG_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1,1:N,  Iflips)] = fit_function(Mxy_overlapped(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1,1:N,  Iflips)] = fit_function(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
end

disp('Input:')
disp(['k_aKG_2HG  '])
disp(num2str(k_aKG_2HG,2))

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
Nparams_fit = length(k_fit)/3;
k_fit = k_fit([1;1 + Nparams_fit;1 + 2*Nparams_fit]); 
disp(['k_aKG_2HG   '])
disp([num2str(k_fit.',2), flip_description_array])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([1;1 + Nparams_fit;1 + 2*Nparams_fit]); 
disp(['k_aKG_2HG   '])
disp([num2str(k_fit.',2), flip_description_array])

figure
subplot(121) , plot(t, squeeze(Sn(1,:,:)))
title('aKG C1 signals')
subplot(122) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':'), hold off
title('2HG+aKG C5 signals and fits (dots=complex fit, dashed=magnitude)')
legend(flip_descripton)




