% Compare kinetic models


% Test values
Tin = 0; Tacq = 48; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; kPL = 0.05; std_noise = 0.01;

        % gamma variate input function - most realistic
        Tarrival = 0;  Tbolus = 12;
        input_function = realistic_input_function(N, TR, Tarrival, Tbolus);
        Mz0 = [0,0]; 


clear params_est params_fixed params_fit params_fit_noise
%params_fit(1,1) = struct([]);
%params_fit_noise(1,1) = struct([]);
        
        % default fitting parameters
        R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;
        params_fixed.R1P = R1P_est;             
        params_est.kPL = kPL_est;params_fixed.R1L = R1L_est;
%        params_est.kPL = kPL_est;params_est.R1L = R1L_est;
        
        % NOTES:
        % - similar for just fitting kPL (10^-7 diffs!)
        % - with T1 too (10^-6 diffs)
fit_fcn = {@fit_kPL, @fit_pyr_kinetics};

flips = ones(2,N)*30*pi/180;  flip_descripton = 'constant, single-band';

kPLs = [0:0.001:0.1];
for IkPL = 1:length(kPLs)
    kPL = kPLs(IkPL);

noise_S = randn([2 N])*std_noise;  % same noise for all flip schedules
    [Mxy, Mz] = simulate_Nsite_model(Mz0, [R1P R1L], [kPL 0], flips, TR, input_function);
    % add noise
    Sn = Mxy + noise_S;

    
    for n = 1:length(fit_fcn)
        params_fit = fit_fcn{n}(Mxy, TR, flips, params_fixed, params_est);
        kPL_fit(IkPL,n) = params_fit.kPL;
        
        params_fit_noise = fit_fcn{n}(Sn, TR, flips, params_fixed, params_est);
        kPL_fit_noise(IkPL,n) = params_fit_noise.kPL;
        
        
    end
    
end

%%
figure(1)
subplot(121)
plot(kPL_fit(:,1), kPL_fit(:,2),'.')

subplot(122)
plot(kPL_fit_noise(:,1), kPL_fit_noise(:,2),'.')


figure(2)
subplot(121)
plot(kPL_fit(:,1), kPL_fit(:,1) - kPL_fit(:,2),'.')

subplot(122)
plot(kPL_fit_noise(:,1), kPL_fit_noise(:,1) - kPL_fit_noise(:,2),'.')
