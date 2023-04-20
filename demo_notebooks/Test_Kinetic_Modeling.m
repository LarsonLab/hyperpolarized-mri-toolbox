%This toolbox contains several kinetic modeling and quantification methods, 
%as well as simulation tools and phantoms, which are demonstrated in this notebook.

%SOMEWHAT WORKS, LINES 213-225 DOES NOT WORK
%RECOMMED DOING IT ONE-BY-ONE
fit_function =@fit_pyr_kinetics;
plot_fits = 1;

%Test values 
Tin = 0;
Tacq = 48;
TR = 3;
N = Tacq/TR;
R1P = 1/25;
R1L = 1/25;
kPL = 0.05; 
Mz0 = [0,0];
std_noise = 0.01;

%setup input function
input_function = zeros(1,N);
input_condition = 1;
switch input_condition
    case 1 
     disp('gamma variate input function - most realistic')
     Tarrival = 0;
     Tbolus = 12;
     input_function = realistic_input_function(N,TR,Tarrival,Tbolus);
    case 2
      disp('boxcar input function')
      Tboulus = 12;
      Tarrival = 0;
      Ibolus = (1:round(Tboulus/TR))+round(Tarrival/TR);
      Rinj = 1/Tboulus;
      Mz0 = [0,0];
      input_function(Ibolus) = Rinj*TR;
    case 3
        disp('no input function')
        Mz0 = [1.5,0];
    case 4
        disp('no input function, delayed start')
end

t = (0:N-1)*TR +Tin;

figure
plot(t,input_function)
xlabel('time(s)'), ylabel('input function ')

% Choose a flip angle scheme

flip_scheme = 1;

switch flip_scheme 
    case 1
        flips = ones(2,N)*30*pi/180;
        flip_description = 'constant, single-band';
    case 2
        flips = repmat([20;35]*pi/180,[1 N]);
        flip_description = 'constant, multi-band';
    case 3 
        k12 = 0.05;
        flips = [vfa_const_amp(N,pi/2,exp(TR(k12)));...
            vfa_opt_signal(N, pi/2, exp(TR*(R1L)))];
        flip_description = 'T1-effective variable flip';
    case 4
        flips = [vfa_const_amp(N, pi/2, exp(-TR *(k12)));...
            vfa_const_amp(N,pi/2,exp(-TR*(-k12)))];
        flip_description = 'saturation recovery';
    case 5
        flips = [vfa_const_amp(N,pi/2,exp(-TR*(k12)));...
            ones(1,N)*pi/2];
        flip_description = 'saturation recovery';
end

% disp(flip_description)
% 
% subplot(121),plot(t,squeeze(flips(1,:,:))*180/pi)
% title('Pyruvate Flips')
% xlabel('time(s)'),ylabel('degrees')
% subplot(122),plot(t,squeeze(flips(2,:,:))*180/pi)
% title('Lactate flips')
% xlabel('time(s)'),ylabel('degrees')
% 
%generate simulated data
noise_S = randn([2 N])*std_noise; %same noise for all flip schedules
[Mxy, Mz] = simulate_Nsite_model(Mz0,[R1P R1L], [kPL 0], flips, TR, input_function);
%add noise
Sn = Mxy + noise_S;
% 
% subplot(211),plot(t, Mz)
% title('M_Z')
% 
% subplot(212),plot(t,Mxy)
% title('M_{XY} (signal)')
% legend('pyruvate','lactate')
% 
% figure 
% plot(t,Sn/std_noise)
% legend('pyruvate','lactate')
% title('signal with noise')
% xlabel('time')
% ylabel('(SNR)')

% %variable estimates 
R1P_est = 1/25;
R1L_est = 1/25;
kPL_est = .02;

%setup params_fixed and params_est structures required by the fitting
%functions 
params_fixed.R1P = R1P_est;
params_fixed.R1L = R1L_est;

%estimated kPL, and use this value as initial guess
params_est.kPL = kPL_est;

%noiseless
[params_fit, Sfit] = fit_function(Mxy, TR, flips, params_fixed, params_est,[],plot_fits);
disp('no noise')

%added noise, real-valued(not magnitude) data
[params_fit_noise ,Sfit_noise] = fit_function(Sn,TR, flips, params_fixed, params_est, [], plot_fits);
disp('added noise, real-valued data')


% added noise, magnitude data
% NOTE that this requires an estimate of the noise
[params_fit_noise_mag, Sfit_noise_mag] = fit_function(abs(Sn), TR, flips, params_fixed, params_est, std_noise, plot_fits);
disp('added noise, magnitude data')
% 
% %parameter estimates 
% R1P_est = 1/25;
% R1L_est = 1/25;
% kPL_est = .02;
% 
clear params_fixed params_est
% setup params_fixed and params_est structures required by the fitting functions
% fix pyruvate relaxation rates
params_fixed.R1P = R1P_est;
% estimate kPL and t1 of lactate
params_est.kPL = kPL_est; 
params_est.R1L = R1L_est;
% set constraints on lactate T1:
params_est.R1L_lb = 1/40;
params_est.R1L_ub = 1/15;
% 
% % [params_fit,Sfit] = fit_function(Mxy,TR,flips,params_fixed,params_est, [],plot_fits);
% % disp('no noise')
% % 
% % 
% % % added noise, real-valued (not magnitude) data
% 
% % [params_fit_noise, Sfit_noise] = fit_function(Sn, TR, flips, params_fixed, params_est, [], plot_fits);
% % disp('added noise, real-valued data')
% 
% 
% % % added noise, magnitude data
% % % NOTE that this requires an estimate of the noise
% %Test Values
% Tin = 0;
% Tacq = 48;
% TR = 3;
% N = Tacq/TR;
% R1P = 1/25;
% R1L = 1/25;
% R1B = 1/15;
% R1A = 1/25;
% kPL = 0.05;
% kPB = 0.03;
% kPA = 0.03;
% std_noise = 0.005;

% generate simulated data
noise_S = randn([2 N])*std_noise;  % same noise for all flip schedules
[Mxy,Mz] = simulate_Nsite_model(Mz0, [R1P R1L], [kPL 0], flips, TR, input_function);
Sn = Mxy + noise_S;

[params_fit_noise_mag, Sfit_noise_mag] = fit_function(abs(Sn), TR, flips, params_fixed, params_est, std_noise, plot_fits);
disp('added noise, magnitude data')

%update flip angles(apply same to all products)
flips_4site = cat(1, flips, repmat(flips(2,:),[2 1 1]));
%update initial magnetization 
Mz0_4site = [Mz0,Mz0(2),Mz0(2)];

%genereate simulated data
noise_S = randn([4 N]*ceil(std_noise));
% [Mxy,Mz] = simulate_Nsite_model(Mz0_4site, [R1P R1L R1B R1A],[kPL 0; kPB 0; kPA 0], ...
%     flips_4site,TR,input_function);
% add noise
%Sn = Mxy + noise_S;

% subplot(311),plot(t,Mz)
% title('M_Z')
% subplot(312),plot(t,Mxy)
% title('M_{XY} (signal)')
% legend('pyruvate','lactate','bicarb','alanine')
% subplot(313),plot(t,Sn)
% title('signal with noise')

% parameter estimates
R1P_est = 1/25; R1L_est = 1/25; R1B_est = 1/15; R1A_est = 1/25;
kPL_est = .02; kPB_est = .02; kPA_est = .02;

% setup params_fixed and params_est structures required by the fitting functions
clear params_fixed params_est 
% fix relaxation rates
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
% estimated kPX, and use this value as initial guess
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;
R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;

% % noiseless
% [params_fit Sfit] = fit_function(Mxy, TR, flips_4site, params_fixed, params_est, [], plot_fits);
% disp('no noise')

% added noise, real-valued (not magnitude) data
% [params_fit_noise, Sfit_noise] = fit_function(Sn, TR, flips_4site, params_fixed, params_est, [], plot_fits);
% disp('added noise, real-valued data')

% added noise, magnitude data
% NOTE that this requires an estimate of the noise

% [params_fit_noise_mag, Sfit_noise_mag] = fit_function(abs(Sn), TR, flips_4site, params_fixed, params_est, std_noise, plot_fits);
% disp('added noise, magnitude data')