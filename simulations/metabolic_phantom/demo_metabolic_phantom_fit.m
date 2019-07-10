
clear all

%% Setup phantom

nx = 16;
ny = 16;
nz = 1;
kTRANS_low  = 0.02;
kTRANS_high = 0.05;
kPL_low     = 0.01;
kPL_high    = 0.03;
linear_kTRANS_gradient = true;
linear_kPL_gradient = true;

[kTRANS, kPL] = metabolic_phantom(nx, ny, nz, kTRANS_low, kTRANS_high, kPL_low, kPL_high, linear_kTRANS_gradient, linear_kPL_gradient);

%% Generate dynamics

% Test values
Tin = 0; Tacq = 75; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25;
std_noise = 1e-4;

t = [1:N]*TR;
Tarrival = 0; Tbolus = 10;
input_function = realistic_input_function(N, TR, Tarrival, Tbolus);% normalize for a total magnetization input = 1
Mz0 = [0,0];

flips(1:2,1:N) = ones(2,N)*20*pi/180;  % constant, single-band

data_nonoise = zeros(nx,ny,nz,2, N);
data = zeros(nx,ny,nz,2, N);
for Ix = 1:nx
    for Iy = 1:ny
        for Iz = 1:nz
            [Mxy Mz] = simulate_Nsite_model(Mz0, [R1P R1L], [kPL(Ix,Iy,Iz) 0], flips, TR, input_function*kTRANS(Ix,Iy,Iz) );
            data_nonoise(Ix,Iy,Iz,:,:) = Mxy;
            % add noise
            noise_S = randn([2 N])*std_noise;
            data(Ix,Iy,Iz,:,:) = Mxy + noise_S;
        end
    end
end

figure('Name', 'Pyruvate data')
plot_voxels(squeeze(data(:,:,1,1,:)))

figure('Name', 'Lactate data')
plot_voxels(squeeze(data(:,:,1,2,:)))

%% Fit data

clear params_fixed params_est params_fit
R1P_est = R1P; R1L_est = R1L;  kPL_est = 0.01;
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est;

% no noise
[params_fit_nonoise Sfit_nonoise] = fit_kPL(data_nonoise, TR, flips, params_fixed, params_est);

%  noisy
[params_fit Sfit] = fit_kPL(data, TR, flips, params_fixed, params_est);

Splot = [0 kPL_high];
figure
subplot(221)
imagesc(kPL, Splot), colorbar
title('Original k_{PL}')
subplot(222)
imagesc(kTRANS), colorbar
title('Original k_{TRANS}')
subplot(223)
imagesc(params_fit_nonoise.kPL, Splot), colorbar
title('Fit k_{PL} (no noise)')
subplot(224)
imagesc(params_fit.kPL, Splot), colorbar
title('Fit k_{PL} (with noise)')
