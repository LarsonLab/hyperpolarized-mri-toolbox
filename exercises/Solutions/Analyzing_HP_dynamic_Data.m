% setup hyperpolarized-mri-toolbox

% IT WORKS

%% Load an imaging dataset
% sample_data available in toolbox
% Easiest: rat or mouse EPI data (single channel RF coil)
% Challenge: Human Brain EPI (multi-channel coil, need to implement coil combination)
clear all;
load("sample_data/Rat Kidneys EPI/exp2_constant/exp2_constant.mat")


% lists all variables
whos

%% Visualize dataset
% for all, find areas that have signal

% 1. display an image from each metabolite at a single time-point

%makes a new figure window, and return its handle.
figure()
%this creates a chart 
subplot(1,2,1)
%it shows 3d arry of pyr,between the range of pyr vector 
imshow(pyr(:,:,8),[min(pyr(:)) max(pyr(:))])
title('pyruvate 7th Time Point')
subplot(122)
imshow(lac(:,:,8),[min(lac(:)) max(lac(:))])
title('Lactate 7th Time Point')
%2. display all slices and/or timepoints for a single metabolite
%permute rearranges the dimensions of pyr/lac and orders is specied by the
%vector
% montage() command is useful, but expects a MxNx1xP matrix for grayscale images
% so use something like this to add a singleton dimension to data:
% data_for_montage = permute(data, [1 2 ndims(data)+1 plot_dimension])
pyr_montage = permute(pyr,[1,2,4,3]);
lac_montage = permute(lac,[1,2,4,3]);

figure()
montage(pyr_montage,'DisplayRange',[min(pyr(:)) max(pyr(:))])
title('Pyrucvate across time points')

figure()
montage(lac_montage,'DisplayRange',[min(lac(:)) max(lac(:))])
title('Lactate across time points')


 % 3. extract the time curves for all metabolites from a single voxel, and
 % plot versus time
voxel = [21,17]; 

%squeeze can remove the demension of length 1
pyr_time = squeeze(pyr(voxel(1),voxel(2),:));
lac_time = squeeze(lac(voxel(1),voxel(2),:));
Number_of_times = size(pyr,3); 
time = 0:Number_of_times-1;

figure()
plot(time,pyr_time,time,lac_time,'LineWidth',2)
legend('Pyruvate','Lactate')
xlabel('time points')
ylabel('|Mxy|')




 %% Excersise:Model-Free Metrics 

 %Plot Lactate/Pyruvate vs time for a single voxel 
voxel = [21,17]; 

%squeeze can remove the demension of length 1
pyr_time = squeeze(pyr(voxel(1),voxel(2),:));
lac_time = squeeze(lac(voxel(1),voxel(2),:));
Number_of_times = size(pyr,3); 
time = 0:Number_of_times-1;

figure()
plot(time,pyr_time,time,lac_time,'LineWidth',2)
legend('Pyruvate','Lactate')
xlabel('time points')
ylabel('|Mxy|')
title('Lactate/Pyruvate')




% Display an image of the AUCratio for each slice
% use compute_AUCratio()
% create a signal mask to only show relevant voxel

%this concatenates lac to the pyr along dimension 3 
S = cat(3,pyr_montage,lac_montage);

%mask out vozels that do not have sufficient pyruvate SNR 
nx = size(pyr,1); 
ny = size(pyr,2);
pyr_AUC = sum(pyr,3);
pyr_max = max(pyr_AUC(:));
ind = find(pyr_AUC >= 0.15*pyr_max);
s_flat = reshape(S,[nx*ny,size(S,3),size(S,4)]);
s_masked = s_flat(ind,:,:);

AUCratio = compute_AUCratio(s_masked);

AUC_map = zeros([nx*ny,1]);
AUC_map(ind) = AUCratio;
AUC_map = reshape(AUC_map,[nx,ny]);

figure()
imagesc(AUC_map)
colormap('hot')
colorbar
title('Lac/Pyr AUC Ratio')

%% Excersise:Kinetic Modeling

%Compute kPL and display kPL maps
%use fit_pyr_kinetics, with fixed T1 vaules
%create a signal mask to only fit and show relevant voxels 

%estimated values for relaxation and rate constants 
R1P_est = 1/20; 
R1L_est = 1/25;
kPL_est = 0.002;

%fixed parameters
params_fixed.R1P = R1P_est;
params_fixed.R1L = R1L_est;

%fit parameters
params_est.kPL = kPL_est;

flips =[flips_pyr'; flips_lac']; %restructure given flip angles 
%mask out voxels that do not have sufficient pyruvvate SNR
nx = size(pyr,1); 
ny = size(pyr,2);
pyr_AUC = sum(pyr,3);
pyr_max = max(pyr_AUC(:));
ind = find(pyr_AUC >= 0.15*pyr_max);
s_flat = reshape(S,[nx*ny,size(S,3),size(S,4)]);
s_flat2 = double(s_flat);
s_masked = s_flat2(ind,:,:);

[params_fit, Sfit, ufit, error_metrics] = fit_pyr_kinetics (s_masked, TR, flips, params_fixed, params_est);
%Construct kPL map
kpl_map = zeros([nx*ny,1]);
kpl_map(ind) = params_fit.kPL;
kpl_map = reshape(kpl_map,[nx ny]);
figure()
imagesc(kpl_map)
colormap('parula')
colorbar

%Compute kPL with different fixed T1 values 
clear params_est
clear params_fixed

%estimated values for relaxation and rate constants
R1P_est = 1 / 25; 
R1L_est = 1 / 30; 
kPL_est = 0.002; 

%fixed parameters
params_fixed.R1P = R1P_est;
params_fixed.R1L = R1L_est;

%fit parameters
params_est.kPL = kPL_est;

[params_fit, Sfit, ufit, error_metrics] = fit_pyr_kinetics(s_masked, TR, flips, params_fixed, params_est);

%construct kPL map
kpl_map = zeros([nx*ny, 1]);
kpl_map(ind) = params_fit.kPL;
kpl_map = reshape(kpl_map, [nx ny]);
figure()
imagesc(kpl_map)
colormap('parula')
colorbar

%Compute kPL while fitting (not fixing) T1 values

clear params_est
clear params_fixed

%estimated values for relaxation and rate constants
R1P_est = 1 / 20; 
R1L_est = 1 / 25; 
kPL_est = 0.002; 

%fixed parameters
params_fixed.R1P = R1P_est; 

%fit parameters -- this time fit lactate T1 as well
params_est.kPL = kPL_est; 
params_est.R1L = R1L_est;

[params_fit, Sfit, ufit, error_metrics] = fit_pyr_kinetics(s_masked, TR, flips, params_fixed, params_est);


%construct kPL map
kpl_map = zeros([nx*ny, 1]);
kpl_map(ind) = params_fit.kPL;
kpl_map = reshape(kpl_map, [nx ny]);
figure()
imagesc(kpl_map)
colormap('parula')
colorbar

%% Excersise:Load and Visualize a HP spectroscopy dataset

%Load a spectroscopic imaging dataset 
%sample_data available in toolbox
%Options human prostate dynamic MRSI,spectral data is [f,x,y,y]

clear all;
load('/Users/ernestodiaz/Desktop/hyperpolarized-mri-toolbox/sample_data/Human Prostate Dynamic MRSI/pc6071_spectra.mat')

whos
% Visualize dataset
% for all, find areas that have signal

% 1. use plot_voxels() to show all spectra at a given time point
% may require using permute() to match dimensionality
s_new = permute(abs(spectra(:,:,:,7)),[2 3 1]);
%whos s_new
figure()
plot_voxels(s_new)

% 2. display a spectrum from a single voxel at a single time-point
single_spectrum_1 = squeeze(abs(spectra(:,4,12,7)));
%whos single_spectrum_1
single_spectrum_2 = squeeze(abs(spectra(:,4,10,7)));
%whos single_spectrum_2
figure()
subplot(211)
plot(1:size(spectra,1),single_spectrum_1,'LineWidth',2)
xlabel('Frequency')
ylabel('|Mxy')
title('Location:[4,12], 7th time point')
subplot(212)
plot(1:size(spectra,1),single_spectrum_2,'LineWidth', 2)
xlabel('frequency')
ylabel('|Mxy|')
title('Location: [4,10], 7th time point')

% Extract metabolite maps

% 1. measure a metabolite (e.g. pyruvate) peak location from a single voxel spectrum
[M,Ipyr] = max(single_spectrum_1);
[M,Ilac] = max(single_spectrum_2);

% 2. compute peak heights for all voxels and time-points
% use max() with an input of the measured peak location + [-5,5] to account for inhomogeneity
pyr_peaks = max(abs(spectra(Ipyr-5:Ipyr+5,:,:,:)), [],1);
lac_peaks = max(abs(spectra(Ilac-5:Ilac+5,:,:,:)), [],1);

% 2. (bonus) integrate peak area from phased spectra
% find_phase_corr() function maybe helpful
Ifit = [Ipyr-25:Ipyr+25, Ilac-25:Ilac+25]; 

%apply phase correction 
phase_corr = zeros(size(spectra,2), size(spectra,3));
whos phase_corr
for a=1:size(spectra,2)
    for b=1:size(spectra,3)
        spec = spectra(Ifit,a,b,:); 
        phase_corr(a,b) = find_phase_corr(spec(:)); 

        spectra(:,a,b,:) = spectra(:,a,b,:)*exp(1i*phase_corr(a,b));
    end
end

%intergrate peak using corrected spectra
pyr_area = abs(sum(spectra(Ifit(1),:,:,:),1));
lac_area = abs(sum(spectra(Ifit(2),:,:,:),1));


% Visualize Images
% 1. display a metabolite image for all time-points
pyr_peaks_montage = permute(pyr_peaks,[2 3 1 4]);
lac_peaks_montage = permute(lac_peaks,[2 3 1 4]);

figure()
montage(pyr_peaks_montage,'DisplayRange',[min(pyr_peaks(:)) max(pyr_peaks(:))]);
title('Pyruvate across time points')

figure()
montage(lac_peaks_montage,'DisplayRange',[min(lac_peaks(:)) max(lac_peaks(:))]);
title('Lactate across time points')


% 2. extract the time curves from a single voxel, and plot
voxel= [4,12];
pyr_time = squeeze(pyr_peaks(:,voxel(1),voxel(2),:));
lac_time = squeeze(lac_peaks(:,voxel(1),voxel(2),:));
time = 0:Nt-1;

figure()
plot(time,pyr_time,time,lac_time,'LineWidth',2)
legend('Pyruvate','Lactate')
xlabel('time points')
ylabel('|Mxy|')






