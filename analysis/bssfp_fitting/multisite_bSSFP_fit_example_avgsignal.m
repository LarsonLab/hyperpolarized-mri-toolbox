%% bSSFP-lactate Data
% Load in and prepare Dataset

clear; close all;

% load data
data_dir = './../../sample_data/Rat Kidneys Spiral bSSFP/';
load(strcat(data_dir,'shot1_bssfp_lac.mat'))
load(strcat(data_dir,'roi.mat'))

%mask and avg image per metabolite
mask = kidneys;
Nt = size(pyr_exp1,4);
Sp = mask_and_avg(pyr_exp1, mask);
Sl = mask_and_avg(lac_exp1, mask);
Sa = mask_and_avg(ala_exp1, mask);
S_complex = [Sp; Sl; Sa] ./ max(Sp); %normalize by max signal
S = abs(S_complex); % fit to magnitude data

% Check average signal before fitting 
figure()
plot(0:Nt-1, S)
legend('pyruvate', 'lactate', 'alanine')
xlabel('time points')
ylabel('average magnitude signal over ROI')

%% PK Model Fitting: bSSFP-lactate Data
%%%% Example: 2 site pyr->lac, pyruvate: 3D GRE, lactate: 3D bSSFP  %%%%%%

% only fit pyruvate and lactate for now
S = S(1:2, :);

% ----- set parameters -----
% flip angles
FAP = 3; % [deg]
FAL = 60; % [deg]
Nz = 16; interleaves = 4;
flips.P = repmat(FAP, [1, Nz]);
flips.L = repmat(FAL, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);
cat_flips = [4, -16, 24, -36, 48, -60]; % bSSFP catalyzation flip angles

% TR
TRP = .18; % [s]
TRL = .01529; % [s] 
TR.P = repmat(TRP, [1, Nz]);
TR.L = repmat(TRL, [1, Nz*interleaves]);
cat_TR = [TRL, TRL, TRL, TRL, TRL, TRL]; % bSSFP catalyzation TR

% relaxation and kPL
% add parameters that are fixed to params_fixed struct
% add parameters to fit to params_est struct
params_fixed.R1P = 1/30; % [1/s]
params_fixed.R1L = 1/25;
params_fixed.R2L = 1/0.8; % set lactate T2 because lactate-bSSFP
params_est.kPL = 0.002;

TempRes = 4; % [s], time between each time point
params_est.Mz0_P = 3; 
params_est.Mz0_L = 0.5;
acq_sequence = ["3DGRE", "3DbSSFP"];
verbose = 1;

% ----- run fitting -----
S = reshape(S, [1 size(S)]);
[fitparams, error, Mxy_fit, Mz_fit, input_fit] = multisite_bSSFP_fit(S, params_fixed, params_est, flips, TR, TempRes, acq_sequence,...
    'cat_flips', cat_flips, 'cat_TR', cat_TR, 'verbose', verbose);

%% GRE-all Data
% Load in and prepare Dataset

clear; close all;

% load data
data_dir = './../../sample_data/Rat Kidneys Spiral bSSFP/';
load(strcat(data_dir,'shot2_gre_all.mat'))
load(strcat(data_dir,'roi.mat'))

%mask and avg image per metabolite
mask = kidneys;
Nt = size(pyr_exp2,4);
Sp = mask_and_avg(pyr_exp2, mask);
Sl = mask_and_avg(lac_exp2, mask);
Sa = mask_and_avg(ala_exp2, mask);
S_complex = [Sp; Sl; Sa] ./ max(Sp); %normalize by max signal
S = abs(S_complex); % fit to magnitude data

% Check average signal before fitting 
figure()
plot(0:Nt-1, S)
legend('pyruvate', 'lactate', 'alanine')
xlabel('time points')
ylabel('average magnitude signal over ROI')

%% PK Model Fitting: GRE-all Data
%%%% Example: 2 site pyr->lac, pyruvate: 3D GRE, lactate: 3D GRE  %%%%%%

% only fit pyruvate and lactate for now
S = S(1:2, :);

% ----- set parameters -----
% flip angles
FAP = 3; % [deg]
FAL = 7.67; % [deg]
Nz = 16; % number of z-encodes
flips.P = repmat(FAP, [1, Nz]);
flips.L = repmat(FAL, [1, Nz]);

% TR
TRP = .125; % [s]
TRL = .125; % [s] 
TR.P = repmat(TRP, [1, Nz]);
TR.L = repmat(TRL, [1, Nz]);

% relaxation and kPL
% add parameters that are fixed to params_fixed struct
% add parameters to fit to params_est struct
params_fixed.R1P = 1/30; 
params_fixed.R1L = 1/25;
params_est.kPL = 0.002;

TempRes = 4; % [s], time between each time point
params_est.Mz0_P = 3; 
params_est.Mz0_L = 0.5; 
acq_sequence = ["3DGRE", "3DGRE"];
verbose = 1;

% ----- run fitting -----
S = reshape(S, [1 size(S)]);
[fitparams, error, Mxy_fit, Mz_fit, input_fit] = multisite_bSSFP_fit(S, params_fixed, params_est, flips, TR, TempRes, acq_sequence,...
    'verbose', verbose);

%% util functions

function S = mask_and_avg(data, mask)
% assumes data and mask are the same size
% assumes last dimension is time and want to take avg across spatial
% dimensions

    sz = size(data);
    masked = double(data) .* mask;
    flattened = reshape(masked, [prod(sz(1:end-1)) , sz(end)]);
    ind = find(flattened(:,1));
    S = mean(flattened(ind,:),1);
end
