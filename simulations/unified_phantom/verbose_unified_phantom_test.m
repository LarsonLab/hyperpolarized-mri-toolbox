clear; close all;
mask = load('util/mask.mat').masks;

create_graphs = false;

% tissue_structure
heart = tissue_structure("heart", mask, ["lv" "rv" "lvmy" "rvmy"]);

% GLOBAL PK PARAMS
TR = 3.6;
Nt = 30;
flips = repmat([20; 30; 30;], 1, Nt) .* (pi/180);

% LV PK PARAMS ------------------
disp("LV stuff")
% - metabolites
pyr = metabolite( ...
    Mz0=0, ...
    R1=1/30, ...
    k=[0,0]);

lac = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=0.0075);

bic = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=0.0011);

% - input function
Tarrival = 7;
Tbolus = 10;
input_function = realistic_input_function(Nt, TR, Tarrival, Tbolus);

% - put it all into pk_params
lv_pk_params = pk_params(...
    substrate = pyr, ...
    products=[lac, bic], ...
    TR=TR, ...
    input_function=input_function, ...
    flips=flips);

% simulate lv metabolite dynamics
lv_met_dynamics = pk_model.generate_met_dynamics(lv_pk_params);

% plotting
if create_graphs
    tpts = 1:30;
    figure; plot(tpts, lv_met_dynamics(1,:))
    hold on
    plot(tpts, lv_met_dynamics(2,:))
    plot(tpts, lv_met_dynamics(3,:))
    hold off
end




% RV PK PARAMS ----------------
disp("RV stuff")
% metabolites
pyr = metabolite( ...
    Mz0=1, ...
    R1=1/30, ...
    k=[0,0]);

lac = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=[0.0045,0]);

bic = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=[0.0005,0]);

% - input function
Tarrival = 0;
Tbolus = 10;
input_function = realistic_input_function(Nt, TR, Tarrival, Tbolus);

% - put it all into pk_params
rv_pk_params = pk_params(...
    substrate = pyr, ...
    products=[lac, bic], ...
    TR=TR, ...
    input_function=input_function, ...
    flips=flips);

% simulate lv metabolite dynamics
rv_met_dynamics = pk_model.generate_met_dynamics(rv_pk_params);

% plotting
if create_graphs
    tpts = 1:30;
    figure; plot(tpts, rv_met_dynamics(1,:))
    hold on
    plot(tpts, rv_met_dynamics(2,:))
    plot(tpts, rv_met_dynamics(3,:))
    hold off
end




% LVMY PK PARAMS ----------------
disp("LVMY stuff")
% metabolites
pyr = metabolite( ...
    Mz0=0, ...
    R1=1/30, ...
    k=[0,0]);

lac = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=[0.06,0]);

bic = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=[0.04,0]);

% - input function
Tarrival = 10;
Tbolus = 10;
input_function = realistic_input_function(Nt, TR, Tarrival, Tbolus);


% - put it all into pk_params
lvmy_pk_params = pk_params(...
    substrate = pyr, ...
    products=[lac, bic], ...
    TR=TR, ...
    input_function=input_function, ...
    flips=flips);

% simulate lv metabolite dynamics
lvmy_met_dynamics = pk_model.generate_met_dynamics(lvmy_pk_params);

% plotting
if create_graphs
    tpts = 1:30;
    figure; plot(tpts, lvmy_met_dynamics(1,:))
    hold on
    plot(tpts, lvmy_met_dynamics(2,:))
    plot(tpts, lvmy_met_dynamics(3,:))
    hold off
end




% RVMY PK PARAMS ----------------
disp("LVMY stuff")
% metabolites
pyr = metabolite( ...
    Mz0=0, ...
    R1=1/30, ...
    k=[0,0]);

lac = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=[0.02,0]);

bic = metabolite( ...
    Mz0=0, ...
    R1=1/25, ...
    k=[0.01,0]);

% - input function
Tarrival = 14;
Tbolus = 10;
input_function = realistic_input_function(Nt, TR, Tarrival, Tbolus);

% - put it all into pk_params
rvmy_pk_params = pk_params(...
    substrate = pyr, ...
    products=[lac, bic], ...
    TR=TR, ...
    input_function=input_function, ...
    flips=flips);

% simulate lv metabolite dynamics
rvmy_met_dynamics = pk_model.generate_met_dynamics(rvmy_pk_params);

% plotting
if create_graphs
    tpts = 1:30;
    figure; plot(tpts, rvmy_met_dynamics(1,:))
    hold on
    plot(tpts, rvmy_met_dynamics(2,:))
    plot(tpts, rvmy_met_dynamics(3,:))
    hold off
end


% PK MODEL ------------------
disp("pk model time")
met_dynamics = cat(3, lv_met_dynamics, rv_met_dynamics, lvmy_met_dynamics, rvmy_met_dynamics);
met_dynamics = permute(met_dynamics, [3,1,2]); % = (tissue, met, time_pt)

met_images = pk_model.generate_met_images(heart, met_dynamics);

% MRI SYSTEM ----------------
disp("mri time")
sample_size = [32 32 11; 16 16 11; 24 24 11];

augmentation_params = struct(...
    "XTranslation", [-1,1], ...
    "YTranslation", [-1,1], ...
    "Scale", [0.95,1.1], ...
    "XReflection", true, ...
    "Rotation", [-5,5]);
augmentation_params = namedargs2cell(augmentation_params);
coil_lim = [0.4, 1.2];
SNR = [150 40 20];
output_size = [32 32 11; 32 32 11; 32 32 11];

met_images = mri_system.augment(met_images, augmentation_params{:});
met_images = mri_system.apply_coil_lim(met_images, coil_lim, heart.Mask);
met_images_mres = mri_system.make_met_images_multires(met_images, sample_size);
met_images_mres = mri_system.add_rician_noise(met_images_mres, SNR);
met_images_mres = mri_system.upsample_to_output_size(met_images_mres, output_size);


%% DISPLAY 
figure;
imagescn(met_images_mres{1}(:,:,5,:), [0, max(met_images_mres{1}(:,:,5,:), [], 'all')])

figure;
imagescn(met_images_mres{2}(:,:,5,:), [0, max(met_images_mres{2}(:,:,5,:), [], 'all')])

figure;
imagescn(met_images_mres{3}(:,:,5,:), [0, max(met_images_mres{3}(:,:,5,:), [], 'all')])
