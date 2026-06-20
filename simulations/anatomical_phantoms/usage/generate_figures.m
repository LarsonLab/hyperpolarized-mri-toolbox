% copying over brainweb
clear; close all;
addpath('../');
addpath('../brainweb_clone');

%% PARAMETERS
disp('setting up...');
t_all = tic;
tic
% tissue
mask = double(load('brainweb_fuzzy.mat').im_mask);

k_trans = [1, 0.2, 0.2;
           3, 0.4, 0.4];

% pk model params
tr = 4;
n_t = 30;
flips = repmat([20; 30; 30;], 1, n_t) .* (pi/180);
r1 = [1/30, 1/25, 1/25]; % pyr, lac, bic
k = [0, 0.05, 0.03; % lac, in order [vasc, gm, wm]
     0, 0.02, 0.01]; % bic
t_arrival = [0, 0, 0]; % vasc, gm, wm
t_bolus = 8;

% input function and mz0
n_tissues = 3;
input_function = zeros(n_tissues, n_t);
for i_tissue = 1:n_tissues
    input_function(i_tissue, :) = realistic_input_function(n_t, tr, t_arrival(i_tissue), t_bolus);
end
mz0 = [input_function(1), input_function(1)*.5, input_function(1)*.5;
       0, input_function(1)*.01, input_function(1)*.01;
       0, input_function(1)*.005, input_function(1)*.005];

% mri
coil_lim = [0.4, 1.2];
augment_params = struct(...
    "XTranslation", [-1,1], ...
    "YTranslation", [-1,1], ...
    "Scale", [0.95,1.1], ...
    "XReflection", true, ...
    "Rotation", [-5,5]);

sample_size = [16 16 8; 16 16 8; 16 16 8];
snr = [150 40 20];
output_size = [32 32 8; 32 32 8; 32 32 8];
disp(['took ', num2str(toc), 's', newline]);

%% RUNNING THE MODEL -----------------------------------------------------------
% tissue
disp('creating tissue structure...'); tic;
brain = tissue_structure("brain", mask, ["vasc", "gm", "wm"], [100,100,100]);
brain = brain.create_k_trans_map(k_trans);
disp(['took ', num2str(toc), 's', newline]);

%% pk model
disp('running pk model...'); tic;
[met_images, ~, met_dynamics, ~, met_images_no_ktrans] = pk_model.run_pk_model(mz0, r1, k, flips, tr, brain, input_function=input_function);
disp(['took ', num2str(toc), 's', newline]);

% mri
disp('running mri system...'); tic;
cell_augment_params = namedargs2cell(augment_params); % unpack the augmentation parameters
met_images_aug = mri_system.augment(met_images, cell_augment_params{:});
[met_images_coil_lim, coil_sens_weights] = mri_system.apply_coil_lim(met_images_aug, coil_lim, brain.Mask);
met_images_mres = mri_system.make_met_images_multires(met_images_coil_lim, sample_size);
met_images_mres_noise = mri_system.add_rician_noise(met_images_mres, snr);

met_images_upsampled = mri_system.upsample_to_output_size(met_images_mres_noise, output_size);
disp(['took ', num2str(toc), 's', newline]);


%% DISPLAY ---------------------------------------------------------------------
disp('plotting...'); tic;

% tissue ----------
% mask
brain.plot_alpha_composite_image(slice=50, order=[2,3,1]); % plot vasculature last
saveas(gcf, 'figures/alpha_composite.png');

%% kTRANS
slices = 10:10:100;
f = display_tiled_images(brain.K_trans_map(:,:,slices), true, 'ktrans', []);
f.Position = [0, 0, 1440, 185];
saveas(gcf, 'figures/ktrans.png');


% pk -------------
% met dynamics
% size(met_dynamics) = [tissue, met, time_pt]
time_pts = 1:size(met_dynamics, 3);

f = figure(Name='met dynamics');
set(gcf,'Color','white');

% construction
t = tiledlayout(size(met_dynamics, 1), 1);
tissue_names = ["Vasculature", "Gray Matter", "White Matter"];
for i_tissue = 1:size(met_dynamics, 1)
    nexttile;
    hold on;
    for i_met = 1:size(met_dynamics, 2)
        plot(time_pts, squeeze(met_dynamics(i_tissue, i_met, :)), LineWidth=1);
    end
    ax = gca;
    ax.Color = 'white';
    ax.XColor = 'black'; ax.YColor = 'black';
    xlabel('time point');
    ylabel('signal')
    leg = legend(["Pyruvate", "Lactate", "Bicarbonate"], TextColor='black', Color='white', EdgeColor='black');
    title(tissue_names(i_tissue), Color='black');
    hold off;
end
f.Position = [0, 0, 1000, 600];
saveas(f, 'figures/met_dynamics.png');


% met images (no ktrans)
slice = round(size(met_images_no_ktrans, 3) / 2);
time_pts = 1:3:n_t;
f = display_tiled_images(met_images_no_ktrans(:,:,slice,:,time_pts), true, 'met images no ktrans', ["Pyruvate", "Lactate", "Bicarbonate"]);
f.Position = [0, 0, 1440, 475];
saveas(gcf, 'figures/1-met_img_no_ktrans.png');

% met images (ktrans)
slice = round(size(met_images, 3) / 2);
time_pts = 1:3:n_t;
f = display_tiled_images(met_images(:,:,slice,:,time_pts), true, 'met images', ["Pyruvate", "Lactate", "Bicarbonate"]);
f.Position = [0, 0, 1440, 475];
saveas(gcf, 'figures/2-met_img_w_ktrans.png');


% met images of various mri steps --

% augmentations
slice = round(size(met_images_aug, 3) / 2);
time_pts = 1:3:n_t;
f = display_tiled_images(met_images_aug(:,:,slice,:,time_pts), true, 'met images (augmentations)', ["Pyruvate", "Lactate", "Bicarbonate"]);
f.Position = [0, 0, 1440, 475];
saveas(gcf, 'figures/3-met_img_w_augs.png');

% coil limits
slice = round(size(met_images_coil_lim, 3) / 2);
time_pts = 1:3:n_t;
f = display_tiled_images(met_images_coil_lim(:,:,slice,:,time_pts), true, 'met images (coil lim)', ["Pyruvate", "Lactate", "Bicarbonate"]);
f.Position = [0, 0, 1440, 475];
saveas(gcf, 'figures/4-met_img_w_coil_lims.png');

% coil limit maps
slices = 1:10:100;
f = display_tiled_images(coil_sens_weights(:,:,slices), true, 'Coil Sensitivity Map', []);
f.Position = [0, 0, 1440, 185];
saveas(gcf, 'figures/coil_sensitivity.png');

%% multiresolution
slice = round(size(met_images_mres{1}, 3) / 2);
time_pts = 1:3:n_t;
f = display_tiled_images(met_images_mres, true, 'met images multires', ["Pyruvate", "Lactate", "Bicarbonate"], slice, time_pts);
f.Position = [0, 0, 1440, 475];
saveas(gcf, 'figures/5-met_images_mres.png');

% noise
slice = round(size(met_images_mres_noise{1}, 3) / 2);
time_pts = 1:3:n_t;
f = display_tiled_images(met_images_mres_noise, true, 'met images noisy', ["Pyruvate", "Lactate", "Bicarbonate"], slice, time_pts);
f.Position = [0, 0, 1440, 475];
saveas(gcf, 'figures/6-met_images_noise.png');


%% final met images (after mri)

slice = round(size(met_images_upsampled{1}, 3) / 2);
time_pts = 1:3:n_t;
f = display_tiled_images(met_images_upsampled, true, 'final met images!', ["Pyruvate", "Lactate", "Bicarbonate"], slice, time_pts);
f.Position = [0, 0, 1440, 475];
saveas(gcf, 'figures/7-met_images_upsampled.png');


%% AUCs
aucs = cell(1,3);

for i = 1:3
    aucs{i} = sum(met_images_upsampled{i}, 4);
    % turn it into (row, col, 1, slice)
    new_size = cat(2, size(aucs{i}, 1,2), 1, size(aucs{i}, 3));
    aucs{i} = reshape(aucs{i}, new_size);
end

f = display_tiled_images(aucs, true, 'aucs', ["Pyruvate", "Lactate", "Bicarbonate"]);
f.Position = [0, 0, 1440, 565];
saveas(gcf, 'figures/aucs.png');

disp(['took ', num2str(toc), 's', newline]);
disp('done! (plots might take a while to load)');
disp(['took a total of ', num2str(toc(t_all)), 's']);


function [fig, I] = display_tiled_images(I, has_colorbar, figurename, labels, slice, time_pts)
    % I = (row, col, slice). 1 row, slice columns
    % I = (row, col, met, time). met rows, time columns
    % I = {met} --> (row, col, slice, time). met rows, time columns. assumes all have the same number of timesteps

    % stored as:
    % I = {met} --> (row, col, time/slice). met rows, time/slice columns. assumes all have the same number of timesteps

    % slice and time_pts is only for the 3rd option because there's not a great way to just remove that dimension
    
    arguments
        I
        has_colorbar
        figurename
        labels
        slice (1,1) = 1 % these two are optional, won't be used unless I is a cell
        time_pts = NaN
    end

    % keep the dimensions consistent
    % in the form {met} --> (row, col, time)
    if iscell(I)
        if any(isnan(time_pts))
            time_pts = 1:size(I{1}, 4);
        end
        for row = 1:numel(I)
            I{row} = squeeze(I{row}(:, :, slice, time_pts));
        end
    else
        % convert (x,y,z) --> (x,y,1,z)
        I = squeeze(I);
        if numel(size(I)) == 3
            new_shape = cat(2, size(I, 1:2), 1, size(I, 3));
            I = reshape(I, new_shape);
        end
        % turn it into cell
        cell_I = cell([1, size(I, 3)]);
        for row = 1:size(I, 3)
            cell_I{row} = squeeze(I(:,:,row,:));
        end

        I = cell_I;
    end


    % setup
    fig = figure(Name=figurename);
    set(gcf,'Color','white');

    if has_colorbar
        n_cols = size(I{1}, 3) + 1;
    else
        n_cols = size(I{1}, 3);
    end


    % construction
    t = tiledlayout(numel(I), n_cols);
    t.Padding = 'none';
    t.TileSpacing = 'none';
    for row = 1:numel(I)
        scale = [0, max(I{row}(:,:,:), [], 'all')];
        for col = 1:size(I{row}, 3)
            nexttile;
            imshow(I{row}(:,:,col), scale);
            if col == 1 && ~isempty(labels)
                ylabel(labels(row), Color='black');
            end
        end

        if has_colorbar
            cb_ax = nexttile();
            axis(cb_ax, 'off');
            cb = colorbar('location','west');
            clim(scale);
            cb.Color = 'black';
        end
    end

    colormap fire;
end
