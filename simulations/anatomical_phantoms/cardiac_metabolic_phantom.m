function [kinetic_maps, kTRANS_map, Mz0_maps, input_function_map, met_images_multires] = cardiac_metabolic_phantom(kinetic_rates, kTRANS_scales, Mz0, input_functions, sample_size, output_size, sim_params, heart_idx)
% CARDIAC_METABOLIC_PHANTOM generates standardized 3-dimensional perfusion
%   and metabolism maps for simulated experiments. Supports 3 chemical pool
%   kinetic rate mapping.
%
%   Toolboxes required: Image Processing
%
%   Arguments:
%       kinetic_rates   = kinetic rates to simulate. size = [n_mets - 1, n_tissues]
%       kTRANS_scales   = volume transfer constants for different tissues.
%                         Tissues are, in order, (LV, RV, LV Myocardium, RV
%                         Myocardium). size = [1, n_tissues] or [2, n_tissues]. 
%                         If 2 rows are provided, the generated kTRANS map 
%                         will have a linear gradient.
%       Mz0             = initial magnetization per metabolite per
%                         compartment. size = [n_mets, n_tissues]
%       input_functions = additional input of substrate per tissue per
%                         time point. size = [n_tissues, n_time_pts]
%       sample_size     = desired sampling ("aquisition") matrix size per
%                         metabolite. size = [1,3] or [n_mets,3].
%                         Columns are, in order, (row, col, slice). If only
%                         one row is provided, all metabolites will have
%                         the same sample size.
%       output_size     = desired output matrix size. size = [1,3]
%       sim_params      = parameters used for the kinetic simulations
%                         (struct). Required parameters are:
%                           - TR    = repetition time (s) (scalar).
%                           - R1    = relaxation rates for each metabolite (1/s).
%                                     size = [1, n_mets]
%                           - flips = flip angles per RF pulse per
%                                     metabolite (rad). size = [n_mets, Nt]
%                           - SNR   = signal-to-noise ratio per metabolite.
%                                     size = [1, n_mets]
%
%   Optional Arguments:
%       heart_idx       = heart mask used (1-20). default = 1;
%
%   Outputs:
%       kinetic_maps        = kinetic rate maps of met1->met2 and met1->met3.
%                             size = [row, col, slice, reaction]
%       kTRANS_map          = perfusion map. size = [row, col, slice]
%       Mz0_maps            = initial magnetization maps. size = [row, col, 
%                             slice, metabolite]
%       input_function_map  = map of additional input of substrate over
%                             time. size = [row, col, slice, time_pt]
%       met_images_multires = multiresolution metabolite images (cell array)
%
% Copyright, 2025

    % parse input arguments
    arguments
        kinetic_rates (2,4) {mustBeNumeric}
        kTRANS_scales (:,4) {mustBeNumeric}
        Mz0 (3,4) {mustBeNumeric}
        input_functions (4,:) {mustBeNumeric}
        sample_size (:,3) {mustBeInteger, mustBePositive}
        output_size (1,3) {mustBeInteger, mustBePositive}
        sim_params struct
        heart_idx {mustBeInteger, mustBePositive} = 1
    end

    % load masks
    current_path = pwd;
    mask_path = fullfile(current_path,'util/_src_cardiac', num2str(heart_idx), 'cardiac_masks.mat');   
    load(mask_path,'im_mask'); 

    if size(im_mask, 4) ~= 4
        error('Unexpected number of tissues present in the imported masks');
    end

    % find n_mets and max_sample_size
    n_mets = size(Mz0, 1);
    if size(sample_size, 1) == 1
        max_sample_size = sample_size;
        sample_size = repmat(sample_size, n_mets, 1);
    else
        [~, max_sample_size_row] = max(sample_size(:,1));
        max_sample_size = sample_size(max_sample_size_row, :);
    end
    

    % GENERATE MAPS
    [kinetic_maps, kTRANS_map, Mz0_maps, input_function_map] = generate_maps(im_mask, max_sample_size, kinetic_rates, kTRANS_scales, Mz0, input_functions);
    
    % SIMULATE METABOLITE DYNAMIC IMAGES
    met_images = simulate_metabolite_dynamic_images(kinetic_maps, kTRANS_map, Mz0_maps, input_function_map, sim_params.R1, sim_params.flips, sim_params.TR);

    % MULTIRESOLUTION
    met_images_multires = make_met_images_multires(met_images, sample_size);
    % add rician noise
    met_images_multires = add_rician_noise(met_images_multires, sim_params.SNR);

    % resize to output size
    new_kTRANS_map = zeros(cat(2, output_size, size(kTRANS_map, 4)));
    new_kinetic_maps = zeros(cat(2, output_size, size(kinetic_maps, 4:5)));
    new_Mz0_maps = zeros(cat(2, output_size, size(Mz0_maps, 4:5)));
    for Itissue = 1:size(kTRANS_map, 4)
        new_kTRANS_map(:,:,:,Itissue) = imresize3(kTRANS_map(:,:,:,Itissue), output_size, 'lanczos3');

        for Imet = 1:size(kinetic_maps, 5)
            new_kinetic_maps(:,:,:,Itissue,Imet) = imresize3(kinetic_maps(:,:,:,Itissue,Imet), output_size, 'lanczos3');
        end
    
        for Imet = 1:size(Mz0_maps, 5)
            new_Mz0_maps(:,:,:,Itissue,Imet) = imresize3(Mz0_maps(:,:,:,Itissue,Imet), output_size, 'lanczos3');
        end
    end
    kinetic_maps = new_kinetic_maps;
    Mz0_maps = new_Mz0_maps;
end




function [kinetic_maps, kTRANS_map, Mz0_maps, input_function_map] = generate_maps(im_mask, max_sample_size, kinetic_rates, kTRANS_scales, Mz0, input_functions)
    % generates kTRANS, kinetic, Mz0, and input function maps
    % Arguments:
    %   im_mask         = tissue masks. [row, col, slice, tissue]
    %   max_sample_size = desired maximum matrix size (maybe that's a
    %                     little unclear). [1, dim]
    %   kinetic_rates   = kinetic rates to simulate. (# of chemical pools
    %                     - 1)x(# of tissues)
    %   kTRANS_scales   = volume transfer constants for different tissues.
    %                     2x(n_tissues) where row=1 is min and row=2 is max. 
    %                     If only 1 row is provided, no kTRANS_map will not
    %                     have a gradient
    %   Mz0             = initial magnetization per metabolite per tissue.
    %                     [metabolite, tissue]
    %   input_functions = additional input of substrate per tissue per
    %                     time point. [tissue, time point]
    %
    % Outputs
    %   kinetic_maps    = kinetic rate maps of met1->met2 and met1->met3.
    %                     [row, col, slice, tissue, reaction]
    %   kTRANS_map      = perfusion map. [row, col, slice, tissue]
    %   Mz0_maps        = initial magnetization maps. [row, col, slice,
    %                     tissue, metabolite]
    %   input_function_map = map of additional input of substrate over
    %                        time. [row, col, slice, time_pt, tissue]

    % verify arguments
    % right now, this isn't very flexible, which I do want to change
    arguments
        im_mask (:,:,:,4) {mustBeNumericOrLogical}
        max_sample_size (1,3) {mustBeInteger, mustBePositive}
        kinetic_rates (2,4) {mustBeNumeric}
        kTRANS_scales (:,4) {mustBeNumeric}
        Mz0 (3,4) {mustBeNumeric}
        input_functions (4,:) {mustBeNumeric}
    end

    if size(kTRANS_scales, 1) > 2
        error('kTRANS_scales must have either 1 row (for no linear kTRANS gradient) or 2 rows (for linear kTRANS gradient');
    end

    % constants
    LV = 1;  % left ventricle idx
    RV = 2;  % right ventricle idx
    LMY = 3; % left myocardium idx
    RMY = 4; % right myocardium idx

    % downsample masks to sample_size
    new_im_mask = zeros(cat(2, max_sample_size, size(im_mask, 4)));
    for tissue = 1:size(im_mask, 4)
        new_im_mask(:,:,:,tissue) = imresize3(im_mask(:,:,:,tissue), max_sample_size, "cubic"); 
    end
    new_im_mask(new_im_mask > 1) = 1;
    new_im_mask(new_im_mask < 0) = 0;
    im_mask = new_im_mask;
    clear new_im_mask;

    mask_size = size(im_mask, 1:3);

    % generate kTRANS maps
    if size(kTRANS_scales, 1) == 2
        kTRANS_grad_lv = generate_linear_gradient(mask_size, kTRANS_scales(1,LV), kTRANS_scales(2,LV));
        kTRANS_grad_rv = generate_linear_gradient(mask_size, kTRANS_scales(1,RV), kTRANS_scales(2,RV));
        kTRANS_grad_lmy = generate_linear_gradient(mask_size, kTRANS_scales(1,LMY), kTRANS_scales(2,LMY));
        kTRANS_grad_rmy = generate_linear_gradient(mask_size, kTRANS_scales(1,RMY), kTRANS_scales(2,RMY));
        
        kTRANS_lv = squeeze(im_mask(:,:,:,LV)) .* kTRANS_grad_lv;
        kTRANS_rv = squeeze(im_mask(:,:,:,RV)) .* kTRANS_grad_rv;
        kTRANS_lmy = squeeze(im_mask(:,:,:,LMY)) .* kTRANS_grad_lmy;
        kTRANS_rmy = squeeze(im_mask(:,:,:,RMY)) .* kTRANS_grad_rmy;
    
        %kTRANS_map = (kTRANS_lv + kTRANS_rv + kTRANS_lmy + kTRANS_rmy) ./ sum_weights;
        kTRANS_map = cat(4, kTRANS_lv, kTRANS_rv, kTRANS_lmy, kTRANS_rmy);
    else
        kTRANS_map = create_map(im_mask, kTRANS_scales);
    end

    % generate kinetic maps
    kinetic_1_2_map = create_map(im_mask, kinetic_rates(1,:));
    kinetic_1_3_map = create_map(im_mask, kinetic_rates(2,:));

    % generate Mz0 maps
    Mz0_1_map = create_map(im_mask, Mz0(1,:));
    Mz0_2_map = create_map(im_mask, Mz0(2,:));
    Mz0_3_map = create_map(im_mask, Mz0(3,:));

    % generate input function maps
    input_function_LV_map = squeeze(im_mask(:,:,:,LV)) ...
        .* reshape(input_functions(LV,:), 1, 1, 1, []);
    input_function_RV_map = squeeze(im_mask(:,:,:,RV)) ...
        .* reshape(input_functions(RV,:), 1, 1, 1, []);
    input_function_LMY_map = squeeze(im_mask(:,:,:,LMY)) ...
        .* reshape(input_functions(LMY,:), 1, 1, 1, []);
    input_function_RMY_map = squeeze(im_mask(:,:,:,RMY)) ...
        .* reshape(input_functions(RMY,:), 1, 1, 1, []);

    % consolidate maps
    kinetic_maps = cat(5, kinetic_1_2_map, kinetic_1_3_map);
    Mz0_maps = cat(5, Mz0_1_map, Mz0_2_map, Mz0_3_map);
    input_function_map = cat(5, input_function_LV_map, input_function_RV_map, input_function_LMY_map, input_function_RMY_map);
end


function [met_images_sp] = simulate_metabolite_dynamic_images(kinetic_maps, kTRANS_map, Mz0_maps, input_function_map, R1, flips, TR)
    % simulates metabolite dynamic images
    % Arguments:
    %   kinetic_maps          = kinetic rates per voxel per metabolite. [row, col,
    %                           slice, tissue, reaction] where reaction = n_mets - 1
    %   kTRANS_map            = volume transfer rate per voxel. [row, col, slice, tissue]
    %   Mz0_maps              = Mz0 values per voxel per metabolite. [row, col, slice,
    %                           tissue, metabolite]
    %   input_function_map    = additional input of substrate per voxel per
    %                           timepoint. [row, col, slice, time_pt, tissue]
    %   R1                    = relaxation times per metabolite. [1, metabolite]
    %   flips                 = flip angles per metabolite per RF pulse. [metabolite,
    %                           time_pt]
    %   TR                    = repetition time (s)
    %
    % Outputs:
    %   met_images_sp         = simulated metabolite dynamic images.
    %                           [row, col, slice, metabolite, time_pt]


    % validate arguments
    arguments
        kinetic_maps (:,:,:,4,2) {mustBeNumeric}
        kTRANS_map (:,:,:,4) {mustBeNumeric}
        Mz0_maps (:,:,:,4,3) {mustBeNumeric}
        input_function_map (:,:,:,:,4) {mustBeNumeric}
        R1 (1,3) {mustBeNumeric}
        flips (3,:) {mustBeNumeric}
        TR (1,1) {mustBeNumeric}
    end

    if ~isequal(size(Mz0_maps, 1:3), size(kinetic_maps, 1:3)) | ...
            ~isequal(size(Mz0_maps, 1:3), size(kTRANS_map, 1:3)) | ...
            ~isequal(size(Mz0_maps, 1:3), size(input_function_map, 1:3))
        error("Mismatched array sizes. First 3 dimensions of Mz0_maps, kinetic_maps, kTRANS_map, and input_function_map must be consistent, representing (row, col, slice).");
    end

    if size(input_function_map, 4) ~= size(flips, 2)
        error("Mismatched array sizes. 2nd dimension of flips and 4th dimension of input_function_map must both have length=Nt.")
    end


    % simulate metabolite images
    Nt = size(flips, 2);
    n_mets = size(Mz0_maps, 5);
    sample_size = size(kTRANS_map, 1:3);
    n_tissues = size(kTRANS_map, 4);
    met_images_sp = zeros(cat(2, sample_size, [n_mets, Nt, n_tissues]));

    for Ix = 1:size(met_images_sp, 1)
        for Iy = 1:size(met_images_sp, 2)
            for Iz = 1:size(met_images_sp, 3)
                for Itissue = 1:size(met_images_sp, 6)
                    Mz0_voxel = squeeze(Mz0_maps(Ix, Iy, Iz, Itissue, :)) .';
                    kinetic_rates_voxel = [kinetic_maps(Ix, Iy, Iz, Itissue, 1) 0;
                                           kinetic_maps(Ix, Iy, Iz, Itissue, 2) 0];
                    kTRANS_voxel = kTRANS_map(Ix, Iy, Iz, Itissue);
                    input_function_voxel = (squeeze(input_function_map(Ix, Iy, Iz, :, Itissue)) .') ...
                        .* kTRANS_voxel;
    
                    [met_images_sp(Ix, Iy, Iz, :, :, Itissue), ~] = simulate_Nsite_model(Mz0_voxel, R1, kinetic_rates_voxel, flips, TR, input_function_voxel);
                end
            end
        end
    end

    % combine tissues
    met_images_sp = sum(met_images_sp, 6);
end

function [met_images_multires] = make_met_images_multires(met_images, sample_size)
    % converts single-resolution dynamic metabolite images to multiresolution
    % Arguments:
    %   met_images      = single-resolution dynamic metabolite images. 
    %                     [row, col, slice, metabolite, time_pt]
    %   sample_size     = desired matrix sizes for each metabolite.
    %                     [metabolite, dim]
    %
    % Outputs:
    %   met_images_multires = multiresolution dynamic metabolite images.
    %                         (n_mets)x1 cell array, where size of each 
    %                         cell = [row, col, slice, time_pt]
    arguments
        met_images (:,:,:,:,:) {mustBeNumeric}
        sample_size (:,3) {mustBeInteger, mustBePositive}
    end

    if size(met_images, 4) ~= size(sample_size, 1)
        error("Mismatched array sizes. 4th dimension of met_images and 1st dimension of sample_size must be equal and have length=n_mets")
    end


    n_mets = size(met_images, 4);
    n_time_pts = size(met_images, 5);
    met_images_multires = cell(n_mets, 1);

    for Imet = 1:n_mets
        met_images_multires{Imet} = zeros( ...
            cat(2, sample_size(Imet, :), n_time_pts) ...
        );
        for time_pt = 1:n_time_pts
            met_images_multires{Imet}(:, :, :, time_pt) = ...
                imresize3(squeeze(met_images(:, :, :, Imet, time_pt)), sample_size(Imet, :), 'box'); % box averages surrounding voxels
        end
    end
end

function [met_images_w_noise] = add_rician_noise(met_images_multires, SNR)
    % adds rician noise to multiresolution dynamic metabolite images
    % Parameters:
    %   met_images_multires = multiresolution dynamic metabolite images. 
    %                         (n_mets)x1 cell array, where size of each 
    %                         cell = [row, col, slice, time_pt]
    %   SNR                 = signal-to-noise ratio per metabolite. [1,
    %                         metabolite]
    %
    % Outputs:
    %   met_images_w_noise  = metabolite images with rician noise. 
    %                         (n_mets)x1 cell array, where size of each 
    %                         cell = [row, col, slice, time_pt]

    % argument validation
    arguments
        met_images_multires cell
        SNR (1,:) {mustBeNumeric}
    end

    if numel(met_images_multires) ~= size(SNR, 2)
        error("Mismatched array sizes. met_images_multires and SNR must be equal in length. Each element represents a metabolite");
    end


    met_images_w_noise = met_images_multires;
    n_mets = size(met_images_multires, 1);
    for Imet = 1:n_mets
        Nt = size(met_images_multires{Imet}, 4);
        sample_size = size(met_images_multires{Imet}, 1:3);

        std_noise = max(sum(met_images_multires{Imet}, 4), [], 'all') ./ (SNR(Imet) * sqrt(Nt));
        noise_R = randn(cat(2, sample_size, Nt)) * std_noise; 
        noise_I = randn(cat(2, sample_size, Nt)) * std_noise;

        met_images_w_noise{Imet} = sqrt((met_images_multires{Imet} + noise_R).^2 + noise_I.^2);
    end
end

function map = create_map(mask, rates)
    % Parameters:
    %   mask = [row, col, slice, tissue]
    %   rates = [1 tissue]
    map = zeros(size(mask));
    for tissue = 1:numel(rates)
        map(:, :, :, tissue) = mask(:, :, :, tissue) .* rates(tissue);
    end
end


function [grad] = generate_linear_gradient(mask_size, kTRANS_low, kTRANS_high)
    x = linspace(-1, 1, mask_size(1));
    y = linspace(-1, 1, mask_size(2));
    z = linspace(-1, 1, mask_size(3));
    [~, Y, ~] = meshgrid(x, y, z);
    grad = 0.5*(kTRANS_high - kTRANS_low)*Y + 0.5*(kTRANS_low + kTRANS_high);
end
