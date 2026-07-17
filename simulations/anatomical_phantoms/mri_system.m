classdef mri_system
    methods (Static)
        function [met_images_mres, coil_sens_weights] = run_mri_system(met_images, sample_size, snr, coil_lim, tissue_mask, output_size, augment_params)
            % Wrapper for easy use of mri system model
            % Parameters:
            %   met_images      = single-resolution metabolite images. size = (row, col, slice, met, tpt)
            %   sample_size     = desired size of multires images. size = (met, dim)
            %   snr             = signal-to-noise ratio. size = (1, met)
            %   coil_lim        = coil limits. [min, max]
            %   tissue_mask     = tissue mask. size = (row, col, slice, tissue)
            %   output_size     = desired output_size. size = (met, dim)
            %   augment_params  = augmentation parameters. struct
            % Outputs:
            %   met_images_mres = multiresolution metabolite images. 
            %                     size = (1, met) cell array. size of each cell = (row, col, slice, time_pt)
            %   coil_sens_weights   = coil sensitivity weights


            % argument validation
            arguments
                met_images (:,:,:,:,:) {mustBeNumeric}
                sample_size (:,3) {mustBeNumeric}
                snr (1,:) {mustBeNumeric}
                coil_lim (1,2) {mustBeNumeric}
                tissue_mask (:,:,:,:) {mustBeNumeric}
                output_size (:,3) {mustBeNumeric} = NaN
                augment_params struct = struct()
            end

            n_mets = size(met_images, 4);
            if size(sample_size, 1) ~= n_mets && size(sample_size, 1) ~= 1
                error('Mismatched array sizes: the 1st dimension of `sample_size` must equal n_mets OR 1');
            end
            if size(snr, 2) ~= n_mets
                error('mismatched array sizes: the 2nd dimension of `snr` must equal n_mets');
            end
            if (size(output_size, 1) ~= n_mets) && (size(output_size, 1) ~= 1) && (all(~isnan(output_size), 'all'))
                error('mismatched array sizes: the 1st dimension of `output_size` must equal n_mets OR 1');
            end

            % run the system
            if ~isempty(fieldnames(augment_params))
                cell_augment_params = namedargs2cell(augment_params); % unpack the augmentation parameters
                met_images = mri_system.augment(met_images, cell_augment_params{:});
            end
            [met_images, coil_sens_weights] = mri_system.apply_coil_lim(met_images, coil_lim, tissue_mask);
            met_images_mres = mri_system.make_met_images_multires(met_images, sample_size);
            met_images_mres = mri_system.add_rician_noise(met_images_mres, snr);

            if all(~isnan(output_size), 'all')
                met_images_mres = mri_system.upsample_to_output_size(met_images_mres, output_size);
            end
        end

        function met_images_multires = make_met_images_multires(met_images, sample_size)
            % Converts single-resolution dynamic metabolite images to multiresolution dynamic metabolite images
            % Parameters:
            %   met_images      = single-resolution dynamic metabolite images. 
            %                     size = (row, col, slice, metabolite, time_pt)
            %   sample_size     = desired matrix sizes for each metabolite.
            %                     size = (metabolite, dim) OR (1, dim)
            %
            % Outputs:
            %   met_images_multires = multiresolution dynamic metabolite images.
            %                         size = (1, met) cell array, where size of each cell = (row, col, slice, time_pt)
            arguments
                met_images (:,:,:,:,:) {mustBeNumeric}
                sample_size (:,3) {mustBeInteger, mustBePositive}
            end

            n_mets = size(met_images, 4);

            if size(sample_size, 1) == 1
                sample_size = repmat(sample_size, [n_mets, 1]);
            elseif size(sample_size, 1) ~= n_mets
                error('Mismatched array sizes: the 1st dimension of `sample_size` must equal n_mets OR 1');
            end
        
            n_mets = size(met_images, 4);
            n_time_pts = size(met_images, 5);
            met_images_multires = cell(1, n_mets);
        
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
        

        function [met_images_w_noise, met_images_w_noise_no_bg] = add_rician_noise(met_images, snr)
            % Adds rician noise to dynamic metabolite images
            % Parameters:
            %   met_images = dynamic metabolite images. either:
            %                - (1, met) cell array, where size of each cell = (row, col, slice, time_pt); OR
            %                - numerical array. size = (met, row, col, slice, time_pt)
            %   snr        = signal-to-noise ratio per metabolite. (1, metabolite)
            %
            % Outputs:
            %   met_images_w_noise  = metabolite images with rician noise. 
            %                         size = (1, met) cell array, where size of each 
            %                         cell = (row, col, slice, time_pt)
            %   met_images_w_noise_no_bg = metabolite images with rician noise, without background noise
            %                         size = (1, met) cell array, where size of each 
            %                         cell = (row, col, slice, time_pt)
        
            % argument validation/input parsing
            arguments
                met_images
                snr (1,:) {mustBeNumeric}
            end

            % convert to cell array if not already a cell array
            if isnumeric(met_images) 
                if numel(size(met_images)) ~= 5
                    error("met_images must either be a numeric array with size = [nmets, row, col, slice, time_pt] OR be a [1, nmets] cell array where each cell = [row, col, slice, time_pt],");
                end

                n_mets = size(met_images, 1);
                cell_met_images = cell(1, n_mets);
                for imet = 1:n_mets
                    single_met_image = met_images(imet,:,:,:,:);
                    single_met_image = reshape(single_met_image, size(single_met_image, 2:numel(size(single_met_image)))); % remove first dimension (without squeeze() lest another dim = 1)

                    cell_met_images{imet} = single_met_image;
                end

                met_images = cell_met_images;
            % error if not either numeric or cell
            elseif ~iscell(met_images)
                error("met_images must either be a numeric array with size = [nmets, row, col, slice, time_pt] OR be a [1, nmets] cell array where each cell = [row, col, slice, time_pt],");
            end

            % error if number of metabolites don't match in met_images and snr
            if numel(met_images) ~= size(snr, 2)
                error("mismatched array sizes. met_images_multires and SNR must both have length = n_mets");
            end
        
            met_images_w_noise = met_images;
            met_images_w_noise_no_bg = met_images;
            n_mets = numel(met_images);
            for imet = 1:n_mets
                [image_w_bg, image_w_out_bg] = mri_system.add_rician_noise_image(met_images{imet}, snr(imet));
                met_images_w_noise{imet} = image_w_bg;
                met_images_w_noise_no_bg{imet} = image_w_out_bg;
            end
        end

        function output_met_images = upsample_to_output_size(met_images_lowres, output_size)
            % Upsamples dynamic metabolite images to output size
            % Parameters:
            %   met_images_lowres   = un-resized metabolite images. 
            %                         size = (1, met) cell array, where each cell = (row, col, slice, time_pt)
            %   output_size         = desired output matrix size. size = (met, dim) OR (1, dim)
            % Outputs:
            %   output_met_images   = upsampled met images
            arguments
                met_images_lowres
                output_size (:,3) {mustBeNumeric}
            end

            n_mets = numel(met_images_lowres);
            if size(output_size, 1) == 1
                output_size = repmat(output_size, [n_mets, 1]);
            elseif size(output_size, 1) ~= n_mets
                error("mismatched array sizes. 1st dimension of `output_size` should both equal n_mets OR 1");
            end

            % upsample
            output_met_images = cell(1, n_mets);
            for i_met = 1:n_mets
                n_tpts = size(met_images_lowres{i_met}, 4);
                met_output_size = output_size(i_met, :);
                output_met_images{i_met} = zeros(cat(2, met_output_size, n_tpts));
                for i_tpt = 1:n_tpts
                    output_met_images{i_met}(:,:,:,i_tpt) = imresize3(met_images_lowres{i_met}(:,:,:,i_tpt), met_output_size);
                end
            end
        end

        function augmented_met_images = augment(met_images, augs)
            % Augments image
            % Positional Parameters:
            %   met_images      = metabolite images. size = (row, col, slice, met, tpt)
            % Name-Value Parameters: see randomAffine2d
            arguments
                met_images (:,:,:,:,:) {mustBeNumeric}
                augs.XTranslation (1,2) {mustBeNumeric} = [0,0]
                augs.YTranslation (1,2) {mustBeNumeric} = [0,0]
                augs.ZTranslation (1,2) {mustBeNumeric} = [0,0]
                augs.Rotation (1,2) {mustBeNumeric} = [0,0]
                augs.Scale (1,2) {mustBeNumeric} = [1,1]
                augs.XReflection (1,1) {mustBeNumericOrLogical} = false
                augs.YReflection (1,1) {mustBeNumericOrLogical} = false
                augs.XShear (1,2) {mustBeNumeric} = [0,0]
                augs.YShear (1,2) {mustBeNumeric} = [0,0]
                augs.Seed (1,1) {mustBeNumeric} = NaN
            end

            % initialize rng
            if ~isnan(augs.Seed)
                rng(augs.Seed);
            else
                rng("shuffle");
            end

            tform = randomAffine2d(...
                XTranslation=augs.XTranslation, ...
                YTranslation=augs.YTranslation, ...
                Rotation=augs.Rotation, ...
                Scale=augs.Scale, ...
                XReflection=augs.XReflection, ...
                YReflection=augs.YReflection, ...
                XShear=augs.XShear, ...
                YShear=augs.YShear ...
            );

            output_view = affineOutputView(size(met_images, [1,2]), tform, BoundsStyle="CenterOutput");
            augmented_met_images = zeros(size(met_images));
            for i_met = 1:size(met_images, 4)
                for i_tpt = 1:size(met_images, 5)
                    augmented_met_images(:,:,:,i_met,i_tpt) = imwarp(met_images(:,:,:,i_met,i_tpt), tform, OutputView=output_view);
                end
            end

            % z translation
            offset = randi(augs.ZTranslation);
            aug_met_img_ztrans = zeros(size(augmented_met_images));
            n_slices = size(met_images,3);
            cutoff = n_slices - abs(offset);
            if offset > 0
                aug_met_img_ztrans(:,:, 1:cutoff, :,:) = augmented_met_images(:,:, (offset + 1):n_slices ,:,:);
                augmented_met_images = aug_met_img_ztrans;
            elseif offset < 0
                aug_met_img_ztrans(:,:, (abs(offset) + 1):n_slices, :,:) = augmented_met_images(:,:, 1:cutoff, :,:);
                augmented_met_images = aug_met_img_ztrans;
            end
        end

        function [met_images_w_coil_lim, coil_sens_weights] = apply_coil_lim(met_images, coil_lim, tissue_mask)
            % Parameters:
            %   met_images  = metabolite images. size = (row, col, slice, met, time_pt)
            %   coil_lim    = coil sensitivity limits. [min, max]
            %   tissue_mask = tissue mask. size = (row, col, slice, tissue)
            % Outputs:
            %   met_images_w_coil_lim   = metabolite images with coil sensitivity. size = (row, col, slice, met, time_pt)
            %   coil_sens_weights       = coil sensitivity weights

            % argument validation
            arguments
                met_images (:,:,:,:,:) {mustBeNumeric}
                coil_lim (1,2) {mustBeNumeric}
                tissue_mask (:,:,:,:) {mustBeNumeric}
            end

            if any(size(met_images, 1:3) ~= size(tissue_mask, 1:3))
                error("mismatched array sizes. `met_images` and `tissue_mask` must have the same number of rows, cols, and slices");
            end

            mask = sum(tissue_mask, 4);
            coil_sens_weights = mri_system.coil_dist_map(mask, coil_lim);

            n_mets = size(met_images, 4);
            n_tpts = size(met_images, 5);
            expanded_coil_sens_weights = repmat(coil_sens_weights, [1, 1, 1, n_mets, n_tpts]);
            
            met_images_w_coil_lim = met_images .* expanded_coil_sens_weights;
        end
    end

    methods (Static, Access = private)
        function [met_image_w_noise, met_image_w_noise_no_bg] = add_rician_noise_image(met_image, snr)
            % Adds rician noise to single dynamic metabolite image
            % Parameters:
            %   met_image = single met image. size = (row, col, slice, time_pt)
            %   snr = signal-to-noise ratio. size = (1,1)
            % Outputs:
            %   met_image_w_noise = metabolite image with noise. size = (row, col, slice, time_pt)
            %   met_image_w_noise_no_bg = metabolite images with noise, without background noise. size = (row, col, slice, time_pt)

            arguments
                met_image (:,:,:,:) {mustBeNumeric}
                snr (1,1) {mustBeNumeric}
            end

            nt = size(met_image, 4);
            sample_size = size(met_image, 1:3);
    
            std_noise = max(sum(met_image, 4), [], 'all') ./ (snr * sqrt(nt));
            noise_R = randn(cat(2, sample_size, nt)) * std_noise; 
            noise_I = randn(cat(2, sample_size, nt)) * std_noise;

            met_image_w_noise = sqrt((met_image + noise_R).^2 + noise_I.^2);

            met_image_w_noise_no_bg = met_image_w_noise;
            met_image_w_noise_no_bg(met_image == 0) = 0;
        end

        % taken from brainweb
        function [weights] = coil_dist_map(mask, lim)
            maskSize = size(mask);
            weights = zeros(maskSize); 

            % create y gradient
            x = linspace(-1, 1, maskSize(1));
            y = linspace(-1, 1, maskSize(2));
            z = linspace(-1, 1, maskSize(3));
            [~, Y, Z] = meshgrid(x, y, z);

            % y gradient
            %lim = [0.6, 1.2];
            grady = 0.5*(lim(2) - lim(1))*Y + 0.5*(lim(1) + lim(2));

            % z gradient
            gradz = (1 - abs(Z).^2) + 0.6;
            gradz = gradz ./ max(gradz, [], 'all');

            grad = grady .* gradz;
            
            for z=1:maskSize(3)
                mask_sl = squeeze(mask(:,:,z));
                
                % get outline/perim of mask
                mask_sl = imfill(bwmorph(bwareaopen(mask_sl,300),"fill"),"holes");
                %figure, imagesc(mask); axis off square;
                bw2 = bwperim(mask_sl);
                %figure, imagesc(bw2)

                %reverse_mask
                mask_rev = imcomplement(mask_sl);
                
                % calculate weights based on distance from mask perim
                w = bwdist(bw2) .^0.5;
                weights(:,:,z) = ((1 - (w ./max(w(:)))) .* grad(:,:,z) .* mask_sl) + mask_rev;
                %figure, imagesc(weights)
            end

            weights(isnan(weights)) = 0;
        end
    end
end
