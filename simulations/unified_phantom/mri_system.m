classdef mri_system
    methods (Static)
        function met_images_mres = run_mri_system(met_images, sample_size, snr)
            % Wrapper for easy use of mri system model
            % Parameters:
            %   met_images      = single-resolution metabolite images. size = (row, col, slice, met, tpt)
            %   sample_size     = desired size of multires images. size = (met, dim)
            %   snr             = signal-to-noise ratio. size = (1, met)
            % Outputs:
            %   met_images_mres = multiresolution metabolite images. 
            %   size = [1, met] cell array. size of each cell = (row, col, slice, time_pt)
            met_images_mres = mri_system.make_met_images_multires(met_images, sample_size);
            met_images_mres = mri_system.add_rician_noise(met_images_mres, snr);
        end

        function met_images_multires = make_met_images_multires(met_images, sample_size)
            % Converts single-resolution dynamic metabolite images to multiresolution dynamic metabolite images
            % Parameters:
            %   met_images      = single-resolution dynamic metabolite images. 
            %                     size = (row, col, slice, metabolite, time_pt)
            %   sample_size     = desired matrix sizes for each metabolite.
            %                     size = (metabolite, dim)
            %
            % Outputs:
            %   met_images_multires = multiresolution dynamic metabolite images.
            %                         size = (1, met) cell array, where size of each cell = (row, col, slice, time_pt)
            arguments
                met_images (:,:,:,:,:) {mustBeNumeric}
                sample_size (:,3) {mustBeInteger, mustBePositive}
            end
        
            if size(met_images, 4) ~= size(sample_size, 1)
                error("Mismatched array sizes. 4th dimension of met_images and 1st dimension of sample_size must be equal and have length=n_mets")
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
        

        function met_images_w_noise = add_rician_noise(met_images, snr)
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
                error("Mismatched array sizes. met_images_multires and SNR must both have length = n_mets");
            end
        
            met_images_w_noise = met_images;
            n_mets = size(met_images, 1);
            for imet = 1:n_mets
                met_images_w_noise{imet} = mri_system.add_rician_noise_image(met_images{imet}, snr(1));
            end
        end
    end

    methods (Static, Access = private)
        function met_image_w_noise = add_rician_noise_image(met_image, snr)
            % Adds rician noise to single dynamic metabolite image
            % Parameters:
            %   met_image = single met image. size = (row, col, slice, time_pt)
            %   snr = signal-to-noise ratio. size = (1,1)
            % Outputs:
            %   met_image_w_noise = metabolite image with noise. size = (row, col, slice, time_pt)

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
        end
    end
end
