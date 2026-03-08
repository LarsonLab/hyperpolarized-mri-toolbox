classdef mri_system
    methods (Static)         
        function met_images_multires = make_met_images_multires(met_images, sample_size)
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
        

        function met_images_w_noise = add_rician_noise(met_images_mres, SNR)
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
                met_images_mres cell
                SNR (1,:) {mustBeNumeric}
            end
        
            if numel(met_images_mres) ~= size(SNR, 2)
                error("Mismatched array sizes. met_images_multires and SNR must be equal in length. Each element represents a metabolite");
            end
        
            met_images_w_noise = met_images_mres;
            n_mets = size(met_images_mres, 1);
            for Imet = 1:n_mets
                Nt = size(met_images_mres{Imet}, 4);
                sample_size = size(met_images_mres{Imet}, 1:3);
        
                std_noise = max(sum(met_images_mres{Imet}, 4), [], 'all') ./ (SNR(Imet) * sqrt(Nt));
                noise_R = randn(cat(2, sample_size, Nt)) * std_noise; 
                noise_I = randn(cat(2, sample_size, Nt)) * std_noise;
        
                met_images_w_noise{Imet} = sqrt((met_images_mres{Imet} + noise_R).^2 + noise_I.^2);
            end
        end
    end
end
