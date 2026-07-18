classdef tissue_structure
    properties
        Name string
        Mask (:,:,:,:) {mustBeNumeric} % size = (row, col, slice, tissue)
        Tissues (1,:) string % size = (1, tissue)
        K_trans_map (:,:,:) {mustBeNumeric} % size = (row, col, slice)
    end
    methods
        function obj = tissue_structure(name, mask, tissues, dims)
            % Parameters:
            %   name    = name of tissue structure. size = (1,1), type = string
            %   mask    = tissue masks. size = (row, col, slice, tissue), type = numeric or logical
            %   tissues = names of tissues. size = (1, tissue), type = string. purely semantic right now
            %   dims    = (optional) desired size of mask. [rows, cols, slices]
            arguments
                name (1,1) string
                mask (:,:,:,:) {mustBeNumericOrLogical}
                tissues (1,:) string
                dims (1,3) {mustBeInteger} = [0,0,0]
            end

            if numel(tissues) ~= size(mask, 4)
                error("Mismatched number of tissues in `mask` and `tissues`")
            end

            obj.Name = name;
            obj.Tissues = tissues;

            if islogical(mask)
                % cast to numeric if mask is logical
                obj.Mask = double(mask);
            else
                obj.Mask = mask;
            end

            if dims ~= [0,0,0]
                obj = obj.downscale_mask(dims);
            end

            obj = obj.normalize_mask();
        end


        function obj = create_k_trans_map(obj, k_trans)
            % Parameters:
            %   k_trans         = k_trans values. size = (1, tissue); OR
            %                     (2, tissue) where row 1 = min k_trans, row 2 = max k_trans. if this is the case, k_trans will follow a linear gradient

            % argument validation
            arguments
                obj
                k_trans (:,:) {mustBeNumeric}
            end

            if size(k_trans, 1) > 2
                error('length of 1st dimension of k_trans must be either 1 or 2');
            end

            if size(k_trans, 2) ~= size(obj.Mask, 4)
                error('mismatched array sizes: 2nd dimension of `k_trans` and 4th dimension of `tissue_structure.Mask` must both equal number of tissues')
            end

            % case where k_trans is constant (not gradient)
            % mask = (row, col, slice, tissue)
            if size(k_trans, 1) == 1
                mask = permute(obj.Mask, [4,1,2,3]);
                obj.K_trans_map = squeeze(pagemtimes(k_trans, mask));
                return
            end

            % case where k_trans is a gradient
            n_tissues = size(k_trans, 2);
            mask_size = size(obj.Mask, 1:3);

            k_trans_map = zeros(size(obj.Mask));
            for i_tissue = 1:n_tissues
                gradient = generate_linear_gradient(mask_size, k_trans(1, i_tissue), k_trans(2, i_tissue));
                k_trans_map(:,:,:,i_tissue) = squeeze(obj.Mask(:,:,:,i_tissue)) .* gradient;
            end

            obj.K_trans_map = sum(k_trans_map, 4);

            % helper function from `brainweb_metabolic_phantom`
            function grad = generate_linear_gradient(maskSize, kTRANS_low, kTRANS_high)
                x = linspace(-1, 1, maskSize(1));
                y = linspace(-1, 1, maskSize(2));
                z = linspace(-1, 1, maskSize(3));
                [~, Y, ~] = meshgrid(x, y, z);
                grad = 0.5*(kTRANS_high - kTRANS_low)*Y + 0.5*(kTRANS_low + kTRANS_high);
            end
        end


        function obj = downscale_mask(obj, dims)
            % downscales mask (because large masks can cause issues with memory requests)
            % Parameters:
            %   dims    = desired dimensions. [rows, cols, slices]
            arguments
                obj
                dims (1,3) {mustBeNumeric} = [100,100,100]
            end

            n_tissues = size(obj.Mask, 4);
            new_mask = zeros(cat(2, dims, n_tissues));

            for i_tissue = 1:n_tissues
                new_mask(:,:,:,i_tissue) = imresize3(obj.Mask(:,:,:,i_tissue), dims, "cubic");
            end

            obj.Mask = new_mask;
        end


        function obj = normalize_mask(obj)
            % normalizes voxels so that the sum of the masks of each tissue in any given voxel is no more than 1
            norm_weights = sum(obj.Mask, 4);
            norm_weights(norm_weights < 1) = 1; % only normalize voxels with a sum > 1
            obj.Mask = obj.Mask ./ repmat(norm_weights, [1,1,1, size(obj.Mask,4)]);
        end


        function obj = apply_transforms(obj, tform2d, z_translation)
            % applies transformations to obj.Mask
            % Parameters:
            %   tform2d         = transform to apply to each slice. type = affinetform2d 
            %   z_translation   = z translation
            arguments
                obj
                tform2d (1,1) affinetform2d
                z_translation (1,1) {mustBeNumeric}
            end
            
            output_view = affineOutputView(size(obj.Mask, [1,2]), tform2d, BoundsStyle="CenterOutput");
            for i_tissue = 1:size(obj.Mask, 4)
                obj.Mask(:,:,:,i_tissue) = imwarp(obj.Mask(:,:,:,i_tissue), tform2d, OutputView=output_view);
            end

            % z translation
            mask_ztrans = zeros(size(obj.Mask));
            n_slices = size(obj.Mask,3);
            cutoff = n_slices - abs(z_translation);
            if z_translation > 0
                mask_ztrans(:,:, 1:cutoff, :,:) = obj.Mask(:,:, (z_translation + 1):n_slices ,:,:);
                obj.Mask = mask_ztrans;
            elseif z_translation < 0
                mask_ztrans(:,:, (abs(z_translation) + 1):n_slices, :,:) = obj.Mask(:,:, 1:cutoff, :,:);
                obj.Mask = mask_ztrans;
            end

        end


        function [rgb, colors, fig] = plot_alpha_composite_image(obj, opts)
            % plots color visualization of mask
            % Optional Parameters:
            %   slice   = slice to plot. size = (1,1)
            %   order   = layer order of each tissue. size = (1, tissue)
            % Outputs:
            %   rgb     = rgb alpha-composited volume. size = (row, col, slice, 3)
            %   colors  = rgb triplets of each tissue. size = (tissue, 3)
            %   fig     = the figure that this creates

            % https://en.wikipedia.org/wiki/Alpha_compositing
            % performs the 'over' operation (basically, laying transparent layers on top of each other): 
            %   a_0 = a_a + a_b * (1 - a_a)
            %   C_0 = (C_a * a_a + C_b * a_b * (1 - a_a)) / a_0
            % where a_0, a_a, and a_b are the alpha values of the pixels
            % and C_0, C_a, and C_b are the color components of the pixels

            arguments
                obj
                opts.slice (1,1) {mustBeNumeric} = NaN
                opts.order (1,:) {mustBeNumeric} = NaN
            end

            n_tissues = size(obj.Mask, 4);

            if any(isnan(opts.order))
                opts.order = 1:n_tissues;
            elseif numel(opts.order) ~= n_tissues
                error('mismatched number of tissues in `order` and `tissue_structure.Mask`');
            end

            % find rgb
            norm_mask = obj.Mask ./ max(sum(obj.Mask, 4), [], 'all');

            colors = lines(n_tissues);
            rgb = zeros(cat(2, size(norm_mask, 1:3), 3)); % rgb version of this
            alpha = ones(size(norm_mask, 1:3)); % alpha values (initialized as 1)

            for i_tissue = opts.order
                mask_alpha = norm_mask(:,:,:,i_tissue);
                alpha_0 = mask_alpha + alpha .* (1 - mask_alpha);
                for color_channel = 1:3
                    C_a = colors(i_tissue, color_channel);
                    C_b = rgb(:,:,:,color_channel);
                    alpha_a = mask_alpha;
                    alpha_b = alpha;
                    rgb(:,:,:,color_channel) = (C_a .* alpha_a + C_b .* alpha_b .* (1 - alpha_a)) ./ alpha_0;
                end
                alpha = alpha_0;
            end

            if isnan(opts.slice)
                return
            end

            % plot
            fig = figure;
            imshow(squeeze(rgb(:,:,opts.slice,:)));

            % legend for this
            fake_legend_lines = zeros(1, n_tissues);
            for i_tissue = 1:n_tissues
                fake_legend_lines(i_tissue) = line(NaN, NaN, 'color', colors(i_tissue, :));
            end
            legend(fake_legend_lines, num2cell(obj.Tissues));
        end
    end
end
