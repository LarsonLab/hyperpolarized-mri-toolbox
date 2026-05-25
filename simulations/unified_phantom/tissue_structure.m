classdef tissue_structure
    properties
        Name string
        Mask (:,:,:,:) {mustBeNumeric}
        Tissues
        K_trans_map (:,:,:,:) {mustBeNumeric}
    end
    methods
        function obj = tissue_structure(name, mask, tissues)
            % Parameters:
            %   name    = name of tissue structure. size = (1,1), type = string
            %   mask    = tissue masks. size = (row, col, slice, tissue), type = numeric or logical
            %   tissues = names of tissues. size = (1, tissue), type = string. purely semantic right now
            arguments
                name (1,1) string
                mask (:,:,:,:) {mustBeNumericOrLogical}
                tissues (1, :) string
            end

            obj.Name = name;

            if islogical(mask)
                % cast to numeric if mask is logical
                obj.Mask = double(mask);
            else
                obj.Mask = mask;
            end

            obj.Tissues = tissues;
        end


        function obj = create_k_trans_map(obj, k_trans)
            % Parameters
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
                [X, Y, Z] = meshgrid(x, y, z);
                grad = 0.5*(kTRANS_high - kTRANS_low)*Y + 0.5*(kTRANS_low + kTRANS_high);
            end
        end
    end
end
