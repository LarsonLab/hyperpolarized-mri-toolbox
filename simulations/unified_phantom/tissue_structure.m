classdef tissue_structure
    properties
        Name string
        Mask {mustBeNumeric}
        Tissues
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
    end
end
