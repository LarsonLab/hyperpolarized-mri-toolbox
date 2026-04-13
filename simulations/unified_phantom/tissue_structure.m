classdef tissue_structure
    properties
        Name string
        Mask {mustBeNumericOrLogical}
        Tissues
    end
    methods
        function obj = tissue_structure(name, mask, tissues)
            % parameters:
            %   name = string
            %   mask = [row, col, slice, tissue]
            obj.Name = name;
            obj.Mask = mask;
            obj.Tissues = tissues;
        end
    end
end
