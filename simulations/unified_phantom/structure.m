classdef structure
    properties
        Name string
        Mask {mustBeNumericOrLogical}
        Tissues
    end
    methods
        function obj = structure(name, mask, tissues)
            % parameters:
            %   name = string
            %   mask = [row, col, slice, tissue]
            obj.Name = name;
            obj.Mask = mask;
            obj.Tissues = tissues;
        end
    end
end
