classdef metabolite
    properties
        Mz0 (1,1) {mustBeNumeric}
        R1 (1,1) {mustBeNumeric}
        Flips (1,:) {mustBeNumeric}
        K (1,:) {mustBeNumeric}
    end
    methods
        function obj = metabolite(args)
            arguments
                args.Mz0 (1,1) {mustBeNumeric} = 0
                args.R1 (1,1) {mustBeNumeric} = 0
                args.flips (1,:) {mustBeNumeric} = 0
                args.k (1,:) {mustBeNumeric} = 0
            end
            obj.Mz0 = args.Mz0;
            obj.R1 = args.R1;
            obj.Flips = args.flips;
            obj.K = args.k;
        end
    end
end
