classdef metabolite
    properties
        Mz0 (1,1) {mustBeNumeric}
        R1 (1,1) {mustBeNumeric}
        K (1,2) {mustBeNumeric}
    end

    methods
        function obj = metabolite(args) % name=value arguments
            % Parameters:
            %   Mz0   = initial magnetization. size = (1,1)
            %   R1    = relaxation rate. size = (1,1)
            %   k     = kinetic rate(s). size = (1,1) for only forward rate, or (1,2) for forward and reverse rate
            arguments
                args.Mz0 (1,1) {mustBeNumeric} = 0
                args.R1 (1,1) {mustBeNumeric} = 0
                args.k (1,:) {mustBeNumeric} = 0
            end

            obj.Mz0 = args.Mz0;
            obj.R1 = args.R1;
            
            % parse k
            if size(args.k, 2) == 1
                obj.K = [args.k, 0];
            elseif size(args.k, 2) == 2
                obj.K = args.k;
            else
                error("size of `k` must be either (1,1) (forward rate only) or (1,2) (both forward and reverse rates)");
            end

        end
    end
end
