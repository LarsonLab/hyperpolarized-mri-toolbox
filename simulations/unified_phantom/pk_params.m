classdef pk_params
    properties
        Substrate (1,1) metabolite
        Products (1,2) metabolite 
        TR (1,1) {mustBeNumeric}
        InputFunction (1,:) {mustBeNumeric}
    end
    methods
        function obj = pk_params(args)
            % substrate = substrate metabolite
            % products = vector of product metabolites (1, nmets - 1)
            % TR = repetition time
            % input_function = vector of additional input of substrate (1,nt)
            arguments
                args.substrate (1,1) metabolite
                args.products (1,2) metabolite
                args.TR (1,1) {mustBeNumeric}
                args.input_function (1,:) {mustBeNumeric}
            end
            obj.Substrate = args.substrate;
            obj.Products = args.products;
            obj.TR = args.TR;
            obj.InputFunction = args.input_function;
        end

        % various getters
        function Mz0 = get_Mz0(pk_params)
            % returns Mz0 = (1, met)
            Mz0 = pk_params.Substrate.Mz0;
            for met = 1:numel(pk_params.Products)
                Mz0 = cat(2, Mz0, pk_params.Products(met).Mz0);
            end
        end
    
        function R1 = get_R1(pk_params)
            % returns R1 = (1, met)
            R1 = pk_params.Substrate.R1;
            for met = 1:numel(pk_params.Products)
                R1 = cat(2, R1, pk_params.Products(met).R1);
            end
        end

        function flips = get_flips(pk_params)
            % returns flips = (met, time_pt)
            flips = pk_params.Substrate.Flips;
            for met = 1:numel(pk_params.Products)
                flips = cat(1, flips, pk_params.Products(met).Flips);
            end
        end

        function k = get_kinetic_rates(pk_params)
            % returns kinetic rates = (met, fw/rv)
            k = zeros(0,2);
            for met = 1:numel(pk_params.Products)
                k = cat(1, k, pk_params.Products(met).K);
            end
        end
    end
end
