classdef pk_params
    properties
        Substrate (1,1) metabolite
        Products (1,2) metabolite 
        TR (1,1) {mustBeNumeric}
        InputFunction (1,:) {mustBeNumeric}
        Flips (:,:) {mustBeNumeric}
    end
    methods
        function obj = pk_params(args)
            % Parameters:
            %   substrate       = substrate metabolite. size = (1,1), type = metabolite
            %   products        = product metabolites. size = (1, n_mets - 1), type = metabolite
            %   TR              = repetition time. size = (1,1)
            %   input_function  = additional input of substrate. size = (1, time_pt)
            %   flips           = flip angles (rad). size = (met, time_pt)
            
            % argument validation
            arguments
                args.substrate (1,1) metabolite
                args.products (1,2) metabolite
                args.TR (1,1) {mustBeNumeric}
                args.input_function (1,:) {mustBeNumeric}
                args.flips (:,:) {mustBeNumeric}
            end

            n_mets = numel(args.substrate) + numel(args.products);
            if size(args.flips, 1) ~= n_mets
                error("Mismatched number of metabolites (substrate + number of products) and rows in `flips`");
            end

            nt = numel(args.input_function);
            if size(args.flips, 2) ~= nt
                error("mismatched number of time points in `input_function` and `flips`");
            end

            obj.Substrate = args.substrate;
            obj.Products = args.products;
            obj.TR = args.TR;
            obj.InputFunction = args.input_function;
            obj.Flips = args.flips;
        end

        % various getters
        function Mz0 = get_Mz0(pk_params)
            % Outputs: Mz0 = (1, met)
            Mz0 = pk_params.Substrate.Mz0;
            for met = 1:numel(pk_params.Products)
                Mz0 = cat(2, Mz0, pk_params.Products(met).Mz0);
            end
        end
    
        function R1 = get_R1(pk_params)
            % Outputs: R1 = (1, met)
            R1 = pk_params.Substrate.R1;
            for met = 1:numel(pk_params.Products)
                R1 = cat(2, R1, pk_params.Products(met).R1);
            end
        end

        function flips = get_flips(pk_params)
            % Outputs: flips = (met, time_pt)
            flips = pk_params.Flips;
        end

        function k = get_kinetic_rates(pk_params)
            % Outputs: kinetic rates = (met, fw/rv)
            k = zeros(0,2);
            for met = 1:numel(pk_params.Products)
                k = cat(1, k, pk_params.Products(met).K);
            end
        end
    end
end
