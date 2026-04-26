classdef pk_model
    methods (Static)
        function met_dynamics = run_pk_model(mz0, r1, k, flips, tr, input_function, t_bolus)
            % Wrapper for easy use of pk model. Accepts either input_function OR t_arrival and t_bolus
            % Parameters:
            %   mz0             = initial magnetization of each metabolite in each tissue. size = (met, tissue)
            %   r1              = relaxation rates. size = (1, met)
            %   k               = forward kinetic rates of product metabolites. size = (met-1, tissue)
            %   flips           = flip angles (rad). size = (met, time_pt)
            %   tr              = temporal resolution. size = (1,1)
            %   input_function  = additional input for substrate per tissue per time point. size = (tissue, tpt)
            %                     if t_bolus is provided, input_function acts as t_arrival.
            %                     t_arrival = arrival time of substrate. size = (1, tissue)
            %   t_bolus         = (optional) time it takes for bolus to enter. size = (1,1)
            % Outputs:
            %   met_dynamics    = metabolite dynamics. size = (tissue, met, time_pt)

            % argument validation
            arguments
                mz0 (:,:) {mustBeNumeric}
                r1 (1,:) {mustBeNumeric}
                k (:,:) {mustBeNumeric}
                flips (:,:) {mustBeNumeric}
                tr (1,1) {mustBeNumeric}
                input_function (:,:) {mustBeNumeric}
                t_bolus (1,1) {mustBeNumeric} = NaN
            end

            % case if t_arrival and t_bolus are provided
            if ~isnan(t_bolus)
                t_arrival = input_function; % just for readability

                % verify t_arrival and t_bolus
                n_tissues = size(mz0, 2);
                if size(t_arrival, 1) ~= 1
                    error("`t_arrival` must have a size = (1, tissue)");
                end
                if size(t_arrival, 2) ~= n_tissues
                    error("mismatched number of tissues in `Mz0` and `t_arrival`");
                end

                % create input_function
                n_tpts = size(flips, 2);
                input_function = zeros(n_tissues, n_tpts);
                for i_tissue = 1:n_tissues
                    input_function(i_tissue, :) = realistic_input_function(n_tpts, tr, t_arrival(i_tissue), t_bolus);
                end
            end
            
            % validate the rest of the arguments
            [n_mets, n_tissues, n_tpts] = pk_model.validate_pk_args(mz0, r1, k, flips, input_function);

            % migrate everything over to compartment-specific pk params
             met_dynamics = zeros(n_tissues, n_mets, n_tpts);
            for i_tissue = 1:n_tissues
                 % get substrate
                 substrate = metabolite( ...
                     Mz0=mz0(1,i_tissue), ...
                     R1=r1(1), ...
                     k=[0,0]);
         
                 % get products
                 products(1, n_mets-1) = metabolite;
                 for i_met = 2:n_mets
                     products(i_met-1) = metabolite( ...
                         Mz0=mz0(i_met, i_tissue), ...
                         R1=r1(i_met), ...
                         k=[k(i_met - 1, i_tissue), 0]);
                 end

                 % put it all into pk_params
                 tissue_pk_params = pk_params(...
                     substrate=substrate, ...
                     products=products, ...
                     TR=tr, ...
                     input_function=input_function(i_tissue, :), ...
                     flips=flips);
         
                 % generate met dynamics
                 met_dynamics(i_tissue,:,:) = pk_model.generate_met_dynamics(tissue_pk_params);
            end
        end


        function met_dynamics = generate_met_dynamics(params)
            % Generates metabolite dynamics from PK parameters
            % Parameters:
            %   params          = pk parameters. size = (1,1), type = pk_params
            % Outputs:
            %   met_dynamics    = metabolite dynamics. size = (met, tpt)
            arguments
                params pk_params
            end
        
            % prep everything for simulate_Nsite_model()
            Mz0 = params.get_Mz0();
            R1 = params.get_R1();
            k = params.get_kinetic_rates();
            flips = params.get_flips();
            TR = params.TR;
            input_function = params.InputFunction;
        
            [met_dynamics, ~] = simulate_Nsite_model(Mz0, R1, k, flips, TR, input_function);
        end


        function met_images = generate_met_images(tissue_struct, met_dynamics)
            % Parameters:
            %   met_dynamics    = metabolite dynamics. size = (tissue, met, time_pt)
            %   tissue_struct   = tissue_structure. size = (1,1), type = tissue_structure
            % Outputs:
            %   met_images      = dynamic metabolite images. 
            %                     size = (row, col, slice, tissue, met, time_pt)
            arguments
                tissue_struct (1,1) tissue_structure
                met_dynamics {mustBeNumeric}
            end
            
            % validate args
            if size(tissue_struct.Mask, 4) ~= size(met_dynamics,1)
                error("mismatched number of tissues in tissue structure and met dynamics");
            end
        
            time_pts = size(met_dynamics, 3);
            n_mets = size(met_dynamics, 2);
            
            % reshape + expand met dynamic
            met_dynamics = reshape(met_dynamics, cat(2,[1,1,1],size(met_dynamics)));
            met_dynamics = repmat(met_dynamics, cat(2, size(tissue_struct.Mask, 1:3), [1,1,1]));
            
            % expand mask
            mask = tissue_struct.Mask; % = (row, col, slice, tissue)
            mask = repmat(mask, cat(2,[1,1,1,1,n_mets,time_pts]));
            
            % multiply
            met_images = met_dynamics .* mask;

            % sum across tissues
            met_images = squeeze(sum(met_images, 4));
        end
    end


    methods (Static, Access = private)
        function [n_mets, n_tissues, n_tpts] = validate_pk_args(mz0, r1, k, flips, input_function, t_arrival)
            % validates arguments for `run_pk_model`
            n_mets = size(mz0, 1);
            if size(r1, 2) ~= n_mets
                error("mismatched number of metabolites in `Mz0` and `R1`");
            end
            if size(flips, 1) ~= n_mets
                error("mismatched number of metabolites in `Mz0` and `flips`");
            end
            if size(k, 1) ~= n_mets - 1
                error("mismatched number of metabolits in `Mz0` and `k`. `k` should have n_mets-1 rows");
            end
            
            n_tissues = size(mz0, 2);
            if size(k, 2) ~= n_tissues
                error("mismatched number of tissues in `Mz0` and `k`");
            end
            if size(input_function, 1) ~= n_tissues
                error("mismatched number of tissues in `Mz0` and `input_function`");
            end
            
            n_tpts = size(flips, 2);
            if size(input_function, 2) ~= n_tpts;
                error("mismatched number of time points in `flips` and `input_function`");
            end
        end
    end
end
