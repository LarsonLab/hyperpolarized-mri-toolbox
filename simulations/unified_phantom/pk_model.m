classdef pk_model
    methods (Static)
        function met_dynamics = run_pk_model(mz0, r1, k, t_arrival, t_bolus, input_function, flips, tr)
            % arguments:
            %   mz0 = [met, tissue]
            %   r1 = [1, met]
            %   k = [met-1, tissue]
            %   t_arrival = [1, tissue]
            %   t_bolus 
            %   input_function = [tissue, tpt]
            %   flips = [met, tpt]
            %   tr

            % argument validation
            arguments
                mz0 (:,:) {mustBeNumeric}
                r1 (1,:) {mustBeNumeric}
                k (:,:) {mustBeNumeric}
                t_arrival (1,:) {mustBeNumeric}
                t_bolus (1,1) {mustBeNumeric}
                input_function (:,:) {mustBeNumeric}
                flips (:,:) {mustBeNumeric}
                tr (1,1)
            end

            [n_mets, n_tissues, n_tpts] = pk_model.validate_pk_args(mz0, r1, k, t_arrival, input_function, flips);

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
            % arguments:
            %   params = pk_params
            % outputs
            %   met_dynamics = [met, tpt]
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
            % parameters:
            %   met_dynamics = numeric (tissue, met, time_pt)
            %   tissue_struct = tissue_structure
            % outputs:
            %   met_images = (row, col, slice, tissue, met, time_pt)
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
        function [n_mets, n_tissues, n_tpts] = validate_pk_args(mz0, r1, k, t_arrival, input_function, flips)
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
            if size(t_arrival, 2) ~= n_tissues
                error("mismatched number of tissues in `Mz0` and `t_arrival`");
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
