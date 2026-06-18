classdef pk_model
    methods (Static)
        function [met_images, dynamics_low_ktrans, dynamics_high_ktrans, images_low_ktrans, images_high_ktrans] = run_pk_model(mz0, r1, k, flips, tr, tissue_struct, opts)
            % Wrapper for easy use of pk model.
            % Parameters:
            %   mz0             = initial magnetization of each metabolite in each tissue. size = (met, tissue)
            %   r1              = relaxation rates. size = (1, met)
            %   k               = forward kinetic rates of product metabolites. size = (met-1, tissue)
            %   flips           = flip angles (rad). size = (met, time_pt)
            %   tr              = temporal resolution. size = (1,1)
            %   tissue_struct   = tissue_structure. size = (1,1), type = tissue_structure
            % Additional options:
            %   input_function  = additional input for substrate per tissue per time point. size = (tissue, tpt). may not be provided if `t_arrival` and `t_bolus` are provided.
            %   t_arrival       = arrival time of substrate. size = (1, tissue). must be provided with `t_bolus`.
            %   t_bolus         = time it takes for bolus to enter. size = (1,1). must be provided with `t_arrival`
            %   plot            = boolean flag to create plots at each step
            %   met_names       = list of met names. only used for plotting. size = (1,met).
            % Outputs:
            %   met_images      = dynamic metabolite images. size = (row, col, slice, tissue, met, time_pt)
            %   dynamics_low_ktrans     = metabolite dynamices with k_trans = 0. size = (tissue, met, time_pt)
            %   dynamics_high_ktrans    = metabolite dynamices with k_trans = 1. size = (tissue, met, time_pt)
            %   images_low_ktrans       = metabolite images with k_trans = 0. size = (row, col, slice, tissue, met, time_pt)
            %   images_high_ktrans      = metabolite images with k_trans = 1. size = (row, col, slice, tissue, met, time_pt)
            arguments
                mz0 (:,:) {mustBeNumeric}
                r1 (1,:) {mustBeNumeric}
                k (:,:) {mustBeNumeric}
                flips (:,:) {mustBeNumeric}
                tr (1,1) {mustBeNumeric}
                tissue_struct (1,1) tissue_structure
                opts.input_function (:,:) {mustBeNumeric} = NaN
                opts.t_arrival (1,:) {mustBeNumeric} = NaN
                opts.t_bolus (1,1) {mustBeNumeric} = NaN
                opts.plot (1,1) {mustBeNumericOrLogical} = false
                opts.met_names (1,:) string = []
            end
            n_tissues = size(mz0, 2);
            dynamics_low_ktrans = pk_model.generate_all_met_dynamics(mz0, r1, k, zeros(1, n_tissues), flips, tr, ...
                input_function=opts.input_function, t_arrival=opts.t_arrival, t_bolus=opts.t_bolus, plot=opts.plot, tissue_names=tissue_struct.Tissues, met_names=opts.met_names);
            dynamics_high_ktrans = pk_model.generate_all_met_dynamics(mz0, r1, k, ones(1, n_tissues), flips, tr, ...
                input_function=opts.input_function, t_arrival=opts.t_arrival, t_bolus=opts.t_bolus, plot=opts.plot, tissue_names=tissue_struct.Tissues, met_names=opts.met_names);

            images_low_ktrans = pk_model.generate_met_images(tissue_struct, dynamics_low_ktrans);
            images_high_ktrans = pk_model.generate_met_images(tissue_struct, dynamics_high_ktrans);

            met_images = pk_model.apply_k_trans(images_low_ktrans, images_high_ktrans, tissue_struct.K_trans_map);
        end


        function met_dynamics = generate_all_met_dynamics(mz0, r1, k, k_trans, flips, tr, opts)
            % Wrapper for easy use of generate_met_dynamics(). Accepts either input_function OR t_arrival and t_bolus
            % Parameters:
            %   mz0             = initial magnetization of each metabolite in each tissue. size = (met, tissue)
            %   r1              = relaxation rates. size = (1, met)
            %   k               = forward kinetic rates of product metabolites. size = (met-1, tissue)
            %   k_trans         = volumetric transfer constants. size = (1, tissue)
            %   flips           = flip angles (rad). size = (met, time_pt)
            %   tr              = temporal resolution. size = (1,1)
            % Additional options:
            %   input_function  = additional input for substrate per tissue per time point. size = (tissue, tpt). may not be provided if `t_arrival` and `t_bolus` are provided.
            %   t_arrival       = arrival time of substrate. size = (1, tissue). must be provided with `t_bolus`.
            %   t_bolus         = time it takes for bolus to enter. size = (1,1). must be provided with `t_arrival`
            %   plot            = boolean flag to create plots at each step
            %   tissue_names    = list of tissue names. only used for plotting. size = (1,tissue).
            %   met_names       = list of met names. only used for plotting. size = (1,met).
            % Outputs:
            %   met_dynamics    = metabolite dynamics. size = (tissue, met, time_pt)

            % argument validation
            arguments
                mz0 (:,:) {mustBeNumeric}
                r1 (1,:) {mustBeNumeric}
                k (:,:) {mustBeNumeric}
                k_trans (1,:) {mustBeNumeric}
                flips (:,:) {mustBeNumeric}
                tr (1,1) {mustBeNumeric}
                opts.input_function (:,:) {mustBeNumeric} = NaN
                opts.t_arrival (1,:) {mustBeNumeric} = NaN
                opts.t_bolus (1,1) {mustBeNumeric} = NaN
                opts.plot (1,1) {mustBeNumericOrLogical} = false
                opts.tissue_names (1,:) = []
                opts.met_names (1,:) = []
            end

            % input function, t_arrival, t_bolus parsing
            provided_opts = [~any(isnan(opts.input_function), 'all'), ~any(isnan(opts.t_arrival), 'all'), ~isnan(opts.t_bolus)]; % e.g. if only input function is provided, this is [1,0,0]

            if provided_opts == [1,0,0] % case when only input_function is provided
                input_function = opts.input_function;
            elseif provided_opts == [0,1,1] % case when t_arrival and t_bolus are provided, but not input_function
                % verify t_arrival and t_bolus
                n_tissues = size(mz0, 2);
                n_tpts = size(flips, 2);
                if size(opts.t_arrival, 1) ~= 1
                    error("`t_arrival` must have a size = (1, tissue)");
                end
                if size(opts.t_arrival, 2) ~= n_tissues
                    error("mismatched number of tissues in `Mz0` and `t_arrival`");
                end
                % create input_function
                input_function = zeros(n_tissues, n_tpts);
                for i_tissue = 1:n_tissues
                    input_function(i_tissue, :) = realistic_input_function(n_tpts, tr, opts.t_arrival(i_tissue), opts.t_bolus);
                end
            elseif provided_opts == [0,0,0] % case when nothing is provided
                n_tissues = size(mz0, 2);
                n_tpts = size(flips, 2);
                input_function = zeros(n_tissues, n_tpts);
            else
                error("Provide EITHER `input_function` OR both `t_arrival` and `t_bolus`")
            end
            
            % validate the rest of the arguments
            [n_mets, n_tissues, n_tpts] = pk_model.validate_pk_args(mz0, r1, k, k_trans, flips, input_function);

            % validate/make up tissue and met names
            if isempty(opts.tissue_names)
                opts.tissue_names = compose("Tissue %d", 1:n_tissues);
            elseif numel(opts.tissue_names) ~= n_tissues
                error("Unexpected number of tissue names provided");
            end

            if isempty(opts.met_names)
                opts.met_names = compose("Met %d", 1:n_mets);
            elseif numel(opts.met_names) ~= n_mets
                n_mets
                numel(opts.met_names)
                error("Unexpected number of met names names provided");
            end


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
                 met_dynamics(i_tissue,:,:) = pk_model.generate_met_dynamics(tissue_pk_params, k_trans(i_tissue));
            end

            % plotting
            % size(met_dynamics) = [tissue, met, time_pt]
            if opts.plot
                tpts = 1:n_tpts;
                figure;
                for i_tissue = 1:n_tissues
                    subplot(n_tissues, 1, i_tissue);
                    hold on;
                    for i_met = 1:n_mets
                        plot(tpts, squeeze(met_dynamics(i_tissue, i_met, :)));
                    end
                    hold off;
                    title(opts.tissue_names(i_tissue));
                    legend(opts.met_names);
                end
            end
        end


        function met_dynamics = generate_met_dynamics(params, k_trans)
            % Generates metabolite dynamics from PK parameters
            % Parameters:
            %   params          = pk parameters. size = (1,1), type = pk_params
            %   k_trans         = volumetric transfer constant. size = (1,1)
            % Outputs:
            %   met_dynamics    = metabolite dynamics. size = (met, tpt)
            arguments
                params (1,1) pk_params
                k_trans (1,1) {mustBeNumeric} = 1;
            end
        
            % prep everything for simulate_Nsite_model()
            Mz0 = params.get_Mz0();
            R1 = params.get_R1();
            k = params.get_kinetic_rates();
            flips = params.get_flips();
            TR = params.TR;
            input_function = params.InputFunction;
        
            [met_dynamics, ~] = simulate_Nsite_model(Mz0, R1, k, flips, TR, input_function .* k_trans);
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

        function met_images = apply_k_trans(met_images_low_ktrans, met_images_high_ktrans, k_trans_map)
            % Parameters:
            %   met_images_low_ktrans   = met_images where k_trans = 0 (no additional input).
            %                             size = (row, col, slice, met, time_pt)
            %   met_images_high_ktrans  = met_images where k_trans = 1 (100% additional input).
            %                             size = (row, col, slice, met, time_pt)
            %   k_trans_map             = map of k_trans values. size = (row, col, slice)
            %  Outputs:
            %   met_images              = met images with k_trans applied. size = (row, col, slice, met time_pt)

            % argument validation
            arguments
                met_images_low_ktrans (:,:,:,:,:) {mustBeNumeric}
                met_images_high_ktrans (:,:,:,:,:) {mustBeNumeric}
                k_trans_map (:,:,:) {mustBeNumeric}
            end
            
            if any(size(met_images_low_ktrans) ~= size(met_images_high_ktrans))
                error("mismatched array sizes. `met_images_low_ktrans` and `met_images_high_ktrans` must have the same size");
            end

            if any(size(k_trans_map) ~= size(met_images_low_ktrans, 1:3))
                error("mismatched array sizes. `k_trans_map` mut have the same number of rows, columns, and slices as `met_images_low_ktrans`");
            end
            
            % expand k_trans_map
            n_mets = size(met_images_low_ktrans, 4);
            n_tpts = size(met_images_low_ktrans, 5);
            full_ktrans_map = repmat(k_trans_map, 1, 1, 1, n_mets, n_tpts);

            % interpolate
            met_images = met_images_low_ktrans + (met_images_high_ktrans - met_images_low_ktrans) .* full_ktrans_map;
        end
    end


    methods (Static, Access = private)
        function [n_mets, n_tissues, n_tpts] = validate_pk_args(mz0, r1, k, k_trans, flips, input_function, t_arrival)
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
            if size(k_trans, 2) ~= n_tissues
                error("mismatched number of tissues in `Mz0` and `k_trans`");
            end
            
            n_tpts = size(flips, 2);
            if size(input_function, 2) ~= n_tpts;
                error("mismatched number of time points in `flips` and `input_function`");
            end
        end
    end
end
