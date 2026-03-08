classdef pk_model
    methods (Static)
        function met_dynamics = generate_met_dynamics(params)
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

        function met_images = generate_met_images(str, met_dynamics)
            % met_dynamics = numeric (tissue, met, time_pt)
            % str = structure
            % outputs:
            %   met_images = (row, col, slice, tissue, met, time_pt)
            arguments
                str (1,1) structure
                met_dynamics {mustBeNumeric}
            end
            
            % validate
            if size(str.Mask, 4) ~= size(met_dynamics,1)
                error("mismatched number of tissues for structure and met dynamics");
            end
        
            time_pts = size(met_dynamics, 3);
            n_mets = size(met_dynamics, 2);
            
            % reshape + expand met dynamic
            met_dynamics = reshape(met_dynamics, cat(2,[1,1,1],size(met_dynamics)));
            met_dynamics = repmat(met_dynamics, cat(2, size(str.Mask, 1:3), [1,1,1]));
            
            % expand mask
            mask = str.Mask; % = (row, col, slice, tissue)
            mask = repmat(mask, cat(2,[1,1,1,1,n_mets,time_pts]));
            
            % multiply
            met_images = met_dynamics .* mask;
        end
    end
end
