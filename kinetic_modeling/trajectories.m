function [Mz u] = trajectories( model, params_fit, params_fixed , TR, Mzscale )

% all using MZ component so can account for variable flip angles
Nmets = size(Mzscale,1); N = size(Mzscale,2);
Mz = zeros(Nmets, N);


params_all = {'kPL', 'kPB', 'kPA', ...
	'R1P', 'R1L', 'R1A', 'R1B', ...
	'S0_L', 'S0_B', 'S0_A', ...
	'Rinj', 'Tarrival', 'Tbolus'};

nfit = 0;
for n = 1:length(params_all)
    if isfield(params_fixed, params_all(n))
        eval([params_all{n} '= params_fixed.(params_all{n});']);
    else
        nfit = nfit+1;
        eval([params_all{n} '= params_fit(nfit);']);
    end
end

% gamma input parameters
A = Tbolus;
B = Tbolus/4; % assumption that leads to reasonable looking shape

% evolution matrix
switch Nmets
    case 2
        Amat = [-R1P - kPL, 0; ...
            kPL, -R1L];
    case 3
        Amat = [-R1P - kPL - kPB, 0, 0; ...
            kPL, -R1L, 0; ...
            kPB, 0, -R1B];
    case 4
        Amat = [-R1P - kPL - kPB - kPA, 0, 0, 0; ...
            kPL, -R1L, 0, 0; ...
            kPB, 0, -R1B, 0; ...
            kPA, 0, 0, -R1A];
end

for It=1:N
    
    if strcmp(model, 'inputless')
        Mz_init = Mz(:,It-1) .* Mzscale(:, It-1);
        % FIX:
         u(It) = ( x1(t+1) - P0*exp((- R1P - kPL)*TR) ) * (R1P + kPL) / (1 - exp((- R1P - kPL)*TR));
    else
        
    t = (It-1)*TR;  % solving for Mz at time t (accounting for previous TR interval)
    
    if It == 1
    	% initial Mz in TR after RF pulse
        Mz_init = zeros(NMets,1);

        if Tarrival < 0
            % account for longer period of bolus prior to imaging
            t_preimaging = [floor(Tarrival/TR):0]*TR; % t = 0
           u(It) = sum( gampdf(t_preimaging-Tarrival,A,B)*Rinj );
        else
           u(It) = gampdf(t-Tarrival,A,B)*Rinj;
        end
        
    else
    	% initial Mz in TR after RF pulse
    	Mz_init = Mz(:,It-1) .* Mzscale(:, It-1);

        u(It) = gampdf(t-Tarrival,A,B)*Rinj;
    end
    end
        
    % solve next time point under assumption of constant input during TR
	Mz(:,It) = expm(Amat*TR)*Mz_init + u; % WRONG WRONG
	    
end

end

