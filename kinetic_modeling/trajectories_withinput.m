function [x1, x2] = trajectories_withinput( params_fit, params_fixed , TR, N, Mzscale )

% all using MZ component so can account for variable flip angles

x1 = zeros(1, N); x2 = zeros(1,N);

params_all = {'kPL', 'R1L', 'R1P', 'Rinj', 'Tarrival', 'Tbolus'};
nfit = 0;
for n = 1:length(params_all)
    if isfield(params_fixed, params_all(n))
        eval([params_all{n} '= params_fixed.(params_all{n});']);
    else
        nfit = nfit+1;
        eval([params_all{n} '= params_fit(nfit);']);
    end
end

for It=1:N
    
    
    t = (It-1)*TR;  % solving for Mz at time t (accounting for previous TR interval
    if t <= Tarrival
        % no arrival during TR
        x1(It) = 0;
        continue
    elseif t-TR >= (Tarrival+Tbolus)
        % after bolus arrival completed
        u = 0;
    elseif t-TR < Tarrival && t > Tarrival
        % TR interval contains arrival time
        u = Rinj * (t - Tarrival)/TR;
    elseif t-TR < (Tarrival+Tbolus) && t > (Tarrival+Tbolus)
        % TR interval contains end of bolus
        u=Rinj * (Tarrival+Tbolus - (t-TR))/TR;
    else
        % TR interval should be completely in bolus
        u = Rinj;
    end
    
    if It == 1
        % assumes experiment starts before signal appears...
        P0 = 0;
        L0 = 0;
    else
        P0 = x1(It-1)*Mzscale(1, It-1);
        L0 = x2(It-1)*Mzscale(2, It-1);
    end
    
    % solve next time point under assumption of constant input during TR
    x2(It) = exp(-R1L*TR)*((L0*R1L*R1P - L0*R1L^2 - kPL*u + L0*R1L*kPL + P0*R1L*kPL)/(R1L*(R1P - R1L + kPL)) + (kPL*u*exp(R1L*TR))/(R1L*(R1P - R1L + kPL))) - exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)));
    
    % solution for next source (pyruvate) time point if needed:
    x1(It) = (exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)))*(R1P - R1L + kPL))/kPL;
        
end

end

