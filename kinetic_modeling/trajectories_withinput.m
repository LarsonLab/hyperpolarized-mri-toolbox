function [x1, x2] = trajectories_withinput( params_fit, params_fixed , TR, N, Mzscale )

% all using MZ component so can account for variable flip angles

x1 = zeros(1, N); x2 = zeros(1,N);

params_all = {'kPL', 'R1L', 'R1P', 'Rinj', 'Tarrival', 'Tend'};
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
    
        if It == 1
            P0 = 0;
            L0 = 0;
        else
        
        
    P0 = x1(It-1)*Mzscale(1, It-1);
    L0 = x2(It-1)*Mzscale(2, It-1);
        end

    t = (It-1)*TR;
    if t < Tarrival
        x1(It) = 0; x2(It) = 0;
        continue
    elseif t < Tend
        u = Rinj;
    else
        u=0;
    end
        
    % solve next time point under assumption of constant input during TR
    x2(It) = exp(-R1L*TR)*((L0*R1L*R1P - L0*R1L^2 - kPL*u + L0*R1L*kPL + P0*R1L*kPL)/(R1L*(R1P - R1L + kPL)) + (kPL*u*exp(R1L*TR))/(R1L*(R1P - R1L + kPL))) - exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)));
    
    % solution for next source (pyruvate) time point if needed:
    x1(It) = (exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)))*(R1P - R1L + kPL))/kPL;
    
    % % Solution without an input function:
    %  x2(t+1) = (-(kPL*exp((- R1P - kPL)*TR) - kPL*exp(-R1L*TR))/(R1P - R1L + kPL))*Mzscale(1, t)*x1(t) ...
    %      +  exp(-R1L*TR)*Mzscale(2, t)*x2(t);
    
end

end

