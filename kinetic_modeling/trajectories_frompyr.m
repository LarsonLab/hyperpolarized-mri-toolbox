function [x2, u] = trajectories_frompyr( params_fit, x1, Mzscale, params_fixed , TR )
% Compute product magnetization (e.g. lactate) using a uni-directional two-site model
% Uses substrate magnetization measurements, estimated relaxation and
% conversion rates
% x1 and x2 are pyruvate and lactate, respectively, longitudinal magnetization (MZ) component estimates in order to account for variable flip angles

N = length(x1);

x2 = zeros(1, N); u = zeros(1,N);

params_all = {'kPL', 'R1L', 'R1P', 'L0_start'};
nfit = 0;
for n = 1:length(params_all)
    if isfield(params_fixed, params_all(n))
        eval([params_all{n} '= params_fixed.(params_all{n});']);
    else
        nfit = nfit+1;
        eval([params_all{n} '= params_fit(nfit);']);
    end
end

x2(1) = L0_start;

for t=1:N-1
    
    P0 = x1(t)*Mzscale(1, t);
    L0 = x2(t)*Mzscale(2, t);
    
    % estimate input, assuming this is constant during TR interval
    u(t) = ( x1(t+1) - P0*exp((- R1P - kPL)*TR) ) * (R1P + kPL) / (1 - exp((- R1P - kPL)*TR));
    
    % solve next time point under assumption of constant input during TR
    x2(t+1) = exp(-R1L*TR)*((L0*R1L*R1P - L0*R1L^2 - kPL*u(t) + L0*R1L*kPL + P0*R1L*kPL)/(R1L*(R1P - R1L + kPL)) + (kPL*u(t)*exp(R1L*TR))/(R1L*(R1P - R1L + kPL))) - exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u(t) + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u(t)*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)));
    
    % % solution for next source (pyruvate) time point if needed:
    % x1(t+1) = (exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u(t) + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u(t)*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)))*(R1P - R1L + kPL))/kPL;
    
    % % Solution without an input function:
    %  x2(t+1) = (-(kPL*exp((- R1P - kPL)*TR) - kPL*exp(-R1L*TR))/(R1P - R1L + kPL))*Mzscale(1, t)*x1(t) ...
    %      +  exp(-R1L*TR)*Mzscale(2, t)*x2(t);
    
end

end

