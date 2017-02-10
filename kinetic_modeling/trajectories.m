function [x2, u] = trajectories( kPL, x1, Mzscale, R1P, R1L , TR )
% Compute product magnetization (e.g. lactate) using a uni-directional two-site model
% Uses substrate magnetization measurements, estimated relaxation and
% conversion rates
% all using MZ component so can account for variable flip angles

N = length(x1);

x2 = zeros(1, N); u = zeros(1,N);

A = [-R1P 0; kPL -R1L];

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

