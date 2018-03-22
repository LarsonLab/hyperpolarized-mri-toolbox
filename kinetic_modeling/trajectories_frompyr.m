function [x2, u] = trajectories_frompyr( params_fit, x1, Mzscale, params_fixed , TR, slice_profile)
% Compute product magnetization (e.g. lactate) using a uni-directional two-site model
% Uses substrate magnetization measurements, estimated relaxation and
% conversion rates
% x1 and x2 are pyruvate and lactate, respectively, longitudinal magnetization (MZ) component estimates in order to account for variable flip angles

if nargin < 6
    slice_profile = 1;
end

N = length(x1);
M = length(slice_profile);

x2 = zeros(M, N); u = zeros(M,N);
x1sliced = zeros(M,N);
x1_profile = zeros(M,N);

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

% initial lactate magnetization
x2(:,1) = L0_start;

% for estimating input at various slice locations
x1_profile(:,1) = slice_profile;% ./ mean(slice_profile);
% first time point slice contributions
x1sliced(:,1) = x1(1) * x1_profile(:,1) ;

% MOVE SLICE PROFILE into Mzscale & Sscale...

    for t=1:N-1

        for m= 1:M

            % scaling to next time point
        P0 = x1(t) * (Mzscale(1, t, m));
        L0 = x2(m,t)* (Mzscale(2, t, m));
        
        % estimate input, assuming this is constant during TR interval
        u(m,t) = ( x1(t+1) - P0*exp((- R1P - kPL)*TR) ) * (R1P + kPL) / (1 - exp((- R1P - kPL)*TR));
        
        % solve next time point under assumption of constant input during TR
        x2(m,t+1) = exp(-R1L*TR)*((L0*R1L*R1P - L0*R1L^2 - kPL*u(m,t) + L0*R1L*kPL + P0*R1L*kPL)/(R1L*(R1P - R1L + kPL)) + (kPL*u(m,t)*exp(R1L*TR))/(R1L*(R1P - R1L + kPL))) - exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u(m,t) + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u(m,t)*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)));
        
        % solution for next source (pyruvate) time point if needed:
        % x1sliced(m,t+1) = (exp(-TR*(R1P + kPL))*((kPL*(P0*R1P - u(m,t) + P0*kPL))/((R1P + kPL)*(R1P - R1L + kPL)) + (kPL*u(m,t)*exp(R1P*TR + kPL*TR))/((R1P + kPL)*(R1P - R1L + kPL)))*(R1P - R1L + kPL))/kPL;
        
        
        % % Solution without an input function:
        %  x2(t+1) = (-(kPL*exp((- R1P - kPL)*TR) - kPL*exp(-R1L*TR))/(R1P - R1L + kPL))*Mzscale(1, t)*x1(t) ...
        %      +  exp(-R1L*TR)*Mzscale(2, t)*x2(t);
        
        end
        
        % estimate evolution of how pyruvate is distributed across slice profile
        % x1_profile(:,t+1) = ( x1sliced(:,t+1) / x1(t+1)); 
    
end

% use averaged signal across slice_profile
u = mean(u,1);
x2 = mean(x2,1);

end

