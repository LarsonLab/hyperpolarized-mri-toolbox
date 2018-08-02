function [Mz_product, u] = trajectories_frompyr_reverse( params_fit, Mz_pyr, Mzscale, params_fixed , TR )
% Compute product magnetization (e.g. lactate) using a uni-directional two-site model
% Uses substrate magnetization measurements, estimated relaxation and
% conversion rates
% x1 and x2 are pyruvate and lactate, respectively, longitudinal magnetization (MZ) component estimates in order to account for variable flip angles

N = length(Mz_pyr);
Mz_all = zeros(2, N);
u = zeros(1,N);

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

Mz_all(1,:) = Mz_pyr;
Mz_all(2,N) = L0_start;

A = [-R1P-kPL, 0
     +kPL, -R1L];
Ap = -R1P-kPL;
 
for It=N:-1:2

Mz_init = Mz_all(:,It);% .* Mzscale(:, It);

% estimate input, assuming this is constant during TR interval
% This calculation could be improved for noise stability?
u(It) = ( Mz_pyr(It-1)*Mzscale(1,It-1) - Mz_init(1)*exp(Ap*-TR) ) * Ap / (exp(Ap*-TR) - 1);

xstar = - inv(A)*[u(It),0].';

% solve next time point under assumption of constant input during TR
Mz_plus = xstar + expm(A*-TR) * (Mz_init - xstar);


Mz_all(:,It-1) = Mz_plus ./ Mzscale(:, It-1);


end

Mz_product = Mz_all(2,:);

end

