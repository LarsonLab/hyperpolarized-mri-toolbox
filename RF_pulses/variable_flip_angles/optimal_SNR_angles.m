function [flips, Mxy, Mz] =optimal_SNR_angles(T,TR,k,R1,Tinj,substrate_flips,dbg);
% flips =optimal_SNR_angles(T1,TR,T,tend,substrate_flips);
%
% Calculates series of flip angles to maximize cumulative signal of
% metabolic product (e.g. lactate) in hyperpolarized MR, for a fixed set of
% defined substrate (e.g. pyruvate) flip angles.
% Uses a gamma variate input model
%
% It can be derived based on this code that, to a very good approximation,
% only the relaxation rate of the product, the total scan time, and the TR
% influence the resulting flip angles
%
%
% INPUTS:
%   T - total scan time in sec
%   TR - repetition time in sec
%   k - conversion rate from product to substrate, 1/sec
%   R1 - relaxation rate of substrate and product (can specify both or one if equal)
%   Tinj - approximate duration of injection
%   substrate_flips (optional) - flip angles of substrate
%   dbg (optional) - print out messages and create plots
%
% OUTPUTS:
%   flips - flip angles for substrate and product
%   Mxy,Mz - simulated magnetizations using optimal scheme and input
%      parameters
%
%
% (c) 2014-2015 The Regents of the University of California
% All Rights Reserved.
%
% Authors: John Maidens, Sonam Machingal

if length(R1) == 1
    R1 = [R1 R1];
end

if length(k) == 1
    % only kPL provided, estimate kTRANS
    k = [k, .03];
end


% define number of acquisitions-
N = max(floor(T/TR),1);

if nargin<6 || isempty(substrate_flips)
    % generate substrate flip angles
    substrate_flips = vfa_const_amp(N,pi/2,exp(-TR*(R1(1) + k(1))));
end

if nargin < 7
    dbg = 0;  % don't print outputs
end

% define model parameters

syms R1P R1L kPL kTRANS t0  A0 alpha_1 beta_1
parameters_of_interest = [kPL kTRANS];
parameters_of_interest_nominal_values = [k(1) k(2)];
nuisance_parameters = [alpha_1 beta_1 A0];
nuisance_parameters_nominal_values = [Tinj/3 Tinj/4 1];  % approximate incorporation of Tinj here, set as known_parameter?
known_parameters = [R1P R1L t0];
known_parameter_values = [R1(1) R1(2) 0];

% define system matrices for differential eq. dx/dt = A*x(t) + B*u(t)
A = [ -kPL-R1P  0   ;
    kPL      -R1L];
B = [kTRANS; 0];

% define input function shape

u = @(t) A0*(t-t0)^alpha_1*exp(-(t-t0)/beta_1);      %gamma variate


% generate SNR-optimal lactate flip angles

[flips, Mxy, Mz, u_sim] = optimal_flip_angle_design_max_product_SNR(A, B, u, TR, N, ...
    parameters_of_interest, nuisance_parameters, ...
    known_parameters, parameters_of_interest_nominal_values, ...
    nuisance_parameters_nominal_values, known_parameter_values, ...
    substrate_flips(:),dbg);

if dbg
    display('=====Metabolite flip angle generation complete=====')
    
    % plot optimal flip angles
    t = [0:N-1]*TR;
    figure
    subplot(3,2,1)
    plot(t,180/pi*flips(2,:), 'mx-')
    title('flip angle scheme for Lactate')
    xlabel('time')
    ylabel('flip angle (degrees)')
    %axis([1 N 0 100])
    
    subplot(3,2,2)
    plot(t,180/pi*flips(1,:), 'mx-')
    title('flip angle scheme for Pyruvate')
    xlabel('time')
    ylabel('flip angle (degrees)')
    %axis([1 N 0 100])
    
    subplot(3,2,3)
    plot(t,Mz(2,:), 'mx-')
    title('Simulated magnetization for Lactate')
    xlabel('time')
    ylabel('signal')
    %axis([1 N 0 max(Mz(:))])
    
    subplot(3,2,4)
    plot(t,Mz(1,:), 'mx-')
    title('Simulated magnetization for pyruvate')
    xlabel('time')
    ylabel('signal')
    %axis([1 N 0 max(Mz(:))])
    
    subplot(3,2,5)
    plot(t,Mxy(2,:), 'mx-')
    title('Simulated signal for Lactate')
    xlabel('time')
    ylabel('signal')
    %axis([1 N 0 max(Mz(:))])
    
    subplot(3,2,6)
    plot(t,Mxy(1,:), 'mx-')
    title('Simulated signal for pyruvate')
    xlabel('time')
    ylabel('signal')
    
    %figure
    %plot(t, u_sim)
end
end


function [thetas_opt, Mxy, Mz, u] = optimal_flip_angle_design_max_product_SNR(A, B, u, TR, N, ...
    parameters_of_interest, nuisance_parameters, ...
    known_parameters, parameters_of_interest_nominal_values, ...
    nuisance_parameters_nominal_values, known_parameter_values, ...
    substrate_flips,dbg)

if dbg
    display('===== Computing optimal flip angles =====')
end

%  Ensure that input arguments are compatible
[i, j] = size(A);
if i ~= j
    error('The argument A must be a square matrix');
end

if size(B, 1) ~= i
    error('The number of rows of B must equal the size of A');
end

if size(B, 2) ~= size(u, 1)
    error('The number of columns of B must equal the number of entries of u')
end

if size(parameters_of_interest) ~= size(parameters_of_interest_nominal_values)
    error('The arguments parameters_of_interest and parameters_of_interest_nominal_values must be the same size');
end

if size(nuisance_parameters) ~= size(nuisance_parameters_nominal_values)
    error('The arguments nuisance_parameters and nuisance_parameters_nominal_values must be the same size');
end

if size(known_parameter_values) ~= size(known_parameters)
    error('The arguments known_parameters and known_parameter_values must be the same size');
end



% compute symbolic system discretization
Ad_sym = expm(TR*A);
Bd_sym = A\((Ad_sym - eye(size(A)))*B);

% evaluate system matrices and input at nominal parameter values
Ad_nom = double(subs(Ad_sym, [parameters_of_interest, nuisance_parameters, known_parameters], ...
    [parameters_of_interest_nominal_values, nuisance_parameters_nominal_values, known_parameter_values]));
Bd_nom = double(subs(Bd_sym, [parameters_of_interest, nuisance_parameters, known_parameters], ...
    [parameters_of_interest_nominal_values, nuisance_parameters_nominal_values, known_parameter_values]));
u_fun = matlabFunction(subs(u,[parameters_of_interest, nuisance_parameters, known_parameters], ...
    [parameters_of_interest_nominal_values, nuisance_parameters_nominal_values, known_parameter_values]));

% define objective function for total product SNR
subindex=@(A) A(2,:); %select second row
obj = @(product_flips) -(sum(sum(subindex(trajectories([substrate_flips, product_flips, pi/2*ones(N,1)], Ad_nom, Bd_nom, u_fun, TR, N)))));

% initialize optimization problem
if dbg
    options = optimset('MaxFunEvals', 50000, 'MaxIter', 1000, 'Display', 'iter');
else
    options = optimset('MaxFunEvals', 50000, 'MaxIter', 1000, 'Display', 'off', 'LargeScale', 'off');
end

initial_product_flips=zeros(N,1); % or vfa_const_amp?
%%

% perform optimization
product_flips_opt = fminunc(obj, initial_product_flips, options);

% add fixed flip angles to result
thetas_opt = [substrate_flips, product_flips_opt].';

[y, Mz, u] = trajectories([thetas_opt.', pi/2*ones(N,1)], Ad_nom, Bd_nom, u_fun, TR, N);
Mxy = y(1:2,:);
end


function [y, x, u] = trajectories(thetas, Ad, Bd, u_fun, TR, N)

n = size(Ad, 1);
x = zeros(n, N);

% compute input and state trajectories
u =  u_fun(TR*(1:N));

x(:, 1) = Bd*u(1);
for k=1:N-1
    x(:, k+1) = Ad * diag(cos(thetas(k,1:n))) * x(:, k) + Bd * u(k+1);
end

% concatenate state and input trajectories
y = sin(thetas').*[x; u];

end



