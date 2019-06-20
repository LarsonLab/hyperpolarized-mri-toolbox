function [Mxy_ev, Mz_ev, Mxy_iv, Mz_iv] = simulate_Nsite_perfused_voxel_model(Mz0, R1, k, kve, vb, flips, TR, VIF)
% [Mxy, Mz] = simulate_Nsite_model(Mz0, R1, k, flips, TR, [input_function])
%
% Simulates the magnetization evolution in a N-site exchange model with
% single substrate and multiple products for varying flip angles.
% Allows for up to N=4
% Applied disretely as flips(:,1), TR, flips(:,2), TR, etc...
%
% INPUTS:
%   Mz0 [1xN] -  model starts from this  magnetization distribution
%       (If a single value is specified for Mz0, then this value is a time (s) for which the model will evolve prior to playing flips.  Magnetization is assumed to start all in substrate (e.g. pyruvate)
%   R1 - relaxation times for N sites (1/s)
%   k [Nx2] - conversion rates between sites (1/s), first column is forward
%   rate, second column is reverse rate
%   flips - flip angles (radians)
%   TR - repetition time (s)
%   input_function (optional) - additional input function to add to
%   pyruvate for all time points.
%
% OUTPUTS:
%   Mxy,Mz - expected magnetization evolution immediately after flips

Nmets = size(flips,1);
N = size(flips,2);

if isempty(R1)
    R1 = zeros(1,Nmets);  % infinite relaxation time
end

if length(R1) == 1
    R1 = R1*ones(1,Nmets);
end

ve = 1-vb;
kvedve = kve/ve;
switch Nmets
    case 2
        A = [-R1(1)-k(1,1)-kvedve +k(1,2)
            +k(1,1) -R1(2)-k(1,2)-kvedve];
    case 3
        A = [-R1(1)-k(1,1)-k(2,1)-kvedve +k(1,2) +k(2,2)
            +k(1,1) -R1(2)-k(1,2)-kvedve 0
            +k(2,1) 0 -R1(3)-k(2,2)-kvedve];
    case 4
        A = [-R1(1)-k(1,1)-k(2,1)-kvedve +k(1,2) +k(2,2) +k(3,2)
            +k(1,1) -R1(2)-k(1,2)-kvedve 0 0
            +k(2,1) 0 -R1(3)-k(2,2)-kvedve 0
            +k(3,1) 0 0 -R1(4)-k(3,2)-kvedve];
end

% Diagonlize to permit matrix integral.
[P,D]=eig(A);
% D is diagonlized matrix (2x2)
% P is matrix of eigenvectors
%   --> A*P=P*D
%   --> A=P*D\P
dD=diag(D); % get diagonal values (1x2)

Mxy_ev(1:Nmets,1) = Mz0(:) .* sin(flips(:,1));
Mz_ev(1:Nmets,1) = Mz0(:) .* cos(flips(:,1));
Mxy_iv = VIF .* sin(flips);
Mz_iv = VIF;

%Calculate signal evolution each TR:

for n = 2:N
    %First account for EV signal already present and its evolution
    Mz_ev_fromev = (exp(dD*TR)).*(P\(Mz_ev(:,n-1)));
    %
    %Now calculate new spins flowing into the system
    %
    %Assume piecewise linear VIF. Diagonalize:
    dff1=P\VIF(:,n-1);    % Diag'd VIF @ start of TR
    dff2=P\VIF(:,n);  % Diag'd VIF @ end of TR
    %Get slope and y-intercept for diagonolized forcing function:
    b=dff1;
    m=(dff2-dff1)/TR;
    %At the end of this TR, inflowing spins will lead to:
    Mz_ev_inflow1=exp(dD*TR).*(-b./dD).*(exp(-dD*TR)-1);
    Mz_ev_inflow2=exp(dD*TR).*m.*(((-TR./dD)-(1./dD./dD)).*exp(-dD*TR)+(1./dD./dD));
    
    %Total signal at end of TR equals IC for next TR:
    Mz_ev_m = P*(Mz_ev_fromev + kvedve*(Mz_ev_inflow1+Mz_ev_inflow2));

    Mxy_ev(:,n) =  Mz_ev_m .* sin(flips(:,n));
    Mz_ev(:,n) = Mz_ev_m .* cos(flips(:,n));
end