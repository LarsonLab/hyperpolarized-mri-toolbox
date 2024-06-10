function [gscale, Mxy, Mz, S] = gradient_scaling_slice_compensation(rf, flips, g, fsim)
% [gscale, Mxy, Mz, S] = gradient_scaling_slice_compensation(rf, flips [, g, fsim])
%   OR
% [gscale, Mxy, Mz, S] = gradient_scaling_slice_compensation(rf, flips [, Nisodelay, fsim])
%
% In a hyperpolarized experiment, the RF pulse slice profiles can become
% substantially distorted in the transition regions due to the uneven
% depletion of the hyperpolarized magnetization.
%
% This function implements a gradient compensation scheme for equalizing
% the summed slice signal across RF pulses with arbitrary flip angles.
% This method was proposed in:
%   Deppe, M. H., Teh, K., Parra-Robles, J., Lee, K. J., and Wild, J. M.
%   Slice profile effects in 2d slice-selective mri of hyperpolarized
%   nuclei. J Magn Reson 202, 2 (Feb 2010), 180?9. 
%   Doi: 10.1016/j.jmr.2009.11.003
%
% INPUTS:
%   rf - RF pulse shape
%   flips - flip angles (radians)
%   g (optional) - gradient waveform
%       OR
%   Nisodelay (optional) - isodelay samples (e.g. is length(rf)/2 for linear
%       phase pulses, such as sinc shapes)
%   fsim (optional) - frequency offset of simulations
%
% OUTPUTS:
%   gscale - gradient scaling factors
%   Mxy, Mz - simulated slice profiles with scaling factors
%   S - simulated total slice signal with scaling factors
%
% (c) 2013-2015 The Regents of the University of California
% All Rights Reserved.
%
% Originally written by: Peder E. Z. Larson

Nrf = length(rf);
N = length(flips);

if nargin < 3 || isempty(g)
    g = [ones(Nrf,1); -round(Nrf/2)];
    rf = [rf(:);0];  Nrf = Nrf+1;
elseif length(g) == 1
    Nisodelay = g;
    g = [ones(Nrf,1); -Nisodelay];
    rf = [rf(:);0];  Nrf = Nrf+1;
end

if nargin < 4
    fsim = 0;
end

g = g * 2*pi /sum(g(:)); 

Sflips = sin(flips) .*cumprod([1 cos(flips(1:end-1))]);
    
for n = 1:N
    rf_all(1:Nrf,n) = rf /sum(rf) * flips(n); % small-tip scaling
end
gscale(1) = 1;

Nsim = 1000;
x = linspace(-2, 2, Nsim);

[a b] = abr(rf_all(1:Nrf,1),gscale(1)*g, x, fsim);

Mxy(1:Nsim,1) = 2*conj(a).*b;
Mz(1:Nsim,1) = 1 - 2*conj(b).*b;

S1 = abs(sum(Mxy))/Nsim;

maxIter = 1000;

for n = 2:N
    gscale(n) = gscale(n-1);
    glim = [gscale(n-1), 10*gscale(n-1)];
    [a b] = abr(rf_all(1:Nrf,n),gscale(n)*g, x, fsim);
    Mxy(1:Nsim,n) = 2*conj(a).*b  .* Mz(1:Nsim, n-1);
    Snew = abs(sum(Mxy(1:Nsim,n)))/Nsim;
    Niter = 0;
    while (abs(S1/Sflips(1) - Snew/Sflips(n)) > S1/Sflips(1)*1e-3) && ...
            (Niter < maxIter)
        if Snew/Sflips(n) > S1/Sflips(1)
            % shrink slice
            gold = gscale(n);
            gscale(n) = mean(glim);
            glim(1) = gold;
        else
            % slice too thin
            gold = gscale(n);
            gscale(n) = mean(glim);
            glim(2) = gold;
        end
        [a b] = abr(rf_all(1:Nrf,n),gscale(n)*g, x, fsim);
        Mxy(1:Nsim,n) = 2*conj(a).*b  .* Mz(1:Nsim, n-1);
        Snew = abs(sum(Mxy(1:Nsim,n)))/Nsim;
        Niter = Niter + 1;
    end
    Mz(1:Nsim,n) = ( 1 - 2*conj(b).*b ) .* Mz(1:Nsim,n-1);
end

figure
%subplot(211), plot(x, real(Mxy), x, imag(Mxy),'--');
[xplot yplot] = meshgrid(1:N, x);
subplot(211)
mesh(xplot, yplot,abs(Mxy)), view([-70 20])
xlabel('Pulse number'), ylabel('Slice profile'), zlabel('|M_{XY}|')
subplot(212)
mesh(xplot, yplot,abs(Mz)), view([-70 20])
xlabel('Pulse number'), ylabel('Slice profile'), zlabel('M_Z')

%pause
