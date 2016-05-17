function [Sscale Mzscale] = flips_scaling_factors(flips, Nt)

% # of flip angles per timepoint (e.g. phase encodes)
Nflips = size(flips,2)/Nt;

Nmets = size(flips,1);

% magnetization and signal scaling factors for each timepoint based on flip
% angles
Mzscale = zeros(Nmets,Nt); Sscale = zeros(Nmets,Nt);
for t= 1:Nt
    Iflips = [1:Nflips] + (t-1)*Nflips;
    Mzscale(:,t) = prod(cos(flips(:,Iflips)),2);
    for n = 1:Nflips
        Sscale(:,t) = Sscale(:,t) + sin(flips(:,Iflips(n))) .* prod(cos(flips(:,Iflips(1:n-1))),2);
    end
end
Sscale = Sscale/Nflips;
