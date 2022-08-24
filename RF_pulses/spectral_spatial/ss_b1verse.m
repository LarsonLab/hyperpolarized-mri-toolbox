function [rfv, gv] = ss_b1verse(g, rf, b1max, gmax, smax, ts, gamma, slew_penalty, dbg);
% [rfv, gv] = ss_b1verse(g, rf, b1max, gmax, smax, ts, gamma, slew_penalty, dbg)
%
% Use VERSE to reduce maximum RF amplitude while retaining an
% identical pulse duration.
% Input RF and gradient must correspond (ie RF should be once VERSE'd if
% gradient is varying).
% A slew rate penalty that reduces the allowable slew rate for larger RF
% amplitudes is also implemented with the slew_penalty variable that
% determines the degree of this penalty.
%
% INPUTS:
%   g - gradient in G/cm
%   rf - scaled as per rftools (sum(rf) = flip)
%   b1max - max RF
%   gmax - max gradient, in G/cm
%   smax - max slew rate, in G/cm/ms
%   ts - sampling interval, in sec
%   gamma - gyromagnetic ratio, in Hz/G
%   slew_penalty (optional) - degree of slew rate reduction
%       default is 0 (no penalty)
%   dbg (optional) - indicates debug level
%
% OUTPUTS:
%   rfv - versed RF
%   gv - new gradient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2011 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:
% - allow negative gradients

if nargin < 8
    slew_penalty = 0;
end

if nargin < 9
  dbg = 0;
end

if length(g) == 1
  g = g*ones(size(rf));
end

if length(b1max) == 1
  b1max = b1max*ones(size(rf));
end

dt = ts*ones(size(rf));
t_unif = ([0:length(rf)-1] + 0.5) * ts;
dt_unif = dt;

S_rfg = 1 / (2*pi*gamma*ts);

% convert to Gauss
if dbg >= 2
  rfg = rf *S_rfg;
  fprintf('Initial Max RF fraction: %f\n', max(abs(rfg) ./ b1max));
end

T = sum(dt);

notdone = 1;

% maxiter may not be necessary because
% the routine should either fail or succeed at some point
niter = 0;
maxiter = 100;

%area = sum(g);

while (niter < maxiter)
    niter = niter + 1;
    
  Imaxed = zeros(size(rf));
  Islewed = zeros(size(rf));
  % fix points with zero gradient
  Islewed(find(g==0)) = 1;
  % fix endpoints
  Islewed([1,end]) = 1;

  if dbg >= 2
    figure(99)
    subplot(311)
    plot(cumsum(dt),g)
    subplot(312)
    plot(cumsum(dt),abs(rf) / (2*pi*gamma*ts))
    subplot(313)
    plot(cumsum(dt),diff(g) ./ ((dt(1:end-1) + dt(2:end))/2))
    pause
  end

  b = ones(1,5)/5; 
  b =firls(8, [0 .03 .06 1], [1 1 0 0]);
  rffilt = filtfilt(b, 1, abs(rf)) / (sum(b)^2);
  rffilt = max(abs(rf), rffilt);
%    rffilt = abs(rf);
  filtmax = max( max(abs(rffilt)), max(b1max)/S_rfg);
  s_mod = (1 - 0.999*rffilt/filtmax) .^slew_penalty;
%  s_mod = (1 - 0.999*abs(rf)/max(abs(rf))) .^slew_penalty;

%    s_mod = (tanh((slew_penalty(1) - rffilt/filtmax) * slew_penalty(2)) + 1 ) / 2;
  smax_mod = smax*1e3 * s_mod;
  
  % Go through gradient forward looking for slew-rate violations
  for k = 2:length(rf)
    slew = (g(k) - g(k-1)) / (dt(k) + dt(k-1)) / 0.5;

    if 1
        smaxk = smax*1e3 * ...
            (1 - (abs(rffilt(k)) + abs(rffilt(k-1)))*S_rfg / (b1max(k)+b1max(k-1))  )^slew_penalty;
    else
        smaxk = smax_mod(k);
    end
    if (slew > smaxk )
      gh = g(k); gl = g(k-1);
      dth = dt(k); dtl = dt(k-1);
        
      a = smax_mod(k) * dth;
      b = smax_mod(k) * dtl + 2*gl;
      c = -2*gh;
    
      scale = (-b + sqrt(b^2-4*a*c))/(2*a);
    
      rf(k) = rf(k) / scale;
      g(k) = g(k) / scale;
      dt(k) = dt(k) * scale;
      Islewed(k) = 1;
    end
  
  end

  % Go through gradient in reverse looking for slew-rate violations
  for k = length(rf)-1:-1:1
    slew = (g(k) - g(k+1)) / (dt(k) + dt(k+1)) / 0.5;
  
    if 1
        smaxk = smax*1e3 * ...
        (1 - (abs(rffilt(k)) + abs(rffilt(k+1)))*S_rfg / (b1max(k) + b1max(k+1)) )^slew_penalty;
    else
        smaxk = smax_mod(k);
    end
    if (slew > smaxk)
      gh = g(k); gl = g(k+1);
      dth = dt(k); dtl = dt(k+1);
      
      a = smax_mod(k) * dth;
      b = smax_mod(k) * dtl + 2*gl;
      c = -2*gh;
    
      scale = (-b + sqrt(b^2-4*a*c))/(2*a);
    
      rf(k) = rf(k) / scale;
      g(k) = g(k) / scale;
      dt(k) = dt(k) * scale;
      Islewed(k) = 1;
    end
  
  end

    
  % Check for violoations of max RF
  for k = 1:length(rf)
    rfscale = abs(rf(k)) * S_rfg / (b1max(k)*.999);
        
    if (rfscale > 1)
      rf(k) = rf(k) / rfscale;
      g(k) = g(k) / rfscale;
      dt(k) = dt(k) * rfscale;
    
      Imaxed(k) = 1;
    end
  end

  
  % compensate for time expansion in removing max RF and slew
  % rate violations by shrinking other time samples and increasing
  % RF and gradient correspondingly
  Ifixed = find(Imaxed | Islewed);
  Ishrink = find(~(Imaxed | Islewed));
  shrink_scale = (T - sum(dt(Ifixed))) / sum(dt(Ishrink));

  % Check if pulse exceeds length requirement
  if shrink_scale < 0
    % error('ss_verse: RF pulse unrealizeable');
    rfv = [];
    gv = [];
    return;
  elseif (abs(shrink_scale-1) < 1e-14) || (isempty(Imaxed) && isempty(Islewed))
    % DONE!
    notdone = 0;
  else
    rf(Ishrink) = rf(Ishrink) / shrink_scale;
    g(Ishrink) = g(Ishrink) / shrink_scale ;
    dt(Ishrink) = dt(Ishrink) * shrink_scale;

    % timepoints at center of little hard pulses
    t = cumsum([0, dt(1:end-1)]) + dt/2;
    
    rf = interp1(t, rf, t_unif, 'spline', 'extrap');
    g = interp1(t, g, t_unif, 'spline', 'extrap');
%    g = g/sum(g) * area;
    dt = dt_unif;  
  end

  niter = niter + 1;

end

% check if gradient limits exceeded
gscale = abs(g) / gmax;
if max(gscale) > 1
  %error('ss_verse: RF pulse unrealizeable');
  rfv = [];
  gv = [];
  return;
end

rfv = rf;
gv = g;

if dbg >= 2
  rfvg = rfv*S_rfg;
  fprintf('Final Max RF fraction: %f\n', max(abs(rfvg) ./ b1max));
end
