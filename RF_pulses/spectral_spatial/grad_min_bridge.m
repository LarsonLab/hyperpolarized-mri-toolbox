function [g, g1, g2, g3] = grad_min_bridge(m0, f, mxg, mxs, ts)
% GRAD_MIN_BRIDGE - Determine gradient trapezoid that gives required area from middle "bridge" section
%
% [g, g1, g2, g3] = grad_min_bridge (m0, f, mxg, mxs, ts)
%
% m0 	- target zeroth moment (G/cm * s)
% f  - fraction of ramp to include in bridge [0..1]    
% mxg	- maximum amplitude (G/cm)
% mxs   - maximum slew rate (G/cm/ms)
% ts    - sample time in s    
%
% g     - Gradient
%

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
%
% $Header: /home/adam/cvsroot/src/ss/grad_min_bridge.m,v 1.5 2013/08/15 03:34:50 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
if (m0 < 0)
	s = -1;
	m0 = -m0;
else
	s = 1;
end;

% Convert mxs to G/cm/s
%
mxs = mxs * 1e3;

% Determine trapezoid parameters
% 
% na - number of constant samples
% nb - number ramp samples
% nc - number samples of ramp in bridge
% A - trapezoid amplitude
%

dg = mxs * ts;				% Max delta in one sample

% Assume triangle at first and see if max amp requirements met
% -- quadratic with aq, bq, cq coeffs
%
if (f ~= 0)
    aq = (2-f)*f;
    bq = f;
    cq = -m0/(dg*ts);
    nb = (-bq + sqrt(bq^2 - 4*aq*cq))/(2*aq);

    nb = ceil(nb);
    nc = max(1,ceil(nb*f));
else
    A = m0 / (2*ts);
    nb = ceil(A / dg);
    nc = 1;
end;

% Test result
%
dg_test = m0 / ((2*nb-nc+1)*nc*ts);
A = nb * dg_test;
if (A <= mxg) && (dg_test < dg),				% This works!
    g1 = s*[1:(nb-nc)] * dg_test;
    g2 = s*[[(nb-nc+1):nb] [nb:-1:(nb-nc+1)]] * dg_test;
    g3 = s*[(nb-nc):-1:0] * dg_test;
    g = [g1 g2 g3];
    if abs((sum(g2)*ts) - s*m0) > 10*eps, 
	fprintf(1,'Area Spec: %f Actual: %f\n', m0, sum(g)*ts);
	error('grad_min_bridge: Area not calculated correctly');
    end;
else  %% Must be trapezoid
    % Subtract area of ramps
    %
    nb = ceil(mxg/dg);
    nc = max(1,ceil(nb*f));
    dg_test = mxg/nb;
    a_ramps = (2*nb-nc+1)*nc * dg_test * ts;
    
    % get number of const samples
    %
    a_const = m0 - a_ramps;
    na = ceil(a_const/ts/mxg);
    
    % Get correct amplitude now
    %
    dg_test = m0 / ( ((2*nb-nc+1)*nc + nb*na) *ts);
    A = nb * dg_test;
    if ((A > mxg) || (dg_test > dg)), 
	error('Amp/Slew being exceeded');
    end;
    g1 = s* [1:(nb-nc)] * dg_test;
    g2 = s* [[(nb-nc+1):nb] nb*ones(1,na) [nb:-1:(nb-nc+1)]] * dg_test;
    g3 = s* [(nb-nc):-1:0] * dg_test;
    g = [g1 g2 g3];
    if abs((sum(g2)*ts) - s*m0) > 10*eps, 
	error('grad_min_bridge: Area not calculated correctly');
    end;
    
end;

return;


