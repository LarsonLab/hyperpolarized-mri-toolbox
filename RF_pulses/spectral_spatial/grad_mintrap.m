function [g] = grad_mintrap (m0, mxg, mxs, ts)
%
% [g] = grad_mintrap (m0, mxg, mxs, ts)
%
% m0 	- target zeroth moment (G/cm * s)
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
% $Header: /home/adam/cvsroot/src/ss/grad_mintrap.m,v 1.4 2013/08/15 03:34:50 adam Exp $
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
% A - trapezoid amplitude
%

dg = mxs * ts;				% Max delta in one sample
nb = ceil(sqrt (m0 / dg / ts));
A = m0/(nb*ts);
if (A <= mxg),
	na = 0;
	dg_act = A/nb;
else
	nb = ceil (mxg / dg);
	dg_act = mxg / nb;
	na = ceil((m0 - (nb^2 * dg_act * ts))/mxg/ts);
	dg_act = m0 / (nb^2 + na*nb)/ts;
	A = nb * dg_act;
end;

% Construct discrete trapezoid --- always end with a zero value
%
g = s * [[1:nb]*dg_act ones(1,na)*A [nb-1:-1:0]*dg_act];

return;


