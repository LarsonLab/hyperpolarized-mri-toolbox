function [gpos,gneg,g1,g2,g3] = grad_ss(m0, n, f, mxg, mxs, ts, equal)
% GRAD_SS - Calculate spectral-spatial bipolar pulse 
%   
% [gpos, gneg, g1, g2, g3] = grad_ss (m0, n, f, mxg, mxs, ts, equal)
%
% m0 	- target zeroth moment (G/cm * s) of one lobe
% n     - total number of samples to use if not []
% f     - fraction of ramp to include in bridge [0..1]    
% mxg	- maximum amplitude (G/cm)
% mxs   - maximum slew rate (G/cm/ms)
% ts    - sample time in s    
% equal - boolean if pos/neg lobes should be same
%
% gpos  - Positive lobe gradient
% gneg  - Negative lobe gradient
% g1, g2, g3 - Ramp up, bridge, ramp down gradients
%              in positive lobe
    
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
% $Header: /home/adam/cvsroot/src/ss/grad_ss.m,v 1.9 2013/08/15 03:34:50 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m0 = abs(m0);				% Must be positive

% Check n is even if equal lobes called for
%
if equal, 
    if bitget(n,1) ~= 0, 
	error('equal lobes specified, but n not even');
    end;
end;

% Convert mxs to G/cm/s
%
mxs_s = mxs * 1e3;
dg = mxs_s * ts;				% Max delta in one sample

% Determine trapezoid parameters
% 
% na - number of constant samples
% nb - number ramp samples
% nc - number samples of ramp in bridge
% A - trapezoid amplitude
%

% Do different things if number of samples is specified
%
[gp_tmp, g1_tmp, g2_tmp, g3_tmp] = grad_min_bridge(m0, f, mxg, mxs, ts);
if equal, 
    gn_tmp = -gp_tmp;
else
    m0_pos = sum(gp_tmp) * ts;
    gn_tmp = grad_mintrap(-m0_pos, mxg, mxs, ts);
end;

if isempty(n), 
    gpos = gp_tmp;
    g1 = g1_tmp;
    g2 = g2_tmp;
    g3 = g3_tmp;
    gneg = gn_tmp;
else
    if (length([gp_tmp gn_tmp]) > n)
	error('grad_ss: Solution not obtained in spec num samples');
    end;

    % Save known solution
    %
    gp_save = gp_tmp;
    g1_save = g1_tmp;
    g2_save = g2_tmp;
    g3_save = g3_tmp;
    gn_save = -gp_save;
    
    nb_save = find(diff(gp_save) == 0, 1, 'first');
    
    % Now keep decreasing number of ramp samples in 
    % positive lobe until "n" exceeded
    %
    spec_met = 1;
    while (spec_met && (nb_save > 1))
	% Get area in ramps
	%
	nb = nb_save - 1;
	nc = max(1,ceil(nb*f));
	a_ramps = (2*nb-nc+1)*nc * dg * ts;
	
	% Get number of constant samples
	%
	a_const = m0 - a_ramps;
	na = max(0,ceil(a_const/(nb*dg*ts)));
	
	% Get correct amplitude, gradients now
	%
	dg_test = m0 / ( ((2*nb-nc+1)*nc + nb*na) *ts);
	A = nb * dg_test;
	if ((A > mxg) || (dg_test > dg)), 
	    spec_met = 0;
	    continue;
	end;
	g1 = [1:(nb-nc)] * dg_test;
	g2 = [[(nb-nc+1):nb] nb*ones(1,na) [nb:-1:(nb-nc+1)]] * dg_test;
	g3 = [(nb-nc):-1:0] * dg_test;
	gp = [g1 g2 g3];
	if abs((sum(g2)*ts) - m0) > 10*eps, 
	    error('grad_ss: Area not calculated correctly');
	end;
	
	if (equal)
	    gn = -gp;
	else
	    gn = grad_mintrap(-sum(gp)*ts, mxg, mxs, ts);
	end;
	
	% See if spec still met
	%
	if length([gp gn]) < n, 
	    spec_met = 1;
	    
	    gp_save = gp;
	    g1_save = g1;
	    g2_save = g2;
	    g3_save = g3;
	    gn_save = gn;
	    nb_save = nb;
	else
	    spec_met = 0;
	end;
    end;
    
    % Fix up result to have "exactly" n samples in it!
    %
    if ~equal
	na = n - length(gn_save) - (2 * nb_save + 1);
    else
	na = (n - 2*(2 * nb_save + 1))/2;
    end;
    nb = nb_save;
    nc = max(1,ceil(nb*f));

    % Get correct amplitude, gradients now
    %
    dg_test = m0 / ( ((2*nb-nc+1)*nc + nb*na) *ts);
    A = nb * dg_test;
    if ((A >= 1.001 * mxg) || (dg_test > 1.001 * dg)), 
	error('Amp/Slew being exceeded');
    end;
    g1 = [1:(nb-nc)] * dg_test;
    g2 = [[(nb-nc+1):nb] nb*ones(1,na) [nb:-1:(nb-nc+1)]] * dg_test;
    g3 = [(nb-nc):-1:0] * dg_test;
    gpos = [g1 g2 g3];
    if abs((sum(g2)*ts) - m0) > 10*eps, 
	error('grad_ss: Area not calculated correctly');
    end;

    if (~equal)
	ratio = sum(-gn_save)/sum(gpos);
	if (ratio < 1-10*eps)  
	    %	warning('grad_ss: Improbable ratio');   % fix problem here
	end;
	gneg = gn_save / ratio;
    else
	gneg = -gpos;
    end;
end;


return;


