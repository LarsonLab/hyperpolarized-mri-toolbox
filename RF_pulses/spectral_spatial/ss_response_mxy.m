function [mxy mz] = ss_response_mxy(g,rf,z,f,ts,gamma, isodelay)
% SS_RESPONSE_MXY - Get Mxy response at given z, f position/frequency
%
% [mxy mz] = ss_response_mxy(g,rf,z,f,ts,gamma)
%
% Input
%    g - gradient, G/cm
%    rf - rf G
%    z - cm
%    f - frequency
%    ts - sampling time (s)
%    gamma - gyromagnetic ratio
%    [isodelay] - Unwind spectral phase shift for given isodelay - default: 0
    
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
% $Header: /home/adam/cvsroot/src/ss/ss_response_mxy.m,v 1.6 2013/08/15 03:34:50 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ( (nargin < 7) || isempty(isodelay) ), 
	isodelay = 0;
    end;
    
    % Convert RF to a rotation in radians
    %
    rf_rot = 2 * pi * gamma * rf(:) * ts;

    % Build gradient that gives rotation  
    % in radians when scaled by "f" (off-resonance in Hz) 
    %
    gf_rot = 2 * pi * ts * ones(size(g(:)));
    
    % Convert gradient to a rotation in radians
    % when scaled by "z"
    %
    gz_rot = 2 * pi * gamma * g(:) * ts;
    
    % Add single samples at end of all waveforms to account for isodelay correction
    %
    rf_rot = [rf_rot;0];
    gf_rot = [gf_rot;(-2 * pi * isodelay)];
    gz_rot = [gz_rot;0];

    
    % Get Mxy now
    %
    [a b] = abr(rf_rot, gz_rot+i*gf_rot, z, f);
    mxy = ab2ex(a, b);
    mz = ab2sat(a, b);

    
    
    
