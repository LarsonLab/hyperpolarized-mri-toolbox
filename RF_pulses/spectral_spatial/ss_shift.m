function rf_shift = ...
    ss_shift(g, rf, z_shift, f_shift);
    
					
% SS_SHIFT - frequency and/or spatial shift of spectral-spatial pulse
%   
% [g, rf, fs_best, z_plot, f_plot, m_plot] = ...
%     ss_design(z_thk, z_tb, z_de, f, a_angs, de,...
%               ptype, z_ftype, s_ftype, ss_type, ...
%               f_off, dbg, no_plot)
%
% INPUTS
%   g - gradient (G/cm)
%   rf - RF (G)
%   z_shift - slice shift [cm]
%   f_shift - frequency shift [Hz]
%
% OUTPUTS
%   rf_shift - updated RF (G)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Peder E. Z. Larson
%
% (c)2024 The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check all inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    if (nargin < 3)
	error(['Usage: ss_shift(z_thk, z_tb, z_d, f, a_angs, d,' ...
	       'ptype, z_fttype, s_ftype, ss_type, foff)']); 
    end
    


    % Check x_shift
    %
    if isempty(z_shift)
	    z_shift = 0;
    end

    % Check f_shift
    %
    if (nargin < 4) || isempty(f_shift) 
    	f_shift = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize globals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ss_globals;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % frequency shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rf_shift = rf .* exp(-1i*2*pi*[0:length(rf)-1]*SS_TS*f_shift);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spatial shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = SS_TS * [1:length(rf)];
    f_g = SS_GAMMA* g *z_shift;
    %rf_shift = rf_shift .* exp(-1i*2*pi* [f_g(2:end),0].*t);

    k_z = SS_GAMMA* (cumsum(g) - sum(g))*SS_TS;
    rf_shift = rf_shift .* exp(-1i*2*pi* k_z*z_shift);

   
