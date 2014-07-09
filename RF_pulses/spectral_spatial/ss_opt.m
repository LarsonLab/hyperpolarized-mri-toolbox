function options = ss_opt(new_options)
% SS_OPT - Set options for spectral-spatial design
%   
% function options = ss_opt(new_options)
%
% options - cell array of current options
% new_options - cell array of new options to set
%
% Pass new_options as [] to reset options
%
% See help ss_globals for options

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
% $Header: /home/adam/cvsroot/src/ss/ss_opt.m,v 1.14 2013/08/15 14:41:09 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ss_globals;

    if (nargin < 1)
	error('Usage: options = ss_opt(new_options)');
    end;
    
    if isempty(new_options)
	SS_INIT_DONE = [];
	ss_globals;
    end;

    new_options = new_options.';
    new_options = new_options(:);
    if (mod(length(new_options), 2) == 1), 
	error('Input parameter (new_options) must be of even length');
    end;

    nopt = length(new_options)/2;
    for idx = 1:nopt, 
	opt = new_options{idx*2-1};
	optparm = new_options{idx*2};
	switch opt,
	 
	 case 'Nucleus'			% Set nucleus and gamma
	  switch optparm,
	   case 'Hydrogen'
	    SS_NUCLEUS = optparm;
	    SS_GAMMA = SS_GAMMA_HYDROGEN;
	   case 'Lithium'
	    SS_NUCLEUS = optparm;;
	    SS_GAMMA = SS_GAMMA_LITHIUM;
	   case 'Carbon'
	    SS_NUCLEUS = optparm;
	    SS_GAMMA = SS_GAMMA_CARBON;
	   case 'Sodium'
	    SS_NUCLEUS = optparm;
	    SS_GAMMA = SS_GAMMA_SODIUM;
	   case 'Phosphorous'
	    SS_NUCLEUS = optparm;
	    SS_GAMMA = SS_GAMMA_PHOSPHOROUS;
	   otherwise
	    error(sprintf('Nucleus: %s not recognized', optparm));
	  end;
	 
	 case 'Max Grad',
	  SS_MXG = optparm;
	 
	 case 'Max Slew', 
	  SS_MXS = optparm;
	  
	 case 'Max B1', 
	  SS_MAX_B1 = optparm;
	  
	 case 'Max Duration', 
	  SS_MAX_DURATION = optparm;
	  
	 case 'Num Lobe Iters', 
	  SS_NUM_LOBE_ITERS = optparm;
	  
	 case 'Equal Lobes', 
	  SS_EQUAL_LOBES = optparm;
	  
	 case 'Verse Fraction', 
	  SS_VERSE_FRAC = optparm;
	  
	 case 'Num Fs Test', 
	  SS_NUM_FS_TEST = optparm;
	  
	 case 'Spect Correct', 
	  SS_SPECT_CORRECT_FLAG = optparm;
	  
	 case 'SLR', 
	  SS_SLR_FLAG = optparm;
	  
	 case 'Min Order', 
	  SS_MIN_ORDER = optparm;
	  
	 case 'B1 Verse', 
	  SS_VERSE_B1 = optparm;
	  
	 case 'Slew Penalty', 
	  SS_SLEW_PENALTY = optparm;

     otherwise
	  SS_INIT_DONE = [];
	  ss_globals;
	  error(sprintf('Option: %s not recognized',opt));
	end;
    end;

    % Fill in options to pass back
    %
    idx = 1;

    options{idx, 1} = 'Nucleus'; 
    options{idx, 2} = SS_NUCLEUS; 
    idx = idx+1;
    
    options{idx, 1} = 'Max Grad';
    options{idx, 2} = SS_MXG;
    idx = idx+1;

    options{idx, 1} = 'Max Slew';
    options{idx, 2} = SS_MXS;
    idx = idx+1;

    options{idx, 1} = 'Sample Time'; 
    options{idx, 2} = SS_TS; 
    idx = idx+1;

    options{idx, 1} = 'Max B1'; 
    options{idx, 2} = SS_MAX_B1; 
    idx = idx+1;

    options{idx, 1} = 'Max Duration'; 
    options{idx, 2} = SS_MAX_DURATION; 
    idx = idx+1;

    options{idx, 1} = 'Num Lobe Iters'; 
    options{idx, 2} = SS_NUM_LOBE_ITERS; 
    idx = idx+1;

    options{idx, 1} = 'Equal Lobes'; 
    options{idx, 2} = SS_EQUAL_LOBES; 
    idx = idx+1;

    options{idx, 1} = 'Verse Fraction'; 
    options{idx, 2} = SS_VERSE_FRAC; 
    idx = idx+1;

    options{idx, 1} = 'Num Fs Test'; 
    options{idx, 2} = SS_NUM_FS_TEST; 
    idx = idx+1;

    options{idx, 1} = 'Spect Correct'; 
    options{idx, 2} = SS_SPECT_CORRECT_FLAG; 
    idx = idx+1;

    options{idx, 1} = 'SLR'; 
    options{idx, 2} = SS_SLR_FLAG; 
    idx = idx+1;

    options{idx, 1} = 'Min Order'; 
    options{idx, 2} = SS_MIN_ORDER; 
    idx = idx+1;

    options{idx, 1} = 'B1 Verse'; 
    options{idx, 2} = SS_VERSE_B1; 
    idx = idx+1;

    options{idx, 1} = 'Slew Penalty'; 
    options{idx, 2} = SS_SLEW_PENALTY; 
    idx = idx+1;


