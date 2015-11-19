%
%
% SS_GLOBALS
%
% Common global definitions for SS package
%
% Parameters
%   NAME             DESCRIPTION                                                    DEFAULT
%   ----------       -------------------------------------------------------        ---------
%   'Nucleus'        Nuclei of interest, possible values include                   'Hydrogen'
%                    'Hydrogen', 'Lithium', 'Carbon', 'Sodium', 'Phosphorous'
%   'Max Grad'       Maximum Gradient Strength (G/cm)                               5
%   'Max Slew'       Maximum Slew Rate (G/cm/ms)                                    20
%   'Sample Time'    Sampling Time (s)                                              4e-6
%   'Max B1'         Maximum B1 amplitude (G)                                       0.2
%   'Max Duration'   Maximum Pulse Duration (s)                                     20e-3
%   'Num Lobe Iters' Spectral Sampling Frequency Iterations                         10
%   'Equal Lobes'    Force Equal positive and negative gradient lobes               0
%   'Verse Fraction' Fraction of ramps to use with VERSE                            0.8
%   'Num Fs Test'    Spectral aliasing frequencys to test                           100
%   'Spect Correct'  Spectral Correction with actual sampling                       0
%   'Spect Correct Reg' Regularization factor for spectral correction inversion 	    0
%   'SLR'            Shinnar-Le Roux Correction for large tip pulses                0
%   'Min Order'      Find minimum order FIR filter                                  1
%   'B1 Verse'       Apply B1-VERSE'ing for reduced peak power                      0
%   'Slew Penalty'   In B1-VERSE, use slew-rate penalty to reduce delay sensitivity 0
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
% $Header: /home/adam/cvsroot/src/ss/ss_globals.m,v 1.14 2014/05/22 20:43:59 peder Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Gamma definitions -- all in Hz/G
%
SS_GAMMA_HYDROGEN = 4257.6;  
SS_GAMMA_LITHIUM = 1654.6;
SS_GAMMA_CARBON = 1070.5;
SS_GAMMA_SODIUM = 1126.2;
SS_GAMMA_PHOSPHOROUS = 1723.5;

% Declare globals
%
global SS_INIT_DONE; 

% Nucleus info
%
global SS_NUCLEUS SS_GAMMA; 

% Gradient/timing parameters
%
global SS_MXG SS_MXS SS_TS;

% RF parameters
%
global SS_MAX_B1 SS_MAX_DURATION;

% Design tolerances, parameters
%
global SS_NUM_LOBE_ITERS SS_EQUAL_LOBES SS_VERSE_FRAC SS_NUM_FS_TEST;
global SS_SPECT_CORRECT_FLAG SS_SPECT_CORRECT_REGULARIZATION SS_SLR_FLAG;
global SS_VERSE_B1 SS_SLEW_PENALTY SS_MIN_ORDER;

% Define globals if not already defined
%
if isempty(SS_INIT_DONE), 
    % Indicate init has been done
    %
    SS_INIT_DONE = 1;
    
    % Nucleus info
    %
    SS_NUCLEUS = 'Hydrogen';
    SS_GAMMA = SS_GAMMA_HYDROGEN;

    % Gradient/timing parameters
    %
    SS_MXG = 5.0;			% G/cm
    SS_MXS = 20;			% G/cm/ms 
    SS_TS = 4e-6;			% Sampling time (s) 

    % RF parameters
    %
    SS_MAX_B1 = 0.2; 			% Gauss
    SS_MAX_DURATION = 20e-3;		% Max allowed duration

    
    % Design tolerances, parameters
    %
    SS_NUM_LOBE_ITERS = 10;
    SS_EQUAL_LOBES = 0;
    SS_VERSE_FRAC = 0.8;
    SS_NUM_FS_TEST = 100;
    SS_SPECT_CORRECT_FLAG = 0;
    SS_SPECT_CORRECT_REGULARIZATION = 0;
    SS_SLR_FLAG = 0;
    SS_MIN_ORDER = 1;
    SS_VERSE_B1 = 0;
    SS_SLEW_PENALTY = 0;
end;
