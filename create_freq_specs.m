function [fspec, a_angs, d] = create_freq_specs(mets)
% CREATE_FREQ_SPECS - Helper function for creating spectral-spatial
% frequency specifications
%
% [fspec, a_angs, d] = create_freq_specs(mets);
%   
% INPUT
%    mets - structure containing:
%       mets(i).f - center frequency of bands (Hz)
%       mets(i).df - bandwidth of bands (Hz)
%       mets(i).ang - flip angle of bands (degress)
%       mets(i).d - ripple of bands (Mxy, default = .01)
%      
% OUTPUT
%    spec, a_angs, d - inputs for ss_design()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2014 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

for n = 1:length(mets)
    a_angs(n) = mets(n).ang*pi/180;
    fspec(2*n-1) = mets(n).f - mets(n).df;
    fspec(2*n) = mets(n).f + mets(n).df;
    if isfield(mets(n), 'd') && ~isempty(mets(n).d)
        d(n) = mets(n).d;
    else
        d(n) = .01;
    end
end

% by default, center frequency specification 
fmid = (mets(1).f+mets(end).f)/2; 
fspec = fspec - fmid;
