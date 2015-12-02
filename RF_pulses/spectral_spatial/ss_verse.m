% rfv = ss_verse(gv,rf)
%
% Computes the versed version of rf for a given time-vayring gradient gv

%  written by John Pauly, 1992
%  Bug fixes by Peder Larson, 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)1992-2011 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rfv = ss_verse(g,rf)

[m n] = size(g);
if m<n, g = g.'; end;
[m n] = size(rf);
if m<n,
  rf = conj(rf');
  [m n] = size(rf);
end;

% reference to middle of gradients
k = cumsum([0;g(1:end-1)]) + g/2;
k = m*k/max(k);

rfv = [];
g = m*g/sum(g);

for j=1:n,
  % interpolate at half integer values
  rft = g.*interp1([1:m]-0.5,rf(:,j),k,'spline','extrap');
  %  rft = g.*interp1([1:m]-0.5,rf(:,j),k,'linear','extrap');
  rfv = [rfv rft];
end;
