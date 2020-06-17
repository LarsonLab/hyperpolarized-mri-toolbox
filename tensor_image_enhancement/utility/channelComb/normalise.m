function [val] = normalise(val)
% Normalise a vector.

% Copyright Chris Rodgers, University of Oxford, 2008-11.
% $Id: normalise.m 4188 2011-05-05 15:48:24Z crodgers $

val = val / norm(val);
