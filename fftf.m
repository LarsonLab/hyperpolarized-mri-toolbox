function xf = fftf(x, N, dim)
% FFTF - Forward Fourier transform, shifted 
%   
% function xf = fftf(x, N, dim)
%
% x  -- time samples
% N  -- dimension of Fourier space 
%        - length of vector if unspecified
%        - length of column if unspecified and xt a matrix
% dim -- dimension to operate over, only used if matrix
%        - column dimension if unspecified
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
% $Header: /home/adam/cvsroot/src/ss/fftf.m,v 1.3 2012/02/01 00:41:22 peder Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    xsize = size(x);
    if length(xsize) > 2, 
	error('Only handles 2D matrices');
    end;
    
    if nargin < 2, 
	[mnval, mndim] = min(xsize);
	if mnval > 1, 
	    dim = 1;
	    N = xsize(dim);
	else
	    dim = 3-mndim;		% Get opposite of mndim
	    N = xsize(dim);
	end;
    elseif nargin < 3, 
	[mnval, mndim] = min(xsize);
	if mnval > 1, 
	    dim = 1;
	    if (N < xsize(dim))
		error('N less than number of rows of x');
	    end;
	else
	    dim = 3-mndim;		% Get opposite of mndim
	    if (N < xsize(dim))
		error('N less than number of elements of x');
	    end;
	end;
    end;

    % Create zeropad array
    % 
    fsize = xsize;
    fsize(dim) = N;
    pad = zeros(fsize);
    
    % Fill in pad array 
    %
    Nd2 = ceil((N+1)/2);
    nxd2 = ceil((xsize(dim)+1)/2);
    sidx = Nd2 - nxd2 + 1;
    if dim == 1, 
	pad(sidx:sidx+xsize(dim)-1, :) = x;
    else
	pad(:, sidx:sidx+xsize(dim)-1) = x;
    end;
    
    xf = fftshift(fft(ifftshift(pad,dim),N,dim), dim);
    
    return;
    