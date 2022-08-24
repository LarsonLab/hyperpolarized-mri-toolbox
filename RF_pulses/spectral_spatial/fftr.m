function x = fftr(xf, N, dim)
% FFTR - Reverse Fourier transform, shifted 
%   
% function x = fftr(xf, N, dim)
%
% xf  -- Fourier samples
% N   -- dimension of x space 
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
% $Header: /home/adam/cvsroot/src/ss/fftr.m,v 1.3 2012/02/01 00:41:22 peder Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xfsize = size(xf);
    if length(xfsize) > 2, 
	error('Only handles 2D matrices');
    end;
    
    if nargin < 2, 
	[mnval, mndim] = min(xfsize);
	if mnval > 1, 
	    dim = 1;
	    N = xfsize(dim);
	else
	    dim = 3-mndim;		% Get opposite of mndim
	    N = xfsize(dim);
	end;
    elseif nargin < 3, 
	[mnval, mndim] = min(xfsize);
	if mnval > 1, 
	    dim = 1;
	    if (N > xfsize(dim))
		error('N greater than number of rows of xf');
	    end;
	else
	    dim = 3-mndim;		% Get opposite of mndim
	    if (N > xfsize(dim))
		error('N greater than number of elements of xf');
	    end;
	end;
    end;

    % Get transform along dim
    %
    xpad = fftshift(ifft(ifftshift(xf,dim),[],dim),dim);

    % Fill in pad array 
    %
    Nd2 = ceil((N+1)/2);
    nxfd2 = ceil((xfsize(dim)+1)/2);
    sidx = nxfd2 - Nd2 + 1;
    if dim == 1, 
	x = xpad(sidx:sidx+N-1, :);
    else
	x = xpad(:,sidx:sidx+N-1);
    end;

    return;
    