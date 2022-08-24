function res = Wavelet_1d(filterType, filterSize, wavScale, dim)
% res = Wavelet_1d(Filtertype, filterSize, wavScale [, dim])

if nargin < 4
	dim  = 1;
end

res.dim = dim;
res.adjoint = 0;
res.qmf = MakeONFilter(filterType, filterSize);
res.wavScale = wavScale;
res = class(res,'Wavelet_1d');
