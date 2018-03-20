function res= TempFFT(dim)
% 
% implements a temporal FFT along the dim dimenison
%
% (c) Ricardo Otazo 2008

res.adjoint = 0;
res.dim = dim;
res = class(res,'TempFFT');
