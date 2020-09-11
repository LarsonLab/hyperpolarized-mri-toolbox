function res = mtimes(a,b)
% (c) Ricardo Otazo 2008
if isa(a,'TempFFT') == 0
    error('In  A*B only A can be TempFFT operator');
end
if a.adjoint
    res=ifft(ifftshift(b,a.dim),[],a.dim)*sqrt(size(b,a.dim));
else
    res=fftshift(fft(b,[],a.dim),a.dim)/sqrt(size(b,a.dim)); 
end
