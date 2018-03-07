function res = mtimes(a,bb)

if isa(a,'Wavelet_1d') == 0
    error('In  A.*B only A can be Wavelet operator');
end

ls = length(size(bb));
bb = permute(bb,circshift([1:ls],[0,-a.dim+1]));,
s = size(bb);
b = reshape(bb,size(bb,1),prod(size(bb))/size(bb,1));

res = b*0;

if a.adjoint
    for n=1:size(b,2)
    res(:,n) = IWT_PO(real(b(:,n)),a.wavScale,a.qmf) + i* IWT_PO(imag(b(:,n)),a.wavScale,a.qmf);
  
  end

else
    

  for n=1:size(b,2)
    res(:,n) = FWT_PO(real(b(:,n)),a.wavScale,a.qmf) + i* FWT_PO(imag(b(:,n)),a.wavScale,a.qmf);
  
  end

end
  
res = reshape(res,s);
res = permute(res, circshift([1:ls],[0,a.dim-1]));


